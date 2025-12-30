import shutil
import functools
import os
import subprocess
import yaml
import glob
import signac
import flow
import sys
import numpy as np
from flow import FlowProject
from MDAnalysis import Universe
import terphenyl_simulations
from natsort import natsorted
from tqdm import tqdm
import mdtraj
from terphenyl_simulations.utils import replace_all_pattern

# Initialize Signac Project


def signac_init():
    # Open parameter file and create a list of statepoints to define
    simulation_statepoints = []
    with open("remd_parameters.yml", "r") as f:
        simulation_parameters = yaml.safe_load(f)

    # Remove replicas and replace with replica_id
    for i in range(simulation_parameters["n_simulations"]):
        sp_i = dict(simulation_parameters)
        sp_i["replica"] = i
        del sp_i["n_simulations"]
        simulation_statepoints.append(sp_i)

    project = signac.get_project()

    # Setup simulation directories
    for sp in simulation_statepoints:
        # Quick check to see if sp exists already
        job = project.open_job(sp)
        if "init" in job.doc.keys():
            continue
        job.doc["init"] = True

        # Setup job directory with template files
        remd_files = glob.glob(
            os.path.join(
                terphenyl_simulations.utils.ROOT_DIR,
                "data/simulation_templates",
                "remd/*",
            )
        )

        for sim_file in remd_files:
            shutil.copy(sim_file, job.path)
        shutil.copy(
            simulation_parameters["build_foldamer"],
            job.fn(simulation_parameters["build_foldamer"]),
        )
        shutil.copy("remd_parameters.yml", job.fn("remd_parameters.yml"))

        with open(job.fn(simulation_parameters["build_foldamer"]), "r") as f:
            job.doc["build_parameters"] = yaml.safe_load(f)
        job.doc["foldamer_name"] = job.doc["build_parameters"]["structure_file"]
        job.doc["system_name"] = "system"


# Decorator to cd into and out of workspace
# before and after operation
def cd_to_job_dir(function):
    @functools.wraps(function)
    def wrap_flow_operation(job):
        top_dir = os.path.abspath(".")
        os.chdir(job.fn(""))
        function(job)
        os.chdir(top_dir)

    return wrap_flow_operation


# FlowProject Operations
@FlowProject.post(lambda job: os.path.exists(job.fn(job.doc["foldamer_name"] + ".pdb")))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def build_foldamer(job):
    foldamer_builder = terphenyl_simulations.build.FoldamerBuilder(
        job.sp["build_foldamer"]
    )
    foldamer_builder.get_foldamer()


@FlowProject.pre.after(build_foldamer)
@FlowProject.post(
    lambda job: glob.glob(job.fn(job.doc["foldamer_name"] + "*.top"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def parameterize_foldamer(job):
    mol_file = job.doc["foldamer_name"] + ".mol"
    pdb_file = job.doc["foldamer_name"] + ".pdb"
    top_generator = terphenyl_simulations.build.MoleculeTopologyGenerator(
        mol_file,
        pdb_file,
        job.sp["build_foldamer"],
        job.doc["build_parameters"]["ff_method"],
    )
    top_generator.get_ff_parameters()
    job.doc["foldamer_topology"] = top_generator.top_file
    job.doc["foldamer_gro"] = top_generator.gro_file


@FlowProject.pre.after(parameterize_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(job.fn("em_" + job.doc["foldamer_name"] + ".pdb"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def minimize_foldamer(job):
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )
    centered_out_name = job.doc["foldamer_gro"].split(".gro")[0] + "_centered.gro"
    gmx_wrapper.center_configuration(job.doc["foldamer_gro"], centered_out_name)
    gmx_wrapper.minimize(
        centered_out_name,
        job.doc["foldamer_topology"],
        prefix="em_" + job.doc["foldamer_name"],
    )
    gmx_wrapper.edit_conf(
        f="em_" + job.doc["foldamer_name"] + ".tpr",
        o="em_" + job.doc["foldamer_name"] + ".pdb",
        conect="yes",
    )
    job.doc["foldamer_gro"] = "em_" + job.doc["foldamer_name"] + ".gro"


@FlowProject.pre.after(minimize_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(job.fn(job.doc["foldamer_name"] + "_system.top"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def build_system(job):
    openff_builder = terphenyl_simulations.build.SystemBuilderOpenFF(
        job.doc["foldamer_name"] + "_charges.sdf",
        job.doc["build_parameters"]["system"]["solvent_smile_str"],
        job.sp["build_foldamer"],
    )
    openff_builder.get_system_topology()
    openff_builder.minimize_system()
    job.doc["foldamer_topology"] = openff_builder.top_file
    job.doc["foldamer_gro"] = openff_builder.gro_file

@FlowProject.pre.after(build_system)
@FlowProject.post(
    lambda job: os.path.exists(job.fn(job.doc["foldamer_name"] + "_system_hmr.top"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def apply_hmr_to_topology(job):
    output_topology = job.doc["foldamer_name"] + "_system_hmr.top"
    sys.argv = ["hmr_topology", "-t", job.doc["foldamer_topology"], "-o",  output_topology, "--hmr_ratio", "3"]
    terphenyl_simulations.scripts.hmr_topology()
    tm = terphenyl_simulations.build.TopologyManager()
    tm.add_topology(output_topology, job.sp["build_foldamer"], "system")
    job.doc["foldamer_topology"] = output_topology
    for submit_file in glob.glob("submit*.slurm"):
        replace_all_pattern("TOPOLOGY_FILE", output_topology, submit_file)


@FlowProject.pre.after(apply_hmr_to_topology)
@FlowProject.post(lambda job: os.path.isdir(job.fn("sim0")))
@FlowProject.operation(directives={"fork" : True})
@cd_to_job_dir
def setup_remd_simulations(job):
    # Setup input arguments
    sys.argv = ["REMD_setup", "-N", str(job.sp["n_replicas"]),
                "--t_range", str(job.sp["t_range"][0]), str(job.sp["t_range"][1]),
                "--sim_id", job.sp["sim_id"],
                "--topology_files", job.doc["foldamer_gro"], job.doc["foldamer_topology"]
            ]

    sys.argv += ["--mdps"] + list(job.sp["mdps"])
    # Run REMD setup
    terphenyl_simulations.scripts.REMD_setup()

# if slurm is an executable
@FlowProject.pre.after(setup_remd_simulations)
@FlowProject.pre(lambda job: shutil.which("sbatch"))
@FlowProject.post(lambda job: os.path.exists(job.fn("sim0/production_npt.log")))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def submit_simulations(job):
    subprocess.Popen(["bash", "submit_all.slurm"], shell = True)
    subprocess.wait()

# if slurm isn't an executable
@FlowProject.pre.after(setup_remd_simulations)
@FlowProject.pre(lambda job: shutil.which("gmx_mpi"))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def run_simulations(job):
    pass


@FlowProject.pre(lambda job: os.path.exists(job.fn("sim0/production_npt.gro")))
@FlowProject.pre(lambda job: os.path.exists(job.fn("sim0/production_npt.xtc")))
@FlowProject.pre(lambda job: shutil.which("gmx"))
@FlowProject.post(lambda job: os.path.exists(job.fn("sim0/production_npt.whole.xtc")))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def apply_pbcs(job):
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )
    for replica_dir in natsorted(glob.glob(job.sp["sim_id"] + "*")):
        for simulation_file in glob.glob(os.path.join(replica_dir, "*.xtc")):
            if "whole.xtc" in simulation_file:
                continue
            output_filename = simulation_file.split(".xtc")[0] + ".whole.xtc"
            tpr_file = simulation_file.split(".xtc")[0] + ".tpr"
            gmx_wrapper.trjconv(
                inputs = [0],
                hide_outputs = False,
                f = simulation_file,
                o = output_filename,
                s = tpr_file,
                pbc = "whole"
            )

@FlowProject.pre(lambda job: os.path.exists(job.fn("sim0/production_npt.whole.xtc")))
@FlowProject.pre(lambda job: shutil.which("gmx"))
@FlowProject.post(lambda job: os.path.isdir(job.fn("clustering_output")))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def cluster_trajectory(job):
    n_lowest_replicas = 10
    production_sim = job.sp["mdps"][-1].split(".")[0]
    top_file = os.path.join(job.sp["sim_id"] + "0", production_sim + ".gro")
    simulation_trajectories = \
        natsorted(glob.glob(os.path.join(job.sp["sim_id"] + "*", production_sim + ".whole.xtc")))

    select_string = "not resname MOL1"
    print(select_string)
    terphenyl_simulations.clustering.clustering_grid_search(
        simulation_trajectories[:n_lowest_replicas],
        top_file,
        select_string,
        n_min_samples=40,
        n_eps=40,
        n_processes=32,
        prefix="grid_search",
        eps_limits=[0.01, 0.2],
        min_sample_limits=[0.005, 0.1],
        plot_filename="ss.png",
        frame_stride=1,
    )

@FlowProject.pre(lambda job: os.path.exists(job.fn("sim0/production_npt.whole.xtc")))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def plot_remd_torsion_distributions(job):

    # Load REMD simulations to file
    production_sim = job.sp["mdps"][-1].split(".")[0]
    structure_universe = Universe(job.sp["sim_id"] + "0/" + production_sim + ".tpr", job.sp["sim_id"] + "0/" + production_sim + ".gro")
    remd_file_list = [job.sp["sim_id"] + str(i) + "/" + production_sim + ".whole.xtc" for i in range(job.sp["n_replicas"])]
    print("Loading REMD trajectory files...")
    remd_trajs = [
        mdtraj.load(xtc_file, top="sim0/berendsen_npt.gro")
        for xtc_file in tqdm(remd_file_list)
    ]
    
    
    # Torsion definitions for first residue
    with open("torsions.yml", "r") as stream:
        monomer_torsions = yaml.safe_load(stream)


    # Torsion Analysis
    output_dir = "torsion_plots"
    terphenyl_simulations.utils.make_path(output_dir)
    for torsion_type in monomer_torsions["torsions"].keys():
        print("Working on", torsion_type, "torsion...")
        torsion_atom_ids = terphenyl_simulations.utils.get_torsion_atom_ids(
            monomer_torsions["torsions"][torsion_type],
            monomer_torsions["offset"],
            job.doc["build_parameters"]["foldamer_length"],
        )

        terphenyl_simulations.plotting.plot_torsions_distributions(
            remd_trajs,
            torsion_atom_ids,
            torsion_type + "Torsion (radians)",
            os.path.join(output_dir, torsion_type + "_remd"),
            torsion_type + " Torsion Plot",
            figsize=[5, 5],
            cbar_params=[250, 450, "Temperature (K)"],
        )


def main():
    if not os.path.isdir("workspace"):
        subprocess.run("signac init".split(" "))
        signac_init()
    FlowProject().main()


if __name__ == "__main__":
    main()
