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
import MDAnalysis as md
from tqdm import tqdm
from flow import FlowProject
from MDAnalysis import Universe
import terphenyl_simulations
from natsort import natsorted
from tqdm import tqdm
import mdtraj
from terphenyl_simulations.utils import replace_all_pattern, get_solvent_structure_file
from terphenyl_simulations.analysis_workflows.labels import *
from terphenyl_simulations.edit_conf import InternalCoordinateEditor
import warnings
warnings.filterwarnings("ignore")

# Initialize Signac Project

def signac_init():
    # Open parameter file and create a list of statepoints to define
    simulation_statepoints = []
    with open("md_parameters.yml", "r") as f:
        simulation_parameters = yaml.safe_load(f)

    # Remove replicas and replace with replica_id
    for i in range(simulation_parameters["n_simulations"]):
        for temp in simulation_parameters["temperatures"]:
            # Replica index
            sp_i = dict(simulation_parameters)
            sp_i["replica"] = i
            del sp_i["n_simulations"]

            # Simulation temperatures
            sp_i["temperature"] = temp
            del sp_i["temperatures"]
            simulation_statepoints.append(sp_i)

    project = signac.get_project()

    # Setup simulation directories
    for sp in simulation_statepoints:
        # Quick check to see if sp exists already
        job = project.open_job(sp)
        if "init" in job.doc.keys():
            continue
        
        job.doc["init"] = True

        # Default file templates
        if not "system" in job.sp.keys():
            job.sp["system"] = "cu_alpine"

        # Setup job directory with template files
        md_files = glob.glob(
            os.path.join(
                terphenyl_simulations.utils.ROOT_DIR,
                "data/simulation_templates",
                "md/",
                job.sp["system"] + "/*"                
            )
        )

        for sim_file in md_files:
            shutil.copy(sim_file, job.path)
            if ".mdp" in sim_file:
                replace_all_pattern(
                    "TEMP",
                    str(job.sp["temperature"]),
                    job.fn(sim_file.split("/")[-1])
                )
        shutil.copy(
            simulation_parameters["build_foldamer"],
            job.fn(simulation_parameters["build_foldamer"]),
        )
        shutil.copy("md_parameters.yml", job.fn("md_parameters.yml"))

        if os.path.exists(simulation_parameters["helix_torsions"]):
            shutil.copy(
                simulation_parameters["helix_torsions"],
                job.fn("helix_torsions.yml"),
            )

        if os.path.exists("torsions.yml"):
            shutil.copy(
                "torsions.yml",
                job.fn("torsions.yml"),
            )

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
    tm = terphenyl_simulations.topology_manager.TopologyManager()
    # print(tm.topology_dictionary)
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
    lambda job: os.path.exists(job.fn("em_" + job.doc["foldamer_name"] + ".tpr"))
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
    gmx_wrapper.trjconv(
        f="em_" + job.doc["foldamer_name"] + ".gro",
        s = "em_" + job.doc["foldamer_name"] + ".tpr",
        o="em_" + job.doc["foldamer_name"] + ".pdb",
        conect="yes",
    )
    job.doc["foldamer_gro"] = "em_" + job.doc["foldamer_name"] + ".gro"
    job.doc["foldamer_pdb"] = "em_" + job.doc["foldamer_name"] + ".pdb"

@FlowProject.pre.after(minimize_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(job.fn("helix_" + job.doc["foldamer_name"] + ".gro"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def set_helix_torsions(job):
    internal_coor_editor = InternalCoordinateEditor(
        job.doc["foldamer_gro"], "em_" + job.doc["foldamer_name"] + ".tpr",
    )

    # Torsion definitions for first residue
    with open("torsions.yml", "r") as yml_read:
        monomer_torsions = yaml.safe_load(yml_read)

    with open("helix_torsions.yml", "r") as stream:
        helix_torsions = dict(yaml.safe_load(stream))


    for torsion_type in monomer_torsions['torsions'].keys():
        print("Adjusting Torsion Type:", torsion_type)
        torsion_atom_ids = terphenyl_simulations.utils.get_torsion_atom_ids(
            monomer_torsions['torsions'][torsion_type],
            monomer_torsions['offset'],
            job.doc["build_parameters"]["foldamer_length"],
        )

        ag = internal_coor_editor.universe.select_atoms("all")
        for torsion_id in torsion_atom_ids:
            torsion_atom_names = [ag.atoms[atom_id].name for atom_id in torsion_id]
            torsion_atom_ids, torsion_values = internal_coor_editor.find_torsions(torsion_atom_names[1:3], positions = [1, 2])
            internal_coor_editor.set_torsion(torsion_atom_names, helix_torsions[torsion_type] * np.pi / 180)
            reversed_torsion = (((helix_torsions[torsion_type] + 360) % 360) - 180)
            reversed_torsion =  helix_torsions[torsion_type]
            internal_coor_editor.set_torsion(torsion_atom_names[::-1], reversed_torsion * np.pi / 180)
            internal_coor_editor.update_internal_coordinates()
    
    internal_coor_editor.write_structure("helix_" + job.doc["foldamer_name"] + ".gro")

    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )

    gmx_wrapper.center_configuration(
        "helix_" + job.doc["foldamer_name"] + ".gro",
        "helix_" + job.doc["foldamer_name"] + "_centered.gro",
    )

    job.doc["foldamer_gro"] = "helix_" + job.doc["foldamer_name"] + "_centered.gro"

    
@FlowProject.pre.after(set_helix_torsions)
@FlowProject.post(
    lambda job: os.path.exists(job.fn("em_helix_" + job.doc["foldamer_name"] + ".gro"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def minimize_helix_foldamer(job):
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )
    centered_out_name = job.doc["foldamer_gro"].split(".gro")[0] + "_centered.gro"
    gmx_wrapper.center_configuration(job.doc["foldamer_gro"], centered_out_name)
    gmx_wrapper.minimize(
        centered_out_name,
        job.doc["foldamer_topology"],
        prefix="em_helix_" + job.doc["foldamer_name"],
    )
    gmx_wrapper.trjconv(
        f="em_helix_" + job.doc["foldamer_name"] + ".gro",
        s="em_helix_" + job.doc["foldamer_name"] + ".tpr",
        o="em_helix_" + job.doc["foldamer_name"] + ".pdb",
        conect="yes",
    )
    job.doc["foldamer_gro"] = "em_helix_" + job.doc["foldamer_name"] + ".gro"
    job.doc["foldamer_pdb"] = "em_helix_" + job.doc["foldamer_name"] + ".pdb"



@FlowProject.pre.after(minimize_helix_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(job.fn(job.doc["foldamer_name"] + "_" + job.doc["system_name"] + "_" + job.doc["build_parameters"]["ff_method"] + ".top"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def build_system(job):
    print(job.fn(job.doc["foldamer_name"] + "_" + job.doc["system_name"] + "_" + job.doc["build_parameters"]["ff_method"] + ".top"))
    ff_names = ["openff-2.0.0", "openff-1.0.0"]
    if job.doc["build_parameters"]["ff_method"] == "bespoke":
        ff_names = glob.glob("*bespoke*.offxml") + ["openff-1.0.0"]
        ff_names = [f.split(".offxml")[0] for f in ff_names]
    openff_builder = terphenyl_simulations.build.SystemBuilderOpenFF(
        job.doc["foldamer_pdb"],
        job.doc["build_parameters"]["system"]["solvent"],
        job.sp["build_foldamer"],
    )
    openff_builder.get_system_structure()

    job.doc["system_gro"] = openff_builder.gro_file
    job.doc["system_pdb"] = openff_builder.pdb_file

    openff_topology_gen = terphenyl_simulations.build.SystemTopologyGenerator(
        job.doc["system_pdb"],
        [job.doc["foldamer_pdb"], get_solvent_structure_file(job.doc["build_parameters"]["system"]["solvent"])],
        [job.doc["foldamer_name"] + "_charges.sdf", None],
        job.doc["foldamer_name"] + "_" + job.doc["system_name"] + "_" + job.doc["build_parameters"]["ff_method"],
        job.doc["build_parameters"]["ff_method"],
        job.sp["build_foldamer"],
        ff_names = ff_names,
        topology_manager = openff_builder.topology_manager
    )
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )
    openff_topology_gen.set_simulation_engine(gmx_wrapper)
    openff_topology_gen.get_parameters()

    job.doc["system_topology"] = openff_topology_gen.top_file
    job.doc["system_gro"] = openff_topology_gen.gro_file


@FlowProject.pre.after(build_system)
@FlowProject.post(
    lambda job: os.path.exists(job.fn(job.doc["foldamer_name"] + "_" + job.doc["system_name"] + "_" + job.doc["build_parameters"]["ff_method"] + "_hmr.top"))
)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def apply_hmr_to_topology(job):
    output_topology = job.doc["system_topology"].split(".top")[0] + "_hmr.top"
    sys.argv = ["hmr_topology", "-t", job.doc["system_topology"], "-o",  output_topology, "--hmr_ratio", "3"]
    terphenyl_simulations.scripts.hmr_topology()
    tm = terphenyl_simulations.build.TopologyManager()
    tm.add_topology(output_topology, job.sp["build_foldamer"], "system")
    job.doc["system_topology"] = output_topology


@FlowProject.pre.after(apply_hmr_to_topology)
@FlowProject.post(lambda job: os.path.exists(job.fn(job.doc["foldamer_name"] + "_" + job.doc["system_name"] + "_" + job.doc["build_parameters"]["ff_method"] + "_hmr_rest.top")))
@FlowProject.post(lambda job: os.path.exists(job.fn("posre.itp")))
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def apply_position_restraints(job):
    output_topology = job.doc["system_topology"].split(".top")[0] + "_rest.top"
    # generate position restraints itp file
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )

    gmx_wrapper.gmx_command(
        "make_ndx",
        {"f" : job.doc["system_gro"]},
        inputs = '2|3\nq\n'
    )

    gmx_wrapper.gmx_command(
        "genrestr",
        {"f" : job.doc["system_gro"], "n" : "index.ndx"},
        inputs = '5\n'
    )

    # Write topology file with posres itp included
    molecule_index = 0
    with open(job.doc["system_topology"], "r") as top_read:
        with open(output_topology, "w") as top_write:
            for line in top_read.readlines():
                if "[ moleculetype ]" in line:
                    molecule_index += 1
                if "[ moleculetype ]" in line and molecule_index == 2:
                    top_write.write('#include "posre.itp"\n\n')
                top_write.write(line)

    job.doc["system_restrained_topology"] = output_topology


@FlowProject.pre.after(apply_position_restraints)
@FlowProject.post(lambda job: "slurm_modified" in job.doc.keys())
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def modify_slurm_files(job):
    for submit_file in glob.glob("submit.berendsen*.slurm"):
        replace_all_pattern("TOPOLOGY_FILE", job.doc["system_restrained_topology"], submit_file)
        replace_all_pattern("INITIAL_STRUCTURE", job.doc["system_gro"], submit_file)
        replace_all_pattern("SIMULATION_NAME", job.doc["build_parameters"]["structure_file"], submit_file)

    for submit_file in ["submit.production.slurm", "submit.finish.slurm", "submit.continue.slurm" ]:
        replace_all_pattern("TOPOLOGY_FILE", job.doc["system_topology"], submit_file)
        replace_all_pattern("INITIAL_STRUCTURE", job.doc["system_gro"], submit_file)
        replace_all_pattern("SIMULATION_NAME", job.doc["build_parameters"]["structure_file"], submit_file)
    job.doc["slurm_modified"] = True


@FlowProject.pre(lambda job: os.path.exists(job.fn("production_npt.gro")))
@FlowProject.pre(lambda job: os.path.exists(job.fn("production_npt.xtc")))
@FlowProject.pre(lambda job: shutil.which("gmx"))
@FlowProject.post(lambda job: os.path.exists(job.fn("production_npt.whole.xtc")))
@FlowProject.pre.after(modify_slurm_files)
@FlowProject.operation(directives={"fork": True})
@cd_to_job_dir
def apply_pbcs(job):
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )
    for simulation_file in glob.glob("*.xtc"):
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



def main():
    if not os.path.isdir("workspace"):
        subprocess.run("signac init".split(" "))
        signac_init()
    FlowProject().main()


if __name__ == "__main__":
    main()