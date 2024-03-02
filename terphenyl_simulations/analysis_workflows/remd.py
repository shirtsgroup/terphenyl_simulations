import shutil
import functools
import os
import subprocess
import yaml
import glob
import signac
from flow import FlowProject
import terphenyl_simulations

# Initialize Signac Project


def signac_init():
    # Open parameter file and create a list of statepoints to define
    simulation_statepoints = []
    with open("remd_parameters.yml", "r") as f:
        simulation_parameters = yaml.safe_load(f)

    # Remove replicas and replace with replica_id
    for i in range(simulation_parameters["n_replicas"]):
        sp_i = dict(simulation_parameters)
        sp_i["replica"] = i
        del sp_i["n_replicas"]
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

def main():
    if not os.path.isdir("workspace"):
        subprocess.run("signac init".split(" "))
        signac_init()
    FlowProject().main()

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
@FlowProject.post(
    lambda job: os.path.exists(
        job.fn(job.doc["foldamer_name"] + ".pdb")
    )
)
@FlowProject.operation
@cd_to_job_dir
def build_foldamer(job):
    foldamer_builder = terphenyl_simulations.build.FoldamerBuilder(
        job.sp["build_foldamer"]
    )
    foldamer_builder.build_foldamer()
    foldamer_builder.write_pdb()
    foldamer_builder.write_mol()


@FlowProject.pre.after(build_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(
        job.fn(job.doc["foldamer_name"] + "_openff-2.0.0.top")
    )
)
@FlowProject.operation
@cd_to_job_dir
def parameterize_foldamer(job):
    mol_file = job.doc["foldamer_name"] + ".mol"
    pdb_file = job.doc["foldamer_name"] + ".pdb"
    top_generator = terphenyl_simulations.build.TopologyGenerator(
        mol_file,
        pdb_file,
        None,
        job.doc["build_parameters"]["ff_method"],
    )
    top_generator.assign_parameters()
    job.doc["foldamer_topology"] = top_generator.top_file
    job.doc["foldamer_gro"] = top_generator.gro_file


@FlowProject.pre.after(parameterize_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(
        job.fn("em.gro")
    )
)
@FlowProject.operation
@cd_to_job_dir
def minimize_foldamer(job):
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(
        job.sp["gromacs_exe"]
    )
    out_name = job.doc["foldamer_gro"].split(".gro")[0] + "_centered.gro"
    gmx_wrapper.center_configuration(job.doc["foldamer_gro"], out_name)
    job.doc["foldamer_gro"] = out_name
    gmx_wrapper.minimize(job.doc["foldamer_gro"], job.doc["foldamer_topology"], prefix = "em_" + job.doc["str"])

@FlowProject.pre.after(minimize_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(
        job.fn("solvated_" + job.doc["foldamer_name"] + ".pdb")
    )
)
@FlowProject.operation
@cd_to_job_dir
def build_system(job):
    packmol_builder = terphenyl_simulations.build.SystemBuilder(
        job.sp["build_foldamer"]
    )
    packmol_builder.build_packmol_inp()
    packmol_builder.solvate_system()

@FlowProject.pre.after(build_system)
@FlowProject.post(
    lambda job: os.path.exists(
        job.fn("system_openff.top")
    )
)
@FlowProject.operation
@cd_to_job_dir
def parameterize_solvated_system(job):
    mol_files = [ job.doc["foldamer_name"] + "_charges.sdf", job.doc["build_parameters"]["system"]["solvent"] + ".pdb"]
    pdb_file = "solvated_" + job.doc["foldamer_name"] + ".pdb"
    top_generator = terphenyl_simulations.build.TopologyGenerator(
        mol_files,
        pdb_file,
        "system_openff",
        "openff-system",
    )
    top_generator.assign_parameters()
    job.doc["foldamer_topology"] = top_generator.top_file
    job.doc["foldamer_gro"] = top_generator.gro_file

if __name__ == "__main__":
    main()
