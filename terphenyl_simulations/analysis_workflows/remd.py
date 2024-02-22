import shutil
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


# FlowProject Operations


@FlowProject.post(
    lambda job: os.path.exists(
        job.fn(job.doc["build_parameters"]["structure_file"] + ".pdb")
    )
)
@FlowProject.operation
def build_foldamer(job):
    foldamer_builder = terphenyl_simulations.build.FoldamerBuilder(
        job.sp["build_foldamer"]
    )
    foldamer_builder.build_foldamer(path=job.fn(""))


@FlowProject.pre.after(build_foldamer)
@FlowProject.operation
def parameterize_foldamer(job):
    pass


@FlowProject.post(
    (
        lambda job: os.path.exists(
            job.fn("solvated_" + job.doc["build_parameters"]["structure_file"] + ".pdb")
        )
    )
)
@FlowProject.operation
def build_system(job):
    top_dir = os.path.abspath(".")
    os.chdir(job.fn(""))
    packmol_builder = terphenyl_simulations.build.SystemBuilder(
        job.sp["build_foldamer"]
    )
    packmol_builder.build_packmol_inp()
    packmol_builder.build_system()
    os.chdir(top_dir)


@FlowProject.pre.after(build_system)
@FlowProject.operation
def parameterize_system(job):
    pass


def main():
    if not os.path.isdir("workspace"):
        subprocess.run("signac init".split(" "))
        signac_init()
    FlowProject().main()


if __name__ == "__main__":
    main()
