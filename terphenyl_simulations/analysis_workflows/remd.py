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
        job.sp["build_foldamer"],
        path=job.fn("")
    )
    foldamer_builder.build_foldamer()
    foldamer_builder.write_pdb()
    foldamer_builder.write_mol()

@FlowProject.pre.after(build_foldamer)
@FlowProject.post(
        lambda job: os.path.exists(
            job.fn("solvated_" + job.doc["build_parameters"]["structure_file"] + ".pdb")
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
    packmol_builder.solvate_system()
    os.chdir(top_dir)

@FlowProject.operation
def parameterize_solvated_system(job):
    pass


@FlowProject.pre.after(build_foldamer)
@FlowProject.post(
    lambda job: os.path.exists(
        job.fn(job.doc["build_parameters"]["structure_file"] + "_openff-2.0.0.top")
    )
)
@FlowProject.operation
def parameterize_foldamer(job):
    top_generator = terphenyl_simulations.build.SimulationTopologyGenerator(
        job.sp["build_foldamer"],
        job.doc["build_parameters"]["ff_method"],
        path = job.fn("")
    )
    top_generator.assign_parameters()
    job.doc['foldamer_topology'] = top_generator.top_file
    job.doc['foldamer_gro'] = top_generator.gro_file

@FlowProject.pre.after(parameterize_foldamer)
@FlowProject.operation
def minimize_foldamer(job):
    top_dir = os.path.abspath(".")
    os.chdir(job.fn(""))
    gmx_wrapper = terphenyl_simulations.gromacs_wrapper.GromacsWrapper(job.sp['gromacs_exe'])
    out_name = job.doc['foldamer_gro'].split('.gro')[0] + '_centered.gro'
    gmx_wrapper.center_configuration(job.doc['foldamer_gro'], out_name)
    job.doc['foldamer_gro'] = out_name
    gmx_wrapper.minimize(job.doc['foldamer_gro'], job.doc['foldamer_topology'])
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
