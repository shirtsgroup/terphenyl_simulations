import shutil
import os
import subprocess
import yaml
import signac
from flow import FlowProject
import terphenyl_simulations

def signac_init():
    # Open parameter file and create a list of statepoints to define
    simulation_statepoints = []
    with open('remd_parameters.yml', 'r') as f:
        simulation_parameters = yaml.safe_load(f)

    # Remove replicas and replace with replica_id
    for i in range(simulation_parameters['n_replicas']):
        sp_i = dict(simulation_parameters)
        sp_i['replica'] = i
        del sp_i['n_replicas']
        simulation_statepoints.append(sp_i)

    project = signac.get_project()

    # Setup simulation directories
    for sp in simulation_statepoints:
        # Quick check to see if sp exists already
        job = project.open_job(sp)
        if "init" in job.doc.keys():
            continue
        job['init'] = True

        # Setup job directory with template files
        shutil.copytree(os.path.join(terphenyl_simulations.utils.ROOT_DIR, 'simulation_templates', 'remd'), job.path)
        foldamer_builder = terphenyl_simulations.build.FoldamerBuilder('mop_tetramer.build')
        foldamer_builder.build_foldamer(path = job.fn(""))

if __name__ == '__main__':
    if not os.path.isdir("workspace"):
        subprocess.run("signac init".split(" "))
        terphenyl_simulations.analysis_workflows.remd.signac_init()

if __name__ == "__main__":
    if not os.path.isdir("workspace"):
        subprocess.run("signac init".split(" "))
    signac_init()
    FlowProject().main()
        