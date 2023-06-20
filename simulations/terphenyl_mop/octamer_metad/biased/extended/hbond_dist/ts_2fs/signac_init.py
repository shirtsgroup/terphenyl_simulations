import signac
import numpy as np
import os
import shutil
import re
import itertools

def replace_all_pattern(pattern, replace, file):
    with open(file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = re.sub(pattern, replace, lines[i])
    with open(file, "w") as f:
        f.writelines(lines)


def main():
    n_walkers = 16
    project = signac.get_project()
    # combinations
    statepoints = [
        {'height' : 2.5,   'sigma' : 0.5, 'bf' : 50},
        {'height' : 2.5,   'sigma' : 0.5, 'bf' : 100},
        {'height' : 2.5,   'sigma' : 0.5, 'bf' : 300},
        {'height' : 2.5,   'sigma' : 0.5, 'bf' : 500},
        {'height' : 0.001, 'sigma' : 0.5, 'bf' : 100000},
        {'height' : 0.003, 'sigma' : 0.5, 'bf' : 100000},
        {'height' : 0.005, 'sigma' : 0.5, 'bf' : 100000},
    ]
    replicas = 1

    for statepoint in statepoints:
        for r in range(replicas):
            statepoint['replica'] = r
            sp_dict = statepoint
            if not len(project.find_jobs(sp_dict).to_dataframe()) == 0:
                print("This statepoint is already defined. Skipping.", sp_dict)
                continue

            job = project.open_job(sp_dict)
            # setup input files for metadynamics simulation
            shutil.copytree("simulation_template", job.path)
            walker_dirs = " ".join(["WALKER" + str(walker_id) for walker_id in range(n_walkers)])
            replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn("submit.berendsen_npt.slurm"))
            replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn("submit.berendsen_nvt.slurm"))
            replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn("submit.production.slurm"))
            replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn("submit.production_finish.slurm"))
            replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn("submit.continue.slurm"))
            for walker_id in range(n_walkers):
                walker_dir = os.path.join(job.path, "WALKER" + str(walker_id))
                shutil.copytree(os.path.join(job.path, "WALKER"), walker_dir)
                replace_all_pattern("height", str(statepoint['height']), os.path.join(walker_dir, "plumed_hbond_dist.dat"))
                replace_all_pattern("sigma", str(statepoint['sigma']), os.path.join(walker_dir, "plumed_hbond_dist.dat"))
                replace_all_pattern("bf", str(statepoint['bf']), os.path.join(walker_dir, "plumed_hbond_dist.dat"))
                replace_all_pattern("bf", str(statepoint['bf']), os.path.join(walker_dir,  "plumed_reweight.dat"))
                replace_all_pattern("sigma", str(statepoint['sigma']), os.path.join(walker_dir, "plumed_reweight.dat"))
                replace_all_pattern("TEMP", str(300), os.path.join(walker_dir, "berendsen_nvt.mdp"))
                replace_all_pattern("TEMP", str(300), os.path.join(walker_dir, "berendsen_npt.mdp"))
                replace_all_pattern("TEMP", str(300), os.path.join(walker_dir, "npt_new.mdp"))
                replace_all_pattern("temp", str(300), os.path.join(walker_dir, "plumed_hbond_dist.dat"))
                replace_all_pattern("temp", str(300), os.path.join(walker_dir, "plumed_reweight.dat"))
            shutil.rmtree(os.path.join(job.path, "WALKER"))
            job.init()

if __name__ == "__main__":
    main()

