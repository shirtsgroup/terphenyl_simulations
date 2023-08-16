import signac
import numpy as np
import os
import shutil
import re
import glob
import itertools

def replace_all_pattern(pattern, replace, file):
    with open(file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = re.sub(pattern, replace, lines[i])
    with open(file, "w") as f:
        f.writelines(lines)


def main():
    n_walkers = 4
    project = signac.get_project()
    # combinations
    statepoints = [
        {
            'height1' : 0,    'height2' : 2.5, 'height3' : 2.5, 'height4' : 2.5,
            'sigma1'  : 0,    'sigma2'  : 1.0, 'sigma3'  : 0.1, 'sigma4'  : 0.4,
            'bf1'     : 10e9, 'bf2'     : 500, 'bf3'     : 500, 'bf4'     : 500,
         }
    ]
    replicas = 3

    for statepoint in statepoints:
        for r in range(replicas):
            statepoint['simulation_rep'] = r
            sp_dict = statepoint

            job = project.open_job(sp_dict)
            # setup input files for metadynamics simulation
            shutil.copytree("simulation_template", job.path)
            walker_dirs = " ".join(["WALKER" + str(walker_id) for walker_id in range(n_walkers)])
            
            # Replace WALKER_DIRS in all slurm files
            for slurm_file in glob.glob(job.fn("submit*.slurm")):
                replace_all_pattern("WALKER_DIRS", walker_dirs, slurm_file)

            # Make walker directories
            for walker_id in range(n_walkers):
                walker_dir = os.path.join(job.path, "WALKER" + str(walker_id))
                shutil.copytree(os.path.join(job.path, "WALKER"), walker_dir)

                # Replace mdp options in WALKER files
                for mdp_file in glob.glob(os.path.join(walker_dir, "*.mdp")):
                    replace_all_pattern("TEMP", str(300), mdp_file)
                
                # Replace values in *.dat files
                for dat_file in glob.glob(os.path.join(walker_dir, "plumed*.dat")):
                    for sp_key in statepoint.keys():
                        replace_all_pattern(sp_key, str(statepoint[sp_key]), dat_file)
            
            # Remove original WALKER template directory
            shutil.rmtree(os.path.join(job.path, "WALKER"))
            job.init()

if __name__ == "__main__":
    main()

