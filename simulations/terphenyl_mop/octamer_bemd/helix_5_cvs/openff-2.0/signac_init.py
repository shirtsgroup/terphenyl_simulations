import signac
import numpy as np
import os
import shutil
import re
import glob
import itertools
import sys

def replace_all_pattern(pattern, replace, file):
    with open(file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = re.sub(pattern, replace, lines[i])
    with open(file, "w") as f:
        f.writelines(lines)


def main():
    n_walkers = 5
    project = signac.get_project()
    # combinations

    cv_limits = [[0, 7], [0, 7], [0, 6], [0, 14], [0, 3]]
    statepoints = [
        {
            'height1' : 0,    'height2' : 0.5, 'height3' : 0.5, 'height4' : 0.5, 'height5' : 0.5,
            'sigma1'  : 0,    'sigma2'  : 0.7, 'sigma3'  : 0.1, 'sigma4'  : 0.4, 'sigma5'  : 0.09,
            'bf1'     : 10e9, 'bf2'     : 500, 'bf3'     : 100, 'bf4'     : 100, 'bf5'     : 100,
            'CV_TEMP' : 250
         }
    ]
    replicas = 3

    for statepoint in statepoints:
        for r in range(replicas):
            statepoint['simulation_rep'] = r
            sp_dict = statepoint

            job = project.open_job(sp_dict)

            # Skip statepoint if it already exists
            if job.path.split("/")[-1] in os.listdir("workspace"):
                continue

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
                    replace_all_pattern("TEMP", str(sp_dict["CV_TEMP"]), mdp_file)
                
                # Replace values in *.dat files
                for dat_file in glob.glob(os.path.join(walker_dir, "plumed*.dat")):
                    for sp_key in statepoint.keys():
                        replace_all_pattern(sp_key, str(statepoint[sp_key]), dat_file)

                replace_all_pattern("WALKER_ID", str(walker_id), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                replace_all_pattern("CV_ID", "cv" + str(walker_id + 1) if walker_id == 0 else "cv" + str(walker_id), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                replace_all_pattern("CV_SIGMA", str(statepoint["sigma" + str(walker_id + 1)]), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                replace_all_pattern("CV_BF", str(statepoint["bf" + str(walker_id + 1)]), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                replace_all_pattern("CV_MIN", str(cv_limits[walker_id][0]), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                replace_all_pattern("CV_MAX", str(cv_limits[walker_id][1]), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                replace_all_pattern("SIM_TEMP", str(sp_dict["CV_TEMP"]), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))


            # Remove original WALKER template directory
            shutil.rmtree(os.path.join(job.path, "WALKER"))
            job.init()

if __name__ == "__main__":
    main()

