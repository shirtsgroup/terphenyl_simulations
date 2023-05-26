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
            replace_all_pattern("height", str(statepoint['height']), job.fn("plumed_hbond_dist.dat"))
            replace_all_pattern("sigma", str(statepoint['sigma']), job.fn("plumed_hbond_dist.dat"))
            replace_all_pattern("bf", str(statepoint['bf']), job.fn("plumed_hbond_dist.dat"))
            replace_all_pattern("bf", str(statepoint['bf']), job.fn("plumed_reweight.dat"))
            replace_all_pattern("sigma", str(statepoint['bf']), job.fn("plumed_reweight.dat"))
            replace_all_pattern("TEMP", str(300), job.fn("berendsen_nvt.mdp"))
            replace_all_pattern("TEMP", str(300), job.fn("berendsen_npt.mdp"))
            replace_all_pattern("TEMP", str(300), job.fn("npt_new.mdp"))
            replace_all_pattern("temp", str(300), job.fn("plumed_hbond_dist.dat"))
            replace_all_pattern("temp", str(300), job.fn("plumed_reweight.dat"))
            job.init()

if __name__ == "__main__":
    main()

