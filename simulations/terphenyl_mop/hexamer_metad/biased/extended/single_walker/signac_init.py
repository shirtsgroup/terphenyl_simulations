import signac
import numpy as np
import os
import shutil
import re


def replace_all_pattern(pattern, replace, file):
    with open(file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = re.sub(pattern, replace, lines[i])
    with open(file, "w") as f:
        f.writelines(lines)


def main():
    project = signac.get_project()
    heights = [1.0, 2.5] # height ranging from 1 kT to 3 kT
    sigmas = [0.5] # A couple different sigmas
    bias_factors = [50, 100, 150, 200, 1e9]
    temperatures = np.array([300])
    replicas = 5

    for h in heights:
        for s in sigmas:
            for bf in bias_factors:
                for r in range(replicas):
                    sp_dict = {'height':float(h), 'sigma':float(s), 'bf':int(bf), 'replica':int(r)}
                    if len(project.find_jobs(sp_dict).to_dataframe()) == 0:
                        job = project.open_job(sp_dict)
                        # setup input files for metadynamics simulation
                        shutil.copytree("simulation_template", job.path)
                        replace_all_pattern("height", str(h), os.path.join(job.path, "plumed.dat"))
                        replace_all_pattern("sigma", str(s), os.path.join(job.path, "plumed.dat"))
                        replace_all_pattern("bf", str(bf), os.path.join(job.path, "plumed.dat"))
                        replace_all_pattern("bf", str(bf), os.path.join(job.path, "plumed_reweight.dat"))
                        replace_all_pattern("sigma", str(bf), os.path.join(job.path, "plumed_reweight.dat"))
                        replace_all_pattern("TEMP", str(300), os.path.join(job.path, "berendsen_nvt.mdp"))
                        replace_all_pattern("TEMP", str(300), os.path.join(job.path, "berendsen_npt.mdp"))
                        replace_all_pattern("TEMP", str(300), os.path.join(job.path, "npt_new.mdp"))
                        job.init()
                    else:
                        print("This statepoint is already defined. Skipping.", sp_dict)

if __name__ == "__main__":
    main()
