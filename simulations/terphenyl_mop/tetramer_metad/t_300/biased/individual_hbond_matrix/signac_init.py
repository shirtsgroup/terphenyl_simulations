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
    heights = np.linspace(2.5*2.5, 2.5*5, 2) # height ranging from 1 kT to 5 kT
    sigmas = np.linspace(0.1, 0.5, 5) # A couple different sigmas

    for h in heights:
        for s in sigmas:
            job = project.open_job({'height':h, 'sigma':s})
            # setup input files for metadynamics simulation
            shutil.copytree("simulation_template", job.path)
            replace_all_pattern("height", str(h), os.path.join(job.path, "plumed.dat"))
            replace_all_pattern("sigma", str(s), os.path.join(job.path, "plumed.dat"))
            job.init()

if __name__ == "__main__":
    main()
