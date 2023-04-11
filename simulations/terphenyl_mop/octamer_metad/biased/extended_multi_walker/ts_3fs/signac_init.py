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
    heights = np.array([0.5, 2.5]) # height ranging from 2.5 kT to 5 kT divided by number of walkers to ensure accumulation of biases are not too large
    sigmas = [0.5] # A couple different sigmas
    bias_factors = [50, 100, 200, 10e9]
    replica = list(range(1))

    for combination in itertools.product(heights, sigmas, bias_factors, replica):
        h, s, bf, r = combination
        sp_dict = {'height':float(h), 'sigma':float(s), 'bf':int(bf), 'replica':r}

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
        replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn("submit.continue.slurm"))
        for walker_id in range(n_walkers):
            walker_dir = os.path.join(job.path, "WALKER" + str(walker_id))
            shutil.copytree(os.path.join(job.path, "WALKER"), walker_dir)
            replace_all_pattern("height", str(h), os.path.join(walker_dir, "plumed.dat"))
            replace_all_pattern("sigma", str(s), os.path.join(walker_dir, "plumed.dat"))
            replace_all_pattern("bf", str(bf), os.path.join(walker_dir, "plumed.dat"))
            replace_all_pattern("bf", str(bf), os.path.join(walker_dir,  "plumed_reweight.dat"))
            replace_all_pattern("sigma", str(bf), os.path.join(walker_dir, "plumed_reweight.dat"))
            replace_all_pattern("TEMP", str(300), os.path.join(walker_dir, "berendsen_nvt.mdp"))
            replace_all_pattern("TEMP", str(300), os.path.join(walker_dir, "berendsen_npt.mdp"))
            replace_all_pattern("TEMP", str(300), os.path.join(walker_dir, "npt_new.mdp"))
            replace_all_pattern("temp", str(300), os.path.join(walker_dir, "plumed.dat"))
            replace_all_pattern("temp", str(300), os.path.join(walker_dir, "plumed_reweight.dat"))
        shutil.rmtree(os.path.join(job.path, "WALKER"))
        job.init()

if __name__ == "__main__":
    main()

