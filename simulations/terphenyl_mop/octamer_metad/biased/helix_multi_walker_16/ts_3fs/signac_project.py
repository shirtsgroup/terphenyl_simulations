import signac
import os
import sys
import glob
import flow
import subprocess
import os
import matplotlib.pyplot as plt
import glob
import plumed
import numpy as np
from tqdm import tqdm
from natsort import natsorted
from flow import FlowProject

# Scripts for running the submit_all.slurm script which submits all simulations
# at once with dependencies linking them. If this fails individual submission
# flow operations can be found below.

@FlowProject.label
def check_submitted(job):

    if job.isfile("status.txt"):
        with open(job.fn("status.txt"), "r") as f:
            lines = f.readlines()
        
        if "SUBMITTED" in lines:
            return True
    return False


@FlowProject.post(check_submitted)
@FlowProject.operation
def submit_all_simulations(job):
    os.chdir(job.path)
    n_jobs_old = len(subprocess.check_output(["squeue", "-u", "tfobe"]).splitlines()) - 1
    subprocess.run(["bash", "submit_all.slurm"])
    n_jobs = len(subprocess.check_output(["squeue", "-u", "tfobe"]).splitlines()) - 1
    jobs_submitted = n_jobs - n_jobs_old
    if jobs_submitted == 5:
        with open(job.fn("status.txt"), "w") as f:
            f.write("SUBMITTED")
    else:
        with open(job.fn("status.txt"), "w") as f:
            f.write("FAILED")

# Labels for individual simulation components

def check_walker_file(job, filename, walker_dirs = ["WALKER0", "WALKER1", "WALKER2", "WALKER3"]):
    walkers = []
    for walker_dir in walker_dirs:
        walkers.append(job.isfile(os.path.join(walker_dir, filename)))
    return all(walkers)

@FlowProject.label
def check_berendsen_nvt_start(job):
    return check_walker_file(job, "berendsen_nvt.log")

@FlowProject.label
def check_berendsen_nvt_finish(job):
    return check_walker_file(job, "berendsen_nvt.gro")

@FlowProject.label
def check_berendsen_npt_start(job):
    return check_walker_file(job, "berendsen_npt.log")

@FlowProject.label
def check_berendsen_npt_finish(job):
    return check_walker_file(job, "berendsen_npt.gro")

@FlowProject.label
def check_production_npt_start(job):
    return check_walker_file(job, "npt_new.log")

@FlowProject.label
def check_production_npt_finish(job):
   return check_walker_file(job, "npt_new.gro")

@FlowProject.post(check_berendsen_nvt_finish)
@FlowProject.operation
def submit_berendsen_nvt_simulations(job):
    os.chdir(job.path)
    subprocess.run(["sbatch", "submit.berendsen_nvt.slurm"])

@FlowProject.pre(check_berendsen_nvt_finish)
@FlowProject.post(check_berendsen_npt_finish)
@FlowProject.operation
def submit_berendsen_npt_simulations(job):
    os.chdir(job.path)
    subprocess.run(["sbatch", "submit.berendsen_npt.slurm"])

@FlowProject.pre(check_berendsen_npt_finish)
@FlowProject.post(check_production_npt_finish)
@FlowProject.operation
def submit_production_npt_simulations(job):
    os.chdir(job.path)
    slurm_id = subprocess.check_output(["sbatch", "submit.production.slurm"])

@FlowProject.pre(check_production_npt_finish)
@FlowProject.operation
def continue_production_npt_simulation(job):
    os.chdir(job.path)
    slurm_id = subprocess.check_output(["sbatch", "submit.continue.slurm"])

# Utility operation to remove all backup files
@FlowProject.label
def has_backup_files(job):
    backup_files = glob.glob(job.fn("WALKER*/#*#"))
    backup_files += glob.glob(job.fn("WALKER*/bck.*"))
    if len(backup_files) == 0:
        return False
    return True

@FlowProject.label
def has_failed_step_files(job):
    step_files = glob.glob(job.fn("WALKER*/step*.pdb"))
    if len(step_files) == 0:
        return False
    return True

@FlowProject.pre(has_backup_files)
@FlowProject.post.not_(has_backup_files)
@FlowProject.operation
def remove_backup_files(job):
    backup_files = glob.glob(job.fn("WALKER*/#*#"))
    backup_files += glob.glob(job.fn("WALKER*/bck.*"))
    for file in backup_files:
        print("Removing", os.path.abspath(file))
        os.remove(file)

@FlowProject.pre(has_failed_step_files)
@FlowProject.post.not_(has_failed_step_files)
@FlowProject.operation
def remove_failed_step_files(job):
    step_files = glob.glob(job.fn("WALKER*/step*.pdb"))
    for file in step_files:
        print("Removing", os.path.abspath(file))
        os.remove(file)

@FlowProject.pre(check_production_npt_start)
@FlowProject.post.isfile("CV_bias_plot.png")
@FlowProject.operation
def plot_CV_bias(job):
    stride = 100
    walker_dirs = glob.glob(job.fn("WALKER*"))
    walker_ids = [int(walker_dir.split("WALKER")[-1]) for walker_dir in walker_dirs]
    plt.figure(dpi=150)
    fig, ax = plt.subplots(len(walker_dirs) + 2, 1, figsize = [10,30])
    ax[0].set_title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    all_CV = []
    for i, walker_dir in tqdm(zip(walker_ids, walker_dirs)):
        walker_data = plumed.read_as_pandas(walker_dir + "/HBOND_SUMS."+str(i))
        all_CV += list(walker_data.values[:,-2][::stride])
        ax[i].plot(1/1000 * walker_data.time[::stride], walker_data.values[:,-2][::stride])
        ax[i].set_ylabel("Replica " + str(i) + " $N_H$")
        ax[len(walker_dirs)].plot(1/1000 * walker_data.time[::stride], walker_data.values[:,-1][::stride])
    hills_data = plumed.read_as_pandas(job.fn("HILLS"))
    ax[len(walker_dirs)].set_ylabel("Bias Energy (kJ/mol)")
    ax[len(walker_dirs)+1].plot(1/1000 * hills_data.time, hills_data.values[:,3], "o", markersize = 0.5)
    ax[len(walker_dirs)+1].set_ylabel("Gaussian Heights (kJ/mol)")
    ax[len(walker_dirs)+1].set_xlabel("Time (ns)")
    plt.savefig(job.fn("CV_bias_plot.png"), dpi = 150)
    plt.close()
    plt.figure(figsize = [5, 2.5], dpi = 150)
    plt.hist(all_CV, bins = 100, density = True)
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("$N_H$")
    plt.ylabel("Probability")
    plt.savefig(job.fn("CV_sampling.png"), dpi = 150, transparent = False)
    plt.close()



@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("sum_hills_FE.png")
@FlowProject.operation
def calculate_sum_hills_FE(job):
    current_dir = os.path.abspath("")
    kt = 300 * 8.314462618 * 10 ** -3
    os.chdir(job.fn(""))
    subprocess.run(["plumed", "sum_hills", "--hills", "HILLS", "--kt", str(kt)]) # Silences output from sum_hills
    fes_data = plumed.read_as_pandas("fes.dat")
    print(fes_data)
    filtered = fes_data[fes_data["n_hbonds"] >= -0.1]
    filtered = filtered[filtered["n_hbonds"] <= 7.2]
    plt.figure(figsize = [5, 2.5])
    plt.xlim([0,7])
    plt.plot(filtered.values[:,0], filtered.values[:,1])
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("$N_{H-bonds}$")
    plt.ylabel("Free Energy (kJ/mol)")
    plt.grid(visible=True, which="both", axis="both")
    plt.savefig("sum_hills_FE.png", transparent = False)
    plt.close()
    os.chdir(current_dir)

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("h_bond_transition_matrix.png")
@FlowProject.operation
def plot_transition_matrix(job):
    stride = 100
    walker_dirs = glob.glob(job.fn("WALKER*"))
    walker_ids = [int(walker_dir.split("WALKER")[-1]) for walker_dir in walker_dirs]
    transition_matrix = np.zeros([8,8])
    for i, walker_dir in tqdm(zip(walker_ids, walker_dirs)):
        walker_data = plumed.read_as_pandas(walker_dir + "/HBOND_SUMS."+str(i))
        h_bond_states = np.rint(walker_data.values[:,-2][::stride])
        for i in range(len(h_bond_states)-1):
            transition_matrix[int(h_bond_states[i]), int(h_bond_states[i+1])] += 1
    
    transition_matrix = transition_matrix / np.sum(transition_matrix, axis = 0)
    fig, ax = plt.subplots(figsize=[5,5])
    ax.matshow(transition_matrix, cmap = 'Greens')

    for (i, j), z in np.ndenumerate(transition_matrix):
        ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
    ax.set_xlabel("$N_{H-bonds}$")
    ax.set_ylabel("$N_{H-bonds}$")

    plt.savefig(job.fn("h_bond_transition_matrix.png"))

@FlowProject.operation
def show_statepoint_table(job):
    print("sp:", job.sp, "dir:", job.fn(""), "status:", check_production_npt_finish(job))


if __name__ == '__main__':
    FlowProject().main()