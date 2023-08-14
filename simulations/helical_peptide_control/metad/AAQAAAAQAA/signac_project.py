import signac
import os
import sys
import glob
import flow
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import glob
from natsort import natsorted
from flow import FlowProject
import terphenyl_simulations as ts

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
    backup_files += glob.glob(job.fn("bck.*"))
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
    backup_files += glob.glob(job.fn("bck.*"))
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
    fig, ax = plt.subplots(len(walker_dirs) + 2, 1, figsize = [20, 5 * len(walker_dirs) + 5])
    ax[0].set_title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    all_CV = []
    for i, walker_dir in zip(walker_ids, walker_dirs):
        filename = walker_dir + "/HBOND_SUMS."+str(i)
        with open(filename) as f:
            headers = f.readline().strip()
        headers = headers.split(" ")[1:]
        walker_data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
        all_CV += list(walker_data.values[:,-2][::stride])
        ax[i].plot(1/1000 * walker_data["time"][::stride], walker_data.values[:,-2][::stride])
        ax[i].set_ylabel("Replica " + str(i) + " $N_H$")
        ax[len(walker_dirs)].plot(1/1000 * walker_data["time"][::stride], walker_data.values[:,-1][::stride])
    hills_file = job.fn("HILLS")
    with open(hills_file) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    hills_data = pd.read_csv(hills_file, skiprows = 3, delim_whitespace=True, header = None, comment="#", names = headers)
    ax[len(walker_dirs)].set_ylabel("Bias Energy (kJ/mol)")
    ax[len(walker_dirs)+1].plot(1/1000 * hills_data["time"], hills_data["height"], "o", markersize = 0.5)
    ax[len(walker_dirs)+1].set_ylabel("Gaussian Heights (kJ/mol)")
    ax[len(walker_dirs)+1].set_xlabel("Time (ns)")
    plt.savefig(job.fn("CV_bias_plot.png"), dpi = 300)
    plt.close()
    plt.figure(figsize = [5, 2.5], dpi = 300)
    plt.hist(all_CV, bins = 100, density = True)
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("$N_H$")
    plt.ylabel("Probability")
    plt.savefig(job.fn("CV_sampling.png"), dpi = 300, transparent = False)
    plt.close()



@FlowProject.pre(check_production_npt_start)
@FlowProject.post.isfile("sum_hills_FE.png")
@FlowProject.operation
def calculate_sum_hills_FE(job):
    current_dir = os.path.abspath("")
    kt = 300 * 8.314462618 * 10 ** -3
    os.chdir(job.fn(""))

    # FE surface as a function of time
    ts.utils.make_path("free_energy")
    os.chdir("free_energy")
    subprocess.run(["plumed", "sum_hills", "--hills", "../HILLS", "--kt", str(kt), "--stride", "160000"]) # Silences output from sum_hills
    fes_files = natsorted(glob.glob("fes*.dat"))
    plt.figure(figsize = [5, 2.5])
    for fes_file in fes_files:
        print(fes_file)
        filename = fes_file
        name = filename.split(".")[0]
        with open(filename) as f:
            headers = f.readline().strip()
        headers = headers.split(" ")[2:]
        fes_data = pd.read_csv(filename, skiprows = 5, delim_whitespace=True, header = None, comment="#", names = headers)
        print(fes_data)
        filtered = fes_data[fes_data["n_hbonds"] >= -0.1]
        filtered = filtered[filtered["n_hbonds"] <= 4.2]
        plt.plot(filtered.values[:,0], filtered.values[:,1])
    plt.xlim([0,4])
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("$N_{H-bonds}$")
    plt.ylabel("Free Energy (kJ/mol)")
    plt.grid(visible=True, which="both", axis="both")
    plt.savefig("sum_hills_all.png", transparent = False, dpi = 300)
    plt.close()

    # Overall FES
    os.chdir("..")
    plt.figure(figsize = [5, 2.5])
    plt.plot(filtered.values[:,0], filtered.values[:,1] - min(filtered.values[:,1]))
    plt.xlim([0,4])
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("$N_{H-bonds}$")
    plt.ylabel("Free Energy (kJ/mol)")
    plt.grid(visible=True, which="both", axis="both")
    plt.savefig("sum_hills_FE.png", transparent = False, dpi = 300)
    plt.close()
    os.chdir(current_dir)

@FlowProject.pre(check_production_npt_start)
@FlowProject.post.isfile("h_bond_transition_matrix.png")
@FlowProject.operation
def plot_transition_matrix(job):
    stride = 100
    walker_dirs = glob.glob(job.fn("WALKER*"))
    transition_matrix = np.zeros([8,8])
    for walker_dir in tqdm(walker_dirs):
        walker_id = int(walker_dir.split("WALKER")[-1])
        filename = walker_dir + "/HBOND_SUMS."+str(walker_id)
        with open(filename) as f:
            headers = f.readline().strip()
        headers = headers.split(" ")[1:]
        walker_data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
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

    plt.savefig(job.fn("h_bond_transition_matrix.png"), dpi = 300)

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("clustering_output/medoid_0.gro")
@FlowProject.operation
def cluster_simulations(job):
    os.chdir(job.fn(""))
    n_walkers = len(glob.glob("WALKER*"))
    ts.clustering.clustering_grid_search(
        ["WALKER" + str(i) + "/npt_new.whole.xtc" for i in range(n_walkers)],
        "WALKER0/berendsen_npt.gro",
        "resname OCT or resname CAP",
        n_min_samples=20,
        n_eps=20,
        n_processes=32,
        prefix="grid_search",
        eps_limits=[0.1, 0.5],
        min_sample_limits=[0.005, 0.15],
        plot_filename = "ss.png",
        frame_stride = 10,
    )

@FlowProject.operation
def show_statepoint_table(job):
    print("sp:", job.sp, "dir:", job.fn(""), "status:", check_production_npt_finish(job))



if __name__ == '__main__':
    FlowProject().main()
