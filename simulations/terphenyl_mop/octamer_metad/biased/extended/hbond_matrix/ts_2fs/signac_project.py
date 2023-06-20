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
import mdtraj as md
from natsort import natsorted
from flow import FlowProject
import heteropolymer_simulations as hs

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


@FlowProject.operation
def remove_backup_files(job):
    backup_files = glob.glob(job.fn("WALKER*/#*#"))
    backup_files += glob.glob(job.fn("WALKER*/bck.*"))
    backup_files += glob.glob(job.fn("bck.*"))
    for file in backup_files:
        print("Removing", os.path.abspath(file))
        os.remove(file)

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
    subprocess.run(["plumed", "sum_hills", "--hills", "HILLS", "--kt", str(kt)]) # Silences output from sum_hills
    filename = "fes.dat"
    with open(filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    fes_data = pd.read_csv(filename, skiprows = 5, delim_whitespace=True, header = None, comment="#", names = headers)
    print(fes_data)
    filtered = fes_data[fes_data["n_hbonds"] >= -0.1]
    filtered = filtered[filtered["n_hbonds"] <= 7.2]
    plt.figure(figsize = [5, 2.5])
    plt.xlim([0,7])
    plt.plot(filtered.values[:,0], filtered.values[:,1] - np.min(filtered.values[:,1]))
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
    hs.clustering.clustering_grid_search(
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

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("WALKER0/npt_new.whole.xtc")
@FlowProject.operation
def write_pbc_whole_rep(job):
    """
    Use Gromacs to write trajectories with whole molecules
    """

    print("Running write_pbc_whole_rep in", job.fn(""))

    # Make sure correct version of gromacs is available
    subprocess.Popen("module load gcc/10.2.0", shell=True).wait()
    subprocess.Popen("module load openmpi/3.1.6-gcc10.2.0", shell=True).wait()
    subprocess.Popen("source /jet/home/tfobe/software/gromacs-2021.5_plumed/bin/GMXRC", shell=True).wait()
    
    os.chdir(job.fn(""))
    print("Converting trajectories...")
    for walker_dir in tqdm(glob.glob("WALKER*")):
        os.chdir(walker_dir)
        subprocess.Popen("mpirun -np 1 gmx_mpi trjconv -f npt_new.xtc -s npt_new.tpr -pbc whole -o npt_new.whole.xtc <<<0", shell=True).wait()
        os.chdir(job.fn(""))

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("WALKER0/HBOND_SUMS")
@FlowProject.operation
def reweight_simulations(job):
    """
    Use PLUMED to unbias the individual simulations 
    """
    current_dir = os.path.abspath("")
    kt = 300 * 8.314462618 * 10 ** -3
    os.chdir(job.fn(""))
    
    print("Reweighting simulations...")
    for walker_dir in tqdm(glob.glob("WALKER*")):
        # navigate to specific WALKER dir
        os.chdir(walker_dir)

        # Run plumed driver to reweight biases of invidual simulations
        # We also use this to get H-bond measures from each frame from npt_new.xtc
        subprocess.run(["plumed", "driver", "--plumed", "plumed_reweight.dat", "--kt", str(kt), "--mf_xtc", "npt_new.xtc", "--igro", "npt_new.gro"]) # Silences output from sum_hills
        os.chdir(job.fn(""))


@FlowProject.pre.isfile("WALKER0/HBOND_SUMS")
@FlowProject.post(lambda job: os.path.isdir(job.fn("hbond_states")))
@FlowProject.operation
def write_hb_state_trajectory(job):
    """
    This operation outputs all structures corresponding to specific states. This
    is a quick check to make sure 0 H-bonds and 7 H-bond structures are the structure
    we expect
    """

    print("Running write_hb_state_trajectory in", job.fn(""), "...")

    current_dir = os.path.abspath("")
    kt = 300 * 8.314462618 * 10 ** -3
    os.chdir(job.fn(""))

    # List where trajectories will aggregate
    n_hbond_trajs = {}
    
    print("Extracting H-bond states...")
    for walker_dir in tqdm(glob.glob("WALKER*")):
        # navigate to specific WALKER dir
        os.chdir(walker_dir)

        # Read reweighted H-bond sums (because these comes are calucated from the output frames)
        filename = "HBOND_SUMS"
        with open(filename) as f:
            headers = f.readline().strip()
        headers = headers.split(" ")[1:]
        walker_data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
        
        # We round states to closest integer
        total_hbonds =  np.rint(walker_data.values[:,-2])
        total_hbonds = np.array(total_hbonds, dtype = np.int64)

        # Load npt_new.whole.xtc
        walker_traj = md.load("npt_new.whole.xtc", top = "npt_new.gro")

        # Accumulate trajectories of each state in n_hbond_trajs dict
        for hb_state in np.unique(total_hbonds):
            hbond_indices = np.where(total_hbonds == hb_state)
            hb_slice_traj = walker_traj[hbond_indices]
            if hb_state in n_hbond_trajs.keys():
                n_hbond_trajs[hb_state] = n_hbond_trajs[hb_state].join(hb_slice_traj)
            else:
                n_hbond_trajs[hb_state] = hb_slice_traj
        
        os.chdir(job.fn(""))


    # Make output directory
    print("Writing trajectories...")
    hs.utils.make_path("hbond_states")

    # Write all hbond states to file
    for hb_state in n_hbond_trajs.keys():
        n_hbond_trajs[hb_state].save_xtc("hbond_states/" + "hb_state_" + str(hb_state) + ".xtc")



@FlowProject.operation
def show_statepoint_table(job):
    print("sp:", job.sp, "dir:", job.fn(""), "status:", check_production_npt_finish(job))


if __name__ == '__main__':
    FlowProject().main()
