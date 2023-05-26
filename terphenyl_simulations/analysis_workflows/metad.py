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
import terphenyl_simulations as ts

# Scripts for running the submit_all.slurm script which submits all simulations
# at once with dependencies linking them. If this fails individual submission
# flow operations can be found below.

@FlowProject.label
def sw_check_nvt_equilibration(job):
    return(os.path.exists(job.fn("berendsen_nvt.gro")))

@FlowProject.label
def sw_check_npt_equilibration(job):
    return(os.path.exists(job.fn("berendsen_npt.gro")))

@FlowProject.label
def sw_check_npt_production(job):
    return(os.path.exists(job.fn("npt_new.gro")))

@FlowProject.post(sw_check_nvt_equilibration)
@FlowProject.post(sw_check_npt_equilibration)
@FlowProject.post(sw_check_npt_production)
@FlowProject.operation
def run_all_simulations(job):
    os.chdir(job.fn(""))
    print("Working on", job.fn(""))
    process = subprocess.Popen("bash run_local.sh".split(" "))
    process.wait()

@FlowProject.pre(sw_check_npt_production)
@FlowProject.operation
def continue_simulations(job):
    os.chdir(job.fn(""))
    print("Working on", job.fn(""))
    subprocess.Popen("gmx convert-tpr -s npt_new.tpr -extend 100000 -o npt_new.tpr", shell=True).wait()
    subprocess.Popen("gmx mdrun -v -deffnm npt_new -cpi npt_new.cpt -plumed plumed_hbond_dist.dat", shell=True).wait()

@FlowProject.pre(sw_check_npt_production)
@FlowProject.post.isfile("npt_new.whole.xtc")
@FlowProject.operation
def write_pbc_whole_rep(job):
    """
    Use Gromacs to write trajectories with whole molecules
    """

    print("Running write_pbc_whole_rep in", job.fn(""))
    os.chdir(job.fn(""))
    subprocess.Popen("gmx trjconv -f berendsen_nvt.xtc -s npt_new.tpr -pbc whole -o berendsen_nvt.whole.xtc<<<0", shell=True, executable="/bin/bash").wait()
    subprocess.Popen("gmx trjconv -f berendsen_npt.xtc -s npt_new.tpr -pbc whole -o berendsen_npt.whole.xtc<<<0", shell=True, executable="/bin/bash").wait()
    subprocess.Popen("gmx trjconv -f npt_new.xtc -s npt_new.tpr -pbc whole -o npt_new.whole.xtc<<<0", shell=True, executable="/bin/bash").wait()

@FlowProject.operation
def remove_backup_files(job):
    backup_files = glob.glob(job.fn("#*#"))
    backup_files += glob.glob(job.fn("bck.*"))
    backup_files += glob.glob(job.fn("bck.*"))
    for file in backup_files:
        print("Removing", os.path.abspath(file))
        os.remove(file)

@FlowProject.pre(sw_check_npt_production)
@FlowProject.post.isfile("CV_bias_plot.png")
@FlowProject.post.isfile("CV_sampling.png")
@FlowProject.operation
def plot_CV_bias(job):
    stride = 100
    plt.figure(dpi=150)
    fig, ax = plt.subplots(3, 1, figsize = [20, 10])
    ax[0].set_title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    all_CV = []
    filename = job.fn("HBOND_SUMS")
    with open(filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[1:]
    data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
    all_CV += list(data.values[:,-2][::stride])
    ax[0].plot(1/1000 * data["time"][::stride], data.values[:,-2][::stride])
    ax[0].set_ylabel("$N_H$")
    ax[-2].plot(1/1000 * data["time"][::stride], data.values[:,-1][::stride])
    hills_file = job.fn("HILLS")
    with open(hills_file) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    hills_data = pd.read_csv(hills_file, skiprows = 3, delim_whitespace=True, header = None, comment="#", names = headers)
    ax[-2].set_ylabel("Bias Energy (kJ/mol)")
    ax[-1].plot(1/1000 * hills_data["time"], hills_data["height"], "o", markersize = 0.5)
    ax[-1].set_ylabel("Gaussian Heights (kJ/mol)")
    ax[-1].set_xlabel("Time (ns)")
    plt.savefig(job.fn("CV_bias_plot.png"), dpi = 300)
    plt.close()
    plt.figure(figsize = [5, 2.5], dpi = 300)
    plt.hist(all_CV, bins = 100, density = True)
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("$N_H$")
    plt.ylabel("Probability")
    plt.savefig(job.fn("CV_sampling.png"), dpi = 300, transparent = False)
    plt.close()

@FlowProject.pre(sw_check_npt_production)
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

@FlowProject.pre.isfile("HBOND_SUMS")
@FlowProject.pre.isfile("npt_new.whole.xtc")
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

    # Read reweighted H-bond sums (because these comes are calucated from the output frames)
    filename = "HBOND_SUMS"
    with open(filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[1:]
    data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
    
    # We round states to closest integer
    total_hbonds =  np.rint(data.values[:,-2])
    total_hbonds = np.array(total_hbonds, dtype = np.int64)

    # Load npt_new.whole.xtc
    traj = md.load("npt_new.whole.xtc", top = "npt_new.gro")
    n_frames = len(traj)
    n_hbond_states = len(total_hbonds)

    stride = (n_hbond_states - 1) / (n_frames - 1)

    if not stride.is_integer():
        print("Structures output from simulation do not align with Plumed H-bond information.")
        print("Please this operation cannot be complete. Ensure nstxout-compressed in the .mdp")
        print("and the STRIDE option in the plumed.dat file are multiples of one another.")
        return
    
    stride = int(stride)
    hbonds_traj = total_hbonds[::stride]

    # Accumulate trajectories of each state in n_hbond_trajs dict
    for hb_state in np.unique(total_hbonds):
        hbond_indices = np.where(hbonds_traj == hb_state)
        hb_slice_traj = traj[hbond_indices]
        if hb_state in n_hbond_trajs.keys():
            n_hbond_trajs[hb_state] = n_hbond_trajs[hb_state].join(hb_slice_traj)
        else:
            n_hbond_trajs[hb_state] = hb_slice_traj
    
    # Make output directory
    print("Writing trajectories...")
    ts.utils.make_path("hbond_states")

    # Write all hbond states to file
    for hb_state in n_hbond_trajs.keys():
        n_hbond_trajs[hb_state].save_xtc("hbond_states/" + "hb_state_" + str(hb_state) + ".xtc")


    
def main():
    FlowProject().main()