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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import terphenyl_simulations as ts

# Helper Functions

def check_walker_file(job, filename, walker_dirs = ["WALKER0", "WALKER1", "WALKER2", "WALKER3"]):
    walkers = []
    for walker_dir in walker_dirs:
        walkers.append(job.isfile(os.path.join(walker_dir, filename)))
    return all(walkers)

def read_plumed_data_file(filename):
    """
    Function for reading output from PRINT operations in plumed
    """
    with open(filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[1:]
    data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
    return data

def read_plumed_hills_file(hills_file):
    """
    Function for reading HILLS output file from plumed driver
    """
    with open(hills_file) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    hills_data = pd.read_csv(hills_file, skiprows = 3, delim_whitespace=True, header = None, comment="#", names = headers)
    return hills_data

def read_plumed_fes_file(fes_filename):
    """
    Function for reading FES output from plumed sum_hills function
    """
    with open(fes_filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    fes_data = pd.read_csv(fes_filename, skiprows = 5, delim_whitespace=True, header = None, comment="#", names = headers)
    return fes_data

def reweight_walker_trajectories(job, plumed_file, kt, gro_file = "npt_new.gro", xtc_file = "npt_new.xtc"):
    # Get original directory path
    current_dir = os.path.abspath("")

    # Change directory to job
    os.chdir(job.fn(""))
    
    print("Reweighting simulations...")
    for walker_dir in tqdm(glob.glob("WALKER*")):
        # navigate to specific WALKER dir
        os.chdir(walker_dir)

        # Run plumed driver to reweight biases of invidual simulations
        # We also use this to get H-bond measures from each frame from npt_new.xtc
        subprocess.run(["plumed", "--no-mpi", "driver", "--plumed", plumed_file, "--kt", str(kt), "--mf_xtc", xtc_file, "--igro", gro_file]) # Silences output from sum_hills
        os.chdir(job.fn(""))
    
    # Return to original directory
    os.chdir(current_dir)

# Signac Labels

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

# Operations for running individual portions of simulation

@FlowProject.post(check_berendsen_nvt_start)
@FlowProject.operation
def submit_all_simulations(job):
    os.chdir(job.path)
    n_jobs_old = len(subprocess.check_output(["squeue", "-u", "tfobe"]).splitlines()) - 1
    subprocess.run(["bash", "submit_all.slurm"])

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

# Analysis operations

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("bemd_state_index.png")
@FlowProject.operation
def plot_state_index_plot(job):
    os.chdir(job.path)
    ts.plotting.plot_bemd_state_index_plot("WALKER0/npt_new.log", "bemd")

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("bemd_transition_matrix.png")
@FlowProject.operation
def plot_transition_matrix(job):
    os.chdir(job.path)
    ts.plotting.plot_bemd_transition_matrix("WALKER0/npt_new.log", "bemd")

@FlowProject.operation
def show_statepoint_table(job):
    print("sp:", job.sp, "dir:", job.fn(""), "status:", check_production_npt_finish(job))

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("CV_sampling.png")
@FlowProject.operation
def plot_multi_CV_plot(job):
    stride = 100
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))
    fig, axs = plt.subplots(nrows = len(walker_dirs), ncols = len(walker_dirs)+1, figsize = [20,10])
    overall_cv1 = []
    overall_cv2 = []
    overall_cv3 = []
    for i, walker_dir in tqdm(enumerate(walker_dirs)):
        colvars_file = os.path.join(walker_dir, "COLVARS."+str(i))
        colvars = read_plumed_data_file(colvars_file)
        overall_cv1 += list(colvars["cv1"][::stride].values)
        overall_cv2 += list(colvars["cv2"][::stride].values)
        overall_cv3 += list(colvars["cv3"][::stride].values)
        axs[0,i].plot(1/1000 * colvars["time"][::stride], colvars["cv1"][::stride], markersize = 0.5)
        axs[1,i].plot(1/1000 * colvars["time"][::stride], colvars["cv2"][::stride], markersize = 0.5)
        axs[2,i].plot(1/1000 * colvars["time"][::stride], colvars["cv3"][::stride], markersize = 0.5)
        axs[0,i].set_xlabel("Time (ns)")
        axs[0,i].set_ylabel("Native Contacts")
        axs[0,i].set_title("Simulation " + str(i))
        axs[1,i].set_xlabel("Time (ns)")
        axs[1,i].set_ylabel("Radius of Gyration (nm)")
        axs[2,i].set_xlabel("Time (ns)")
        axs[2,i].set_ylabel("End-to-end distance (nm)")
        axs[3,i].plot(1/1000 * colvars["time"][::stride], colvars.values[:,-1][::stride], markersize = 0.5, color = "red")
        axs[3,i].set_xlabel("Time (ns)")
        axs[3,i].set_ylabel("Bias (kJ/mol)")
    axs[0,-1].hist(overall_cv1, bins = 50, density = True)
    axs[0,-1].set_xlabel("Native Contacts")
    axs[0,-1].set_ylabel("Probability Density")
    axs[1,-1].hist(overall_cv2, bins = 50, density = True)
    axs[1,-1].set_xlabel("Radius of Gyration (nm)")
    axs[1,-1].set_ylabel("Probability Density")
    axs[2,-1].hist(overall_cv3, bins = 50, density = True)
    axs[2,-1].set_xlabel("End-to-end distance (nm)")
    axs[2,-1].set_ylabel("Probability Density")
    axs[0,-1].set_title("Overall PDF")
    axs[-1, -1].remove()
    

    plt.tight_layout()
    plt.savefig(job.fn("CV_sampling.png"), dpi = 300)

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("multi_sum_hills_FE.png")
@FlowProject.operation
def sum_hills_FE(job):
    # Navigate to working directory
    current_dir = os.path.abspath("")
    
    # Parameters about simulations
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))
    n_simulations = len(walker_dirs)
    kt = 300 * 8.314462618 * 10 ** -3
    x_labels = ["$N_{H-bonds}$", "$R_{G}$", "$d_{end-to-end}$"]
    xy_limits = [[[0,7], [0, 1000]], [[0.75, 3], [0, 500]], [[0, 6], [0, 500]]]

    
    # Setup subplots
    n_row_columns = int(np.ceil(np.sqrt(n_simulations)))
    fig, axes = plt.subplots(nrows = n_row_columns, ncols = n_row_columns, figsize = [10, 10])
    axes = np.roll(axes.flatten(), 1, axis=0).tolist()

    for i, (walker_dir, ax) in enumerate(zip(walker_dirs, axes)):
        if "WALKER0" in walker_dir:
            fig.delaxes(ax)
            continue
        os.chdir(walker_dir)
        walker_id = int(walker_dir.split("WALKER")[-1])
        subprocess.run(["plumed", "--no-mpi", "sum_hills", "--hills", "HILLS." + str(walker_id), "--kt", str(kt), "--mintozero"]) # Silences output from sum_hills
        filename = "fes.dat"
        fes_data = read_plumed_fes_file(filename)
        ax.plot(fes_data.values[:,0], fes_data.values[:,1])
        ax.set_xlabel(x_labels[i - 1])
        ax.set_ylabel("Free Energy (kJ/mol)")
        ax.set_xlim(xy_limits[i-1][0])
        ax.set_ylim(xy_limits[i-1][1])
        ax.set_title("Simulation " + str(walker_id))
        ax.grid(visible=True, which="both", axis="both")
    os.chdir(current_dir)
    plt.savefig(job.fn("multi_sum_hills_FE.png"), transparent = False, dpi = 300)

@FlowProject.pre(check_production_npt_finish)
@FlowProject.operation
def plot_2D_FE_surface(job):
    # Need to modify HILLS filename to HILL.x
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))
    for i, walker_dir in enumerate(walker_dirs):
        plumed_file = os.path.join(walker_dir, "plumed_multi_cv_reweight.dat")
        with open(plumed_file, "r") as fn:
            for line in fn.readlines():
                if "FILE=HILLS" in line and line.strip().split("=")[-1] == "HILLS":
                    ts.utils.replace_all_pattern("HILLS", "HILLS." + str(i), os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
                if "ARG=@replicas" in line:
                    cv_label = "cv" + str(i)
                    if i == 0:
                        cv_label = "cv" + str(i + 1)
                    ts.utils.replace_all_pattern("@replicas:cv1,cv1,cv2,cv3", cv_label, os.path.join(walker_dir, "plumed_multi_cv_reweight.dat"))
 

    # Aggregate all reweighted samples with their bias
    kt = 300 * 8.314462618 * 10 ** -3
    agg_data = None
    reweight_walker_trajectories(job, "plumed_multi_cv_reweight.dat", kt)
    for i, walker_dir in enumerate(walker_dirs):
        reweight_file = os.path.join(walker_dir, "COLVARS_REWEIGHT")
        if agg_data is None:
            agg_data = read_plumed_data_file(reweight_file)
        else:
            agg_data = pd.concat([agg_data, read_plumed_data_file(reweight_file)])
    
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])

    hist, xedges, yedges = np.histogram2d(agg_data.cv1, agg_data.cv2, weights=agg_data.values[:,-1], bins = [15, 30], density = True)
    fe = -kt * np.log(hist)
    axes[0,0].imshow(fe, extent = [xedges.min(), xedges.max(), yedges.min(), yedges.max()])
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("$R_{G}$")
    divider = make_axes_locatable(axes[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(hist[3], cax=cax, orientation='vertical')
    sys.exit()

    hist = axes[1,0].hist2d(agg_data.cv1, agg_data.cv3, weights=agg_data.values[:,-1], bins = [15, 30], density = True)
    axes[1,0].set_xlabel("$N_{H-bonds}$")
    axes[1,0].set_ylabel("$d_{end-to-end}$")
    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(hist[3], cax=cax, orientation='vertical')


    hist = axes[0,1].hist2d(agg_data.cv3, agg_data.cv2, weights=agg_data.values[:,-1], bins = [30, 30], density = True)
    axes[0,1].set_ylabel("$R_{G}$")
    axes[0,1].set_xlabel("$d_{end-to-end}$")
    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(hist[3], cax=cax, orientation='vertical')

    plt.delaxes(axes[1,1])

    plt.savefig(job.fn("2D_FE_surfaces.png"))



def main():
    FlowProject().main()
