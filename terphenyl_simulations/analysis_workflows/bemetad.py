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

def reweight_colvars(job, plumed_file, colvar_file):
    current_dir = os.path.abspath("")


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
    print("sp:", job.sp, "\n", "dir:", job.fn("").split("/")[-2], "\n", "status:", check_production_npt_finish(job), "\n")

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
    xy_limits = [[[0,7], [0, 500]], [[0.75, 3], [0, 100]], [[0, 6], [0, 150]]]
    plumed_bins = [[0, 7], [0, 6], [0, 14]]

    
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
@FlowProject.post.isfile("WALKER0/COLVARS_REWEIGHT")
@FlowProject.operation
def unbias_simulations(job):

    kt = 300 * 8.314462618 * 10 ** -3
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3

    reweight_walker_trajectories(job, "plumed_multi_cv_reweight.dat", kt)


@FlowProject.pre.isfile("WALKER0/COLVARS_REWEIGHT")
@FlowProject.operation
def plot_1D_FE_surface(job):
    # Navigate to working dir
    os.chdir(job.fn(""))
    ts.utils.make_path("1D_FE_plots")
    walker_dirs = natsorted(glob.glob("WALKER*"))


    # Calculate kT for specific temperatures
    kt = 300 * 8.314462618 * 10 ** -3
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3

    # Generate common bins for all CVs
    n_bins = 30

    cv1_bins = np.linspace(0, 7, n_bins)
    cv1_centers = [(cv1_bins[i] + cv1_bins[i+1])/2 for i in range(len(cv1_bins)-1)]
    cv2_bins = np.linspace(0, 4, n_bins)
    cv2_centers = [(cv2_bins[i] + cv2_bins[i+1])/2 for i in range(len(cv2_bins)-1)]
    cv3_bins = np.linspace(0, 13, n_bins)
    cv3_centers = [(cv3_bins[i] + cv3_bins[i+1])/2 for i in range(len(cv3_bins)-1)]


    cv1_hist_sum = np.zeros([n_bins - 1])
    cv2_hist_sum = np.zeros([n_bins - 1])
    cv3_hist_sum = np.zeros([n_bins - 1])

    for i, walker_dir in enumerate(walker_dirs):
        # Extract CVs and unbiasing weights of individual simulation
        reweight_file = os.path.join(walker_dir, "COLVARS_REWEIGHT")
        data = read_plumed_data_file(reweight_file)
        
        # Calculate biasing weights from each BEMD replica


        # Generate histograms for each CV

        if not "WALKER0" in walker_dir:
            v_rescale = np.array([ v - np.max(data.values[:, -1]) for v in data.values[:, -1]])
            ws = np.exp(v_rescale / kt)

            # Generate histograms weighted by unbiasing weights
            cv1_hist, bin_edges = np.histogram(data.cv1, weights=ws, bins = cv1_bins, density = True)
            cv2_hist, bin_edges = np.histogram(data.cv2, weights=ws, bins = cv2_bins, density = True)
            cv3_hist, bin_edges = np.histogram(data.cv3, weights=ws, bins = cv3_bins, density = True)
        else:
            cv1_hist, bin_edges = np.histogram(data.cv1, bins = cv1_bins, density = True)
            cv2_hist, bin_edges = np.histogram(data.cv2, bins = cv2_bins, density = True)
            cv3_hist, bin_edges = np.histogram(data.cv3, bins = cv3_bins, density = True)
            


        cv1_hist_sum = cv1_hist_sum + cv1_hist
        cv2_hist_sum = cv2_hist_sum + cv2_hist
        cv3_hist_sum = cv3_hist_sum + cv3_hist

    cv1_hist_sum = cv1_hist_sum / 4
    cv2_hist_sum = cv2_hist_sum / 4
    cv3_hist_sum = cv3_hist_sum / 4

    # Plot aggregate probability distributions
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])

    axes[0,0].plot(cv1_centers, cv1_hist_sum)
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("Probability Density")
    axes[0,0].grid(visible=True, which="both", axis="both")


    axes[0,1].plot(cv2_centers, cv2_hist_sum)
    axes[0,1].set_xlabel("$R_{g}$")
    axes[0,1].set_ylabel("Probability Density")
    axes[0,1].grid(visible=True, which="both", axis="both")


    axes[1,0].plot(cv3_centers, cv3_hist_sum)
    axes[1,0].set_xlabel("$d_{end-to-end}$")
    axes[1,0].set_ylabel("Probability Density")
    axes[1,0].grid(visible=True, which="both", axis="both")


    plt.delaxes(axes[1,1])
    fig.suptitle("CV Aggregate Probabilities")
    plt.savefig("1D_FE_plots/CV_reweight_proabality_dist.png")
    plt.close()

    # Plot aggregate free energy surfaces    
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])
    
    fe_cv1 = -kt * np.log(cv1_hist_sum)
    axes[0,0].plot(cv1_centers, fe_cv1)
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("Free Energy (kJ/mol)")
    axes[0,0].grid(visible=True, which="both", axis="both")


    fe_cv2 = -kt * np.log(cv2_hist_sum)
    axes[0,1].plot(cv2_centers, fe_cv2)
    axes[0,1].set_xlabel("$R_{g}$")
    axes[0,1].set_ylabel("Free Energy (kJ/mol)")
    axes[0,1].grid(visible=True, which="both", axis="both")


    fe_cv3 = -kt * np.log(cv3_hist_sum)
    axes[1,0].plot(cv3_centers, fe_cv3)
    axes[1,0].set_xlabel("$d_{end-to-end}$")
    axes[1,0].set_ylabel("Free Energy (kJ/mol)")
    axes[1,0].grid(visible=True, which="both", axis="both")


    plt.delaxes(axes[1,1])
    fig.suptitle("CV Free Energy Surfaces")
    plt.savefig("1D_FE_plots/CV_free_energy_surf.png")
    plt.close()




@FlowProject.pre.isfile("WALKER0/COLVARS_REWEIGHT")
@FlowProject.operation
def plot_2D_FE_surface(job):
    # Need to modify HILLS filename to HILL.x
    os.chdir(job.fn(""))
    ts.utils.make_path("2D_FE_plots")
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))

    kt = 300 * 8.314462618 * 10 ** -3
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3

    if not os.path.exists(job.fn("WALKER0/COLVARS_REWEIGHT")):
        reweight_walker_trajectories(job, "plumed_multi_cv_reweight.dat", kt)

    cv1_cv2_hist = None
    cv2_cv3_hist = None
    cv1_cv3_hist = None

    n_bins = 15

    cv1_bins = np.linspace(0, 7,n_bins)
    cv1_centers = [(cv1_bins[i] + cv1_bins[i+1])/2 for i in range(len(cv1_bins)-1)]
    cv2_bins = np.linspace(0, 4, n_bins)
    cv2_centers = [(cv2_bins[i] + cv2_bins[i+1])/2 for i in range(len(cv2_bins)-1)]
    cv3_bins = np.linspace(0,13,n_bins)
    cv3_centers = [(cv3_bins[i] + cv3_bins[i+1])/2 for i in range(len(cv3_bins)-1)]

    cv1_cv2_hist = np.zeros([n_bins - 1, n_bins - 1])
    cv1_cv3_hist = np.zeros([n_bins - 1, n_bins - 1])
    cv2_cv3_hist = np.zeros([n_bins - 1, n_bins - 1])


    for i, walker_dir in enumerate(walker_dirs):
        # Extract CVs of individual simulation
        reweight_file = os.path.join(walker_dir, "COLVARS_REWEIGHT")
        data = read_plumed_data_file(reweight_file)
        
        # Calculate weights of each simulation
        print(walker_dir.split("/")[-1])
        ln_ws = np.array([ v / np.sum(data.values[:,-1]) for v in data.values[:, -1]])
        ws = np.exp(ln_ws)
    
        # Generate histogram for each combination of points
        hist1, xedges, yedges = np.histogram2d(data.cv1, data.cv2, weights=ws, bins = [cv1_bins, cv2_bins], density = True)
        hist2, xedges, yedges = np.histogram2d(data.cv1, data.cv3, weights=ws, bins = [cv1_bins, cv3_bins], density = True)
        hist3, xedges, yedges = np.histogram2d(data.cv2, data.cv3, weights=ws, bins = [cv2_bins, cv3_bins], density = True)

        if  "WALKER0" in walker_dir:
            hist1, xedges, yedges = np.histogram2d(data.cv1, data.cv2, bins = [cv1_bins, cv2_bins], density = True)
            hist2, xedges, yedges = np.histogram2d(data.cv1, data.cv3, bins = [cv1_bins, cv3_bins], density = True)
            hist3, xedges, yedges = np.histogram2d(data.cv2, data.cv3, bins = [cv2_bins, cv3_bins], density = True)
        
        # Plot individual FE surfaces

        fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])

        cv1_x, cv2_y = np.meshgrid(cv1_centers, cv2_centers)
        fe_cv1_cv2 = -kt * np.log(hist1.T)
        axes[0,0].contourf(cv1_x, cv2_y, fe_cv1_cv2)
        axes[0,0].set_xlabel("$N_{H-bonds}$")
        axes[0,0].set_ylabel("$R_{G}$")

        cv1_x, cv3_y = np.meshgrid(cv1_centers, cv3_centers)
        fe_cv1_cv3 = -kt * np.log(hist2.T)
        axes[0,1].contourf(cv1_x, cv3_y, fe_cv1_cv3)
        axes[0,1].set_xlabel("$N_{H-bonds}$")
        axes[0,1].set_ylabel("$d_{end-to-end}$")

        cv2_x, cv3_y = np.meshgrid(cv2_centers, cv3_centers)
        fe_cv2_cv3 = -kt * np.log(hist3.T)
        axes[1,0].contourf(cv2_x, cv3_y, fe_cv2_cv3)
        axes[1,0].set_xlabel("$R_{G}$")
        axes[1,0].set_ylabel("$d_{end-to-end}$")

        plt.delaxes(axes[1,1])
        fig.suptitle(walker_dir.split("/")[-1])
        plt.savefig(job.fn("2D_FE_plots/2D_FE_sim_" + str(i) + ".png"))
        plt.close()

        # Aggregate histograms
        cv1_cv2_hist = cv1_cv2_hist + hist1
        cv1_cv3_hist = cv1_cv3_hist + hist2
        cv2_cv3_hist = cv2_cv3_hist + hist3

    cv1_cv2_hist = cv1_cv2_hist / 4
    cv1_cv3_hist = cv1_cv3_hist / 4
    cv2_cv3_hist = cv2_cv3_hist / 4

    # Plot averaged FE surface
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])
    
    # Need to think about how to do this math...
    # Maybe ask Wei-Tse
    # FE calculation:
    # Averaging FE surface from for all 4 simulations
    cv1_x, cv2_y = np.meshgrid(cv1_centers, cv2_centers)
    fe_cv1_cv2 = -kt * np.log(cv1_cv2_hist.T)

    # Plot contours
    plot = axes[0,0].contourf(cv1_x, cv2_y, cv1_cv2_hist.T)
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("$R_{G}$")
    cbar = plt.colorbar(plot, ax = axes[0,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")

    # CV1 CV2 FE
    cv1_x, cv3_y = np.meshgrid(cv1_centers, cv3_centers)
    fe_cv1_cv3 = -kt * np.log(cv1_cv3_hist.T)

    # Plot contours
    plot = axes[0,1].contourf(cv1_x, cv3_y, cv1_cv3_hist.T)
    axes[0,1].set_xlabel("$N_{H-bonds}$")
    axes[0,1].set_ylabel("$d_{end-to-end}$")
    cbar = plt.colorbar(plot, ax = axes[1,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")


    # CV2 CV3 FE
    cv2_x, cv3_y = np.meshgrid(cv2_centers, cv3_centers)
    fe_cv2_cv3 = -kt * np.log(cv2_cv3_hist.T)

    plot = axes[1,0].contourf(cv2_x, cv3_y, cv2_cv3_hist.T)
    axes[1,0].set_xlabel("$R_{G}$")
    axes[1,0].set_ylabel("$d_{end-to-end}$")
    cbar = plt.colorbar(plot, ax = axes[0,1])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")


    plt.delaxes(axes[1,1])
    fig.suptitle("2D Probability Distributions")

    plt.savefig(job.fn("2D_FE_plots/2D_CV_hist.png"))
    plt.close()

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])
    
    # Plot contours
    plot = axes[0,0].contourf(cv1_x, cv2_y, fe_cv1_cv2)
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("$R_{G}$")
    cbar = plt.colorbar(plot, ax = axes[0,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")

    plot = axes[0,1].contourf(cv1_x, cv3_y, fe_cv1_cv3)
    axes[0,1].set_xlabel("$N_{H-bonds}$")
    axes[0,1].set_ylabel("$d_{end-to-end}$")
    cbar = plt.colorbar(plot, ax = axes[1,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")



    plot = axes[1,0].contourf(cv2_x, cv3_y, fe_cv2_cv3)
    axes[1,0].set_xlabel("$R_{G}$")
    axes[1,0].set_ylabel("$d_{end-to-end}$")
    cbar = plt.colorbar(plot, ax = axes[0,1])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")


    plt.delaxes(axes[1,1])
    fig.suptitle("2D Free Energy Surfaces")
    plt.savefig(job.fn("2D_FE_plots/2D_FE_surfaces.png"))





def main():
    FlowProject().main()
