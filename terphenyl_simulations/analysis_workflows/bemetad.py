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
import yaml
from .utils import *
from .labels import *

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
    params = read_yaml_parameter_file("analysis_parameters.yml")
    stride = params["CV_plots"]["CV_traj_stride"]
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))
    fig, axs = plt.subplots(nrows = len(walker_dirs), ncols = len(walker_dirs)+1, figsize = [20,10])
    
    overall_cvs = [[] for _ in range(len(params["CV_plots"]["x_labels"]))]
    
    for i, walker_dir in tqdm(enumerate(walker_dirs)):
        colvars_file = os.path.join(walker_dir, "COLVARS."+str(i))
        colvars = read_plumed_data_file(colvars_file)
        axs[0,i].set_title("Simulation " + str(i))
        for j in range(len(params["CV_plots"]["x_labels"])):
            cv_label = params["CV_plots"]["cv_names"][j]
            overall_cvs[j] += list(colvars[cv_label][::stride].values)
            axs[j,i].plot(1/1000 * colvars["time"][::stride], colvars[cv_label][::stride], markersize = 0.5)
            axs[j,i].set_xlabel("Time (ns)")
            axs[j,i].set_ylabel(params["CV_plots"]["x_labels"][j])
        axs[-1, i].plot(1/1000 * colvars["time"][::stride], colvars.values[:,-1][::stride], markersize = 0.5, color = "red")
        axs[-1, i].set_xlabel("Time (ns)")
        axs[-1, i].set_ylabel("Bias Energy (kJ/mol)")

    for j in range(len(params["CV_plots"]["x_labels"])):
        axs[j, -1].hist(overall_cvs[j], bins = 50, density = True)
        axs[j, -1].set_xlabel(params["CV_plots"]["x_labels"][j])
        axs[j, -1].set_ylabel("Probability Density")
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
    params = read_yaml_parameter_file("analysis_parameters.yml")
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))
    n_simulations = len(walker_dirs)
    kt = 300 * 8.314462618 * 10 ** -3
    T = 300
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3
        T = job.sp["TEMP"]
    x_labels = params["sum_hills_plots"]["x_labels"]
    xy_limits = params["sum_hills_plots"]["fe_xy_limits"]
    plumed_bins = params["sum_hills_plots"]["plumed_bins"]
    cv_bfs = params["sum_hills_plots"]["cv_bf_labels"]

    # Setup subplots
    n_row_columns = int(np.ceil(np.sqrt(n_simulations - 1)))
    fig, axes = plt.subplots(nrows = n_row_columns, ncols = n_row_columns, figsize = [10, 10])
    axes = axes.flatten().tolist()

    for i, (walker_dir, ax) in enumerate(zip(walker_dirs[1:], axes)):

        # Read HILLS file
        os.chdir(walker_dir)
        walker_id = int(walker_dir.split("WALKER")[-1])
        subprocess.run(["plumed", "--no-mpi", "sum_hills", "--hills", "HILLS." + str(walker_id), "--kt", str(kt), "--mintozero", "--bin", str(plumed_bins[i][2]), "--min", str(plumed_bins[i][0]), "--max", str(plumed_bins[i][1])]) # Silences output from sum_hills
        filename = "fes.dat"
        fes_data = read_plumed_fes_file(filename)

        # Calculating correction scaling for FE surface
        delta_T = job.sp[cv_bfs[i]] * T - T
        prefactor = (T + delta_T) / delta_T

        print("T + dT / dT =", prefactor)

        # Plot FES
        ax.plot(prefactor * fes_data.values[:,0], fes_data.values[:,1])
        ax.set_xlabel(x_labels[i])
        ax.set_ylabel("Free Energy (kJ/mol)")
        ax.set_xlim(xy_limits[i][0])
        ax.set_ylim(xy_limits[i][1])
        ax.set_title("Simulation " + str(walker_id))
        ax.grid(visible=True, which="both", axis="both")
    os.chdir(current_dir)
    fig.suptitle(job.sp, wrap=True)
    plt.savefig(job.fn("multi_sum_hills_FE.png"), transparent = False, dpi = 300)

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("multi_sum_hills_FE_stride.png")
@FlowProject.operation
def sum_hills_FE_stride(job):
    params = read_yaml_parameter_file("analysis_parameters.yml")

    # Get original directory
    current_dir = os.path.abspath("")
    
    # Read in analysis parameters
    walker_dirs = natsorted(glob.glob(job.fn("WALKER*")))
    n_simulations = len(walker_dirs)
    kt = 300 * 8.314462618 * 10 ** -3
    T = 300
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3
        T = job.sp["TEMP"]
    x_labels = params["sum_hills_plots"]["x_labels"]
    xy_limits = params["sum_hills_plots"]["fe_xy_limits"]
    plumed_bins = params["sum_hills_plots"]["plumed_bins"]
    cv_bfs = params["sum_hills_plots"]["cv_bf_labels"]

    # Setup subplots
    n_row_columns = int(np.ceil(np.sqrt(n_simulations - 1)))
    fig, axes = plt.subplots(nrows = n_row_columns, ncols = n_row_columns, figsize = [10, 10])
    axes = axes.flatten().tolist()

    for i, (walker_dir, ax) in enumerate(zip(walker_dirs[1:], axes)):
        stride = 4000

        # Read HILLS file
        os.chdir(walker_dir)
        walker_id = int(walker_dir.split("WALKER")[-1])
        subprocess.run(["plumed", "--no-mpi", "sum_hills", "--hills", "HILLS." + str(walker_id), "--kt", str(kt), "--bin", str(plumed_bins[i][2]), "--min", str(plumed_bins[i][0]), "--max", str(plumed_bins[i][1]), "--stride", str(stride)]) # Silences output from sum_hills
        fes_files = natsorted(glob.glob("fes*"))
        cmap = plt.get_cmap("viridis")
        colors = [cmap(i) for i in np.linspace(0, 1, len(fes_files))]
        for j, fes_file in enumerate(fes_files):
            fes_data = read_plumed_fes_file(fes_file)

            # Calculating correction scaling for FE surface
            delta_T = job.sp[cv_bfs[i]] * T - T
            prefactor = (T + delta_T) / delta_T

            # Plot FES
            ax.plot(prefactor * fes_data.values[:,0], fes_data.values[:,1], lw = 1, color = colors[j])
        ax.set_xlabel(x_labels[i])
        ax.set_ylabel("Free Energy (kJ/mol)")
        ax.set_xlim(xy_limits[i][0])
        # ax.set_ylim(xy_limits[i][1])
        ax.set_title("Simulation " + str(walker_id))
        ax.grid(visible=True, which="both", axis="both")
    os.chdir(current_dir)
    fig.suptitle(job.sp, wrap = True)
    plt.tight_layout()
    plt.savefig(job.fn("multi_sum_hills_FE_stride.png"), transparent = False, dpi = 300)


@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("WALKER3/COLVARS_REWEIGHT")
@FlowProject.operation
def unbias_simulations(job):

    kt = 300 * 8.314462618 * 10 ** -3
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3

    reweight_colvars(job, "plumed_multi_cv_reweight.dat")


@FlowProject.pre.isfile("WALKER0/COLVARS_REWEIGHT")
# @FlowProject.post.isfile("1D_FE_plots/CV_free_energy_surf.png")
@FlowProject.operation
def plot_1D_FE_surface(job):
    # Read parameter file
    params = read_yaml_parameter_file("analysis_parameters.yml")
    current_dir = os.path.abspath("")


    # Navigate to working dir
    os.chdir(job.fn(""))
    ts.utils.make_path("1D_FE_plots")
    walker_dirs = natsorted(glob.glob("WALKER*"))


    # Calculate kT for specific temperatures
    kt = 300 * 8.314462618 * 10 ** -3
    if "TEMP" in job.sp.keys():
        kt = job.sp["TEMP"] * 8.314462618 * 10 ** -3

    # Generate common bins for all CVs

    cv_bins = []
    cv_centers = []

    for n_bins, cv_range in zip(params["1D_FE_plots"]["n_bins"], params["1D_FE_plots"]["cv_ranges"]):
        bins = np.linspace(cv_range[0], cv_range[1], n_bins)
        cv_bins.append(bins)
        cv_centers.append([(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)])

    cv_hist_sums = [np.zeros(len(cv_centers[i])) for i in range(len(cv_centers))]

    for i, walker_dir in enumerate(walker_dirs):
        # Extract CVs and unbiasing weights of individual simulation
        reweight_file = os.path.join(walker_dir, "COLVARS_REWEIGHT")
        data = read_plumed_data_file(reweight_file)
        
        # Calculate biasing weights from each BEMD replica
        # Generate histograms for each CV
        ws = None
        if not "WALKER0" in walker_dir:
            v_rescale = np.array([ v - np.max(data.values[:, -1]) for v in data.values[:, -1]])
            ws = np.exp(v_rescale / kt)

        # Generate histograms weighted by unbiasing weights
        for j in range(len(cv_bins)):
            cv_hist, bin_edges = np.histogram(getattr(data, "cv" +str(j + 1)), weights=ws, bins =cv_bins[j], density = True)
            cv_hist_sums[j] += cv_hist

    cv_hist_sums = [hist / len(walker_dirs) for hist in cv_hist_sums]

    # Plot aggregate probability distributions
    fig, axes = plt.subplots(nrows = int(np.ceil(np.sqrt(len(cv_bins)))), ncols = int(np.ceil(np.sqrt(len(cv_bins)))), figsize = [10, 10])
    axes = axes.flatten()

    # Plot 1D distributions
    for i in range(len(cv_bins)):
        axes[i].plot(cv_centers[i], cv_hist_sums[i])
        axes[i].set_xlabel(params["1D_FE_plots"]["x_labels"][i])
        axes[i].set_ylabel("Probability Density")
        axes[i].grid(visible=True, which="both", axis="both")
    
    # Remove unused plots
    for i in range(len(cv_bins), len(axes)):
        plt.delaxes(axes[i])

    plt.savefig("1D_FE_plots/CV_reweight_proabality_dist.png")
    plt.close()

    # Plot aggregate free energy surfaces   
    fig, axes = plt.subplots(nrows = int(np.ceil(np.sqrt(len(cv_bins)))), ncols = int(np.ceil(np.sqrt(len(cv_bins)))), figsize = [10, 10])
    axes = axes.flatten()

    # Plot 1D FE surfaces
    for i in range(len(cv_bins)):
        axes[i].plot(cv_centers[i], -kt * np.log(cv_hist_sums[i]))
        axes[i].set_xlabel(params["1D_FE_plots"]["x_labels"][i])
        axes[i].set_ylabel("Free Energy (kJ/mol)")
        axes[i].grid(visible=True, which="both", axis="both")
    
    # Remove unused plots
    for i in range(len(cv_bins), len(axes)):
        plt.delaxes(axes[i])

    plt.savefig("1D_FE_plots/CV_free_energy_surf.png")
    plt.close()

    os.chdir(current_dir)





@FlowProject.pre.isfile("WALKER0/COLVARS_REWEIGHT")
@FlowProject.post.isfile("2D_FE_plots/2D_FE_surfaces.png")
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

    n_bins = 20

    cv1_bins = np.linspace(0, 7,n_bins)
    cv1_centers = [(cv1_bins[i] + cv1_bins[i+1])/2 for i in range(len(cv1_bins)-1)]
    cv2_bins = np.linspace(0, 4, n_bins)
    cv2_centers = [(cv2_bins[i] + cv2_bins[i+1])/2 for i in range(len(cv2_bins)-1)]
    cv3_bins = np.linspace(0,13,n_bins)
    cv3_centers = [(cv3_bins[i] + cv3_bins[i+1])/2 for i in range(len(cv3_bins)-1)]

    cv1_cv2_hist = np.zeros([n_bins - 1, n_bins - 1])
    cv1_cv3_hist = np.zeros([n_bins - 1, n_bins - 1])
    cv2_cv3_hist = np.zeros([n_bins - 1, n_bins - 1])

    # Free energy levels to plot
    levels = np.linspace(-10, 20, 40)

    for i, walker_dir in enumerate(walker_dirs):
        # Extract CVs of individual simulation
        reweight_file = os.path.join(walker_dir, "COLVARS_REWEIGHT")
        data = read_plumed_data_file(reweight_file)
        
        # Calculate weights of each simulation
        print(walker_dir.split("/")[-1])
        v_rescale = np.array([ v - np.max(data.values[:, -1]) for v in data.values[:, -1]])
        ws = np.exp(v_rescale / kt)
    
        # Generate histogram for each combination of points
        hist1, xedges, yedges = np.histogram2d(data.cv1, data.cv2, weights=ws, bins = [cv1_bins, cv2_bins], density = True)
        hist2, xedges, yedges = np.histogram2d(data.cv1, data.cv3, weights=ws, bins = [cv1_bins, cv3_bins], density = True)
        hist3, xedges, yedges = np.histogram2d(data.cv2, data.cv3, weights=ws, bins = [cv2_bins, cv3_bins], density = True)

        if  not "WALKER0" in walker_dir:
            hist1, xedges, yedges = np.histogram2d(data.cv1, data.cv2, weights=ws, bins = [cv1_bins, cv2_bins], density = True)
            hist2, xedges, yedges = np.histogram2d(data.cv1, data.cv3, weights=ws, bins = [cv1_bins, cv3_bins], density = True)
            hist3, xedges, yedges = np.histogram2d(data.cv2, data.cv3, weights=ws, bins = [cv2_bins, cv3_bins], density = True)   
        else:
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

    # Plot probability distributions
    plot = axes[0,0].contourf(cv1_x, cv2_y, cv1_cv2_hist.T)
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("$R_{G}$")
    cbar = plt.colorbar(plot, ax = axes[0,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")

    # CV1 CV2 FE
    cv1_x, cv3_y = np.meshgrid(cv1_centers, cv3_centers)
    fe_cv1_cv3 = -kt * np.log(cv1_cv3_hist.T)

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
    plt.tight_layout()
    fig.suptitle(job.sp, wrap=True)

    plt.savefig(job.fn("2D_FE_plots/2D_CV_hist.png"))
    plt.close()

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])
    
    # Plot free energy surfaces
    plot = axes[0,0].contourf(cv1_x, cv2_y, fe_cv1_cv2, levels = levels)
    axes[0,0].set_xlabel("$N_{H-bonds}$")
    axes[0,0].set_ylabel("$R_{G}$")
    cbar = plt.colorbar(plot, ax = axes[0,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")

    plot = axes[0,1].contourf(cv1_x, cv3_y, fe_cv1_cv3, levels = levels)
    axes[0,1].set_xlabel("$N_{H-bonds}$")
    axes[0,1].set_ylabel("$d_{end-to-end}$")
    cbar = plt.colorbar(plot, ax = axes[1,0])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")



    plot = axes[1,0].contourf(cv2_x, cv3_y, fe_cv2_cv3, levels = levels)
    axes[1,0].set_xlabel("$R_{G}$")
    axes[1,0].set_ylabel("$d_{end-to-end}$")
    cbar = plt.colorbar(plot, ax = axes[0,1])
    cbar.ax.set_ylabel("Free Energy (kJ/mol)")


    plt.delaxes(axes[1,1])
    fig.suptitle(job.sp, wrap=True)
    plt.savefig(job.fn("2D_FE_plots/2D_FE_surfaces.png"))





def main():
    FlowProject().main()
