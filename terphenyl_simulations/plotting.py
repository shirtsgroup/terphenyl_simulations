import matplotlib.pyplot as plt
import matplotlib
from .utils import make_path, GromacsLogFile
import mdtraj as md
import os
import numpy as np
import seaborn as sns
from tqdm import tqdm
from .observables import get_torsions
import itertools
import sys


def plot_grid_search(metric_matrix, x_ticks, y_ticks, x_label, y_label, filename, cbar_label):
    """
    Function for plotting metrics against a grid search of parameters
    """

    # Create a directory if prefix has "/" in it
    if "/" in filename:
        make_path(filename)

    # combination metric figure
    fig = plt.figure(dpi=300, figsize = [10, 10])
    ax = fig.add_subplot(111)
    cax = ax.matshow(metric_matrix)
    fig.colorbar(cax).set_label( label = cbar_label, size=18)
    xaxis = np.arange(len(x_ticks))
    yaxis = np.arange(len(y_ticks))
    ax.set_xticks(xaxis)
    ax.set_yticks(yaxis)
    ax.set_xticklabels(x_ticks, rotation = 90)
    ax.set_yticklabels(y_ticks)
    ax.set_xlabel(x_label, fontsize = 18)
    ax.set_ylabel(y_label, fontsize = 18)
    plt.savefig(filename)
    plt.close()


def plot_neighbor_dependent_2D_histogram(
    traj_obj,
    torsions_atom_names_pair,
    torsion_ids,
    prefix,
    n_bins=50,
    legend=None,
    mirror_sym=False,
):

    # Create a directory if prefix has "/" in it
    if "/" in prefix:
        make_path(prefix)

    # Get 2D distriubtion from pairs
    torsions = [
        get_torsions(traj_obj, torsions_atom_names_pair[i])
        for i in range(len(torsions_atom_names_pair))
    ]

    # Bins and centers for distributuion
    bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    if mirror_sym:
        bin_edges = np.linspace(0, np.pi, n_bins + 1)
    bin_centers = np.array(
        [(bin_edges[i] + bin_edges[i + 1]) * 0.5 for i in range(len(bin_edges) - 1)]
    )

    # Figure setup
    n_rows = len(torsions_atom_names_pair)
    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_rows, dpi=300, figsize=[2 * n_rows, 2 * n_rows]
    )

    for i in range(len(torsions_atom_names_pair)):
        for j in range(len(torsions_atom_names_pair)):
            if i == j:
                # 1D distribution
                hist, bin_edges_out = np.histogram(
                    np.array(torsions[i]), bins=bin_edges, density=True
                )
                axes[i, j].plot(bin_centers, hist)
                axes[i, j].set_xlabel(torsion_ids[i])
                axes[i, j].set_ylabel("Density")
            if i < j:
                # 2D distribution
                h = axes[i, j].hist2d(
                    torsions[j],
                    torsions[i],
                    bins=bin_edges,
                    density=True,
                    cmap=plt.cm.viridis,
                )
                axes[i, j].set_xlabel(torsion_ids[j])
                axes[i, j].set_ylabel(torsion_ids[i])
                # plt.colorbar(cax = )

            if i > j:
                axes[i, j].remove()

    plt.tight_layout()
    plt.savefig(prefix + "_2d_torsions.png")
    plt.close()


def plot_torsions_distributions(
    traj_obj_list,
    torsion_atom_names,
    x_axis,
    prefix,
    title,
    n_bins=50,
    legend=None,
    mirror_sym=False,
    offsets=None,
    figsize=None
):
    """
    Function for plotting 1D torsion distributions using MDTraj objects

    Parameters
    ----------
    traj_obj_list : list
        List of mdtraj.trajectory objects to plot torsions for
    torsion_atom_names : list
        List of lists of strings of atom name strings used to define atoms to use when
        getting torsions from trajectory. If multiple list of lists are provided,
        the function will iterate over each list of lists for a given multiple provided
        trajectory files.
    x_axis : string
        Label for x axis. y axis defaults to Density
    prefix : string
        Prefix for file name. Can include directories.
    title: str
        Plot title
    n_bins : int
        Number of bins to use across total angle space
    legend : list
        List of legend entries for plot. No legend plotted if left as None
    mirror_sym : bool
        Collapse 360 degrees histogram to 180 degree histogram. Useful for symetric torsions
    offsets : list
        If any distributions need to be shifted by a number of degrees, providing a
        list of offset values will apply the offset a given dataset based on index.
    """

    if type(traj_obj_list) == md.Trajectory:
        traj_obj_list = [traj_obj_list]

    # Setup figure
    plt.figure(dpi=300)
    sns.set_palette("plasma", n_colors=len(traj_obj_list))

    # Bins and centers for distributuion
    bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    if mirror_sym:
        bin_edges = np.linspace(0, np.pi, n_bins + 1)
    bin_centers = np.array(
        [(bin_edges[i] + bin_edges[i + 1]) * 0.5 for i in range(len(bin_edges) - 1)]
    )

    # Get torsions, bin and plot
    for i, traj_obj in enumerate(traj_obj_list):
        if type(torsion_atom_names[0]) is str:
            torsions = get_torsions(traj_obj, [torsion_atom_names], mirror_sym=mirror_sym)
        if type(torsion_atom_names[0]) is list:
            if type(torsion_atom_names[0][0]) is str:
                torsions = get_torsions(
                    traj_obj, torsion_atom_names, mirror_sym=mirror_sym
                )
            if type(torsion_atom_names[0][0]) is list:
                torsions = get_torsions(
                    traj_obj, torsion_atom_names[i], mirror_sym=mirror_sym
                )        
        if offsets is not None:
            torsions += offsets[i]
        hist, bin_edges_out = np.histogram(
            np.array(torsions), bins=bin_edges, density=True
        )
        plt.figure(figsize = figsize)
        plt.plot(bin_centers, hist)
        plt.xlabel(x_axis)
        plt.ylabel("Density")
        plt.title(title)
    if legend is not None:
        plt.legend(legend)
    
    plt.savefig(prefix + "_torsions.png")
    plt.close()


def plot_torsion_timeseries(traj_obj, torsion_atom_names, filenames, titles = None, time_per_frame = 0.2, output_dir = "torsion_timeseries", degrees = False):
    """
    Function for plotting time series of a set of torsions

    Parameters
    ----------
    traj_obj : mdtraj.Trajectory
        mdtraj trajectory object to pull torsions from
    torsion_atom_names : list
        list of torsions to create time series for
    filename : str
        filename of output plot
    time_per_frame : float
        number of ns per frame
    titles : list (None)
        plot titles for each individual torsion timeseries
    output_dir : str (torsion_timeseries)
        string of directory to output all plot files
    degrees : bool (False)
        Whether to plot torsion in degrees or radians
    """
    time = np.array(list(range(traj_obj.n_frames))) * time_per_frame
    if titles is None:
        titles = ["Plot " + str(i) for i in range(len(torsion_atom_names))]

    make_path(output_dir)

    for i, torsion_id in enumerate(torsion_atom_names):
        torsions = get_torsions(traj_obj, [torsion_id])
        if degrees is True:
            torsions *= 180 / np.pi
        plt.figure(dpi = 300)
        plt.scatter(time, torsions, s=5)
        plt.xlabel("Time (ns)")
        if degrees is True:
            plt.ylabel("Torsion (Degrees)")
            plt.ylim([-180, 180])
        else:
            plt.ylabel("Torsion (Radians)")
            plt.ylim([-np.pi, np.pi])
        plt.title(titles[i])
        plt.savefig(os.path.join(output_dir, filenames[i]))
        plt.close()

def plot_ramachandran_plot(traj_file, top_file, prefix = "remd", bins = 50, title = None, scatter_points_files = [], legend = None):
    """
    Create a ramachandran plot for a peptide system
    """

    traj_object = md.load_xtc(traj_file, top = top_file)

    
    phi_inds, phi_angles = md.compute_phi(traj_object)
    psi_inds, psi_angles = md.compute_psi(traj_object)

    phi_angles = phi_angles.flatten()
    psi_angles = psi_angles.flatten()

    # Get temperature from mdp file
    
    plt.figure(figsize=[5,5])
    plt.hist2d(phi_angles, psi_angles, bins = np.linspace(-np.pi, np.pi, bins + 1))
    plt.xlabel("$\Phi$ Angle (Radians)")
    plt.ylabel("$\Psi$ Angle (Radians)")

    for scatter_points_file in scatter_points_files:
        frame = md.load(scatter_points_file)
        phi_i, phi_angle = md.compute_phi(frame)
        psi_i, psi_angle = md.compute_psi(frame)
        plt.scatter(phi_angle.flatten(), psi_angle.flatten(), marker="x")
    if legend is None:
        plt.legend([fn.split("/")[-1].split(".")[0] for fn in scatter_points_files])
    else:
        plt.legend(legend)

    if title is not None:
        plt.title(title)

    plt.savefig(prefix + "_ramachandran.png")
    plt.close()

def plot_bemd_state_index_plot(log_file, prefix, stride = 500):
    log_file_obj = GromacsLogFile(log_file)
    fig, axes = plt.subplots(
        int(np.ceil(np.sqrt(log_file_obj.n_states))),
        int(np.ceil(np.sqrt(log_file_obj.n_states))),
        figsize = [10,10]
    )
 
    ax_indices = list(itertools.product(range(int(np.ceil(np.sqrt(log_file_obj.n_states)))), repeat = 2))

    for i in range(log_file_obj.n_states):
        state_traj = [state_i[0] for state_i in log_file_obj.states]
        axes[ax_indices[i][0], ax_indices[i][1]].plot(list(range(len(state_traj)))[::stride], state_traj[::stride], linewidth = 0.5)
        axes[ax_indices[i][0], ax_indices[i][1]].set_xlabel("Steps")
        axes[ax_indices[i][0], ax_indices[i][1]].set_ylabel("State Index")
        axes[ax_indices[i][0], ax_indices[i][1]].set_title("Replica " + str(i))
    
    plt.tight_layout()
    plt.savefig(prefix + "_state_index.png", dpi = 150)

def plot_bemd_transition_matrix(log_file, prefix):
    log_file_obj = GromacsLogFile(log_file)
    fig, ax = plt.subplots()
    ax.matshow(log_file_obj.transition_matrix, cmap="YlGn")
    for (i, j), z in np.ndenumerate(log_file_obj.transition_matrix):
        ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center')
    plt.savefig(prefix + "_transition_matrix.png", dpi = 150)



