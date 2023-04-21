import matplotlib.pyplot as plt
from .utils import make_path
import mdtraj as md
import os
import numpy as np
import seaborn as sns
from tqdm import tqdm
from .observables import get_torsions


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
    fig.colorbar(cax, label = cbar_label)
    xaxis = np.arange(len(x_ticks))
    yaxis = np.arange(len(y_ticks))
    ax.set_xticks(xaxis)
    ax.set_yticks(yaxis)
    ax.set_xticklabels(x_ticks, rotation = 90)
    ax.set_yticklabels(y_ticks)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
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
        plt.plot(bin_centers, hist)
        plt.xlabel(x_axis)
        plt.ylabel("Density")
        plt.title(title)
    if legend is not None:
        plt.legend(legend)
    make_path(prefix + "_torsions.png")
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
