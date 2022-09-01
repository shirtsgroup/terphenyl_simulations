import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tqdm import tqdm
from .observables import get_torsions

def plot_grid_search(metric_matrix, x_ticks, y_ticks, x_label, y_label, filename):
    # combination metric figure
    fig = plt.figure(dpi = 300)
    ax = fig.add_subplot(111)
    cax = ax.matshow(metric_matrix)
    fig.colorbar(cax)
    xaxis = np.arange(len(x_ticks))
    yaxis = np.arange(len(y_ticks))
    ax.set_xticks(xaxis)
    ax.set_yticks(yaxis)
    ax.set_xticklabels(x_ticks)
    ax.set_yticklabels(y_ticks)
    ax.set_xlabel(x_label)
    ax.set_ylabel(x_label)
    plt.savefig(filename)

def plot_neighbor_dependent_2D_histogram(traj_obj, torsions_atom_names_pair, torsion_ids, prefix, n_bins = 50, legend = None, mirror_sym = False):
    # Get 2D distriubtion from pairs
    torsions = [get_torsions(traj_obj, torsions_atom_names_pair[i]) for i in range(len(torsions_atom_names_pair)) ]

    # Bins and centers for distributuion
    bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    if mirror_sym:
        bin_edges = np.linspace(0, np.pi, n_bins + 1)
    bin_centers = np.array([(bin_edges[i] + bin_edges[i+1]) * 0.5 for i in range(len(bin_edges) - 1) ])
    
    # Figure setup
    n_rows = len(torsions_atom_names_pair)
    fig, axes = plt.subplots(nrows = n_rows, ncols = n_rows, dpi = 300, figsize = [2 * n_rows, 2 * n_rows])

    print("Plotting probability distributions...")
    for i in range(len(tqdm(torsions_atom_names_pair))):
        for j in range(len(torsions_atom_names_pair)):
            if i == j:
                # 1D distribution
                hist, bin_edges_out = np.histogram(np.array(torsions[i]), bins = bin_edges, density = True)
                axes[i,j].plot(bin_centers, hist)
                axes[i,j].set_xlabel(torsion_ids[i])
                axes[i,j].set_ylabel("Density")
            if i < j:
                # 2D distribution
                h = axes[i,j].hist2d(torsions[j], torsions[i], bins = bin_edges, density = True, cmap = plt.cm.viridis)
                axes[i,j].set_xlabel(torsion_ids[j])
                axes[i,j].set_ylabel(torsion_ids[i])
                # plt.colorbar(cax = )

            if i > j:
                axes[i,j].remove()

    plt.tight_layout()
    plt.savefig(prefix + "_2d_torsions.png")

def plot_torsions_distributions(traj_obj_list, torsion_atom_names, torsion_id, prefix, n_bins = 50, legend = None, mirror_sym = False):
    # Setup figure
    plt.figure(dpi = 300)
    sns.set_palette("plasma", n_colors = len(traj_obj_list))

    # Bins and centers for distributuion
    bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    if mirror_sym:
        bin_edges = np.linspace(0, np.pi, n_bins + 1)
    bin_centers = np.array([(bin_edges[i] + bin_edges[i+1]) * 0.5 for i in range(len(bin_edges) - 1) ])

    # Get torsions, bin and plot
    for traj_obj in tqdm(traj_obj_list):
        torsions = get_torsions(traj_obj, torsion_atom_names, mirror_sym = mirror_sym)
        hist, bin_edges_out = np.histogram(np.array(torsions), bins = bin_edges, density = True)
        plt.plot(bin_centers, hist)
        plt.xlabel(torsion_id)
        plt.ylabel("Density")
    if legend is not None:
        plt.legend(legend)
    plt.savefig(prefix + "_torsions.png")