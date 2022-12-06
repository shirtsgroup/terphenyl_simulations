import heteropolymer_simulations as hs
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from itertools import repeat
from sklearn import metrics
import time


def run_clustering(rmsd_matrix, eps, ms, parallel):
    dbscan, labels = hs.clustering.DBSCAN_RMSD_clustering(rmsd_matrix, eps, ms, parallel = parallel, write_to_file = False)
    if len(np.unique(labels)) > 1:
        ss = metrics.silhouette_score(rmsd_matrix, labels)
    else:
        ss = np.nan

    return (len(np.unique(labels)), ss)


def main():
    t1 = time.time()
    n_min_samples = 40
    n_eps = 40
    n_proccesses = 64
    total_cores = 64


    # Load trajectory
    traj = md.load("sim0/npt_new.whole.xtc", top = "sim0/berendsen_npt.gro")
    
    # Remove solvent
    top = traj.topology
    selection = top.select("resname OCT")
    oct_traj = traj.atom_slice(selection)

    # Construct RMSD matrix once
    rmsd_matrix = hs.clustering.construct_rmsd_matrix(oct_traj)

    # Grid search for DBSCAN hyperparameters
    max_rmsd = np.max(rmsd_matrix)
    eps_values = np.linspace(0.01 * max_rmsd, max_rmsd/4, n_eps)
    total_frames = oct_traj.n_frames
    min_sample_values = np.linspace(0.01 * total_frames, 0.15 * total_frames, n_min_samples, dtype = int)
    min_samples_values = [int(a) for a in min_sample_values]
    grid_shape = (len(eps_values), len(min_sample_values))
    X, Y = np.meshgrid(min_sample_values, eps_values)
    
    # Use pool to fill array
    with mp.Pool(processes = n_proccesses) as pool:
        outputs = pool.starmap(run_clustering, zip(repeat(rmsd_matrix), Y.reshape(-1), X.reshape(-1), repeat(int(total_cores/n_proccesses))))
        n_clusters = [a[0] for a in outputs]
        ss = [a[1] for a in outputs]
        n_clusters = np.array(n_clusters).reshape(grid_shape)
        ss = np.array(ss).reshape(grid_shape)


    # SS figure
    fig = plt.figure(dpi = 300)
    ax = fig.add_subplot(111)
    cax = ax.matshow(ss)
    fig.colorbar(cax)
    xaxis = np.arange(len(min_samples_values))
    yaxis = np.arange(len(eps_values))
    ax.set_xticks(xaxis)
    ax.set_yticks(yaxis)
    ax.set_xticklabels(min_samples_values)
    ax.set_yticklabels([round(eps, 2) for eps in eps_values])
    ax.set_xlabel("Min. Samples")
    ax.set_ylabel('$\epsilon_{DBSCAN}$')
    plt.savefig("grid_search_silhouette_scores.png")

    # SS figure
    fig = plt.figure(dpi = 300)
    ax = fig.add_subplot(111)
    cax = ax.matshow(n_clusters)
    fig.colorbar(cax)
    xaxis = np.arange(len(min_samples_values))
    yaxis = np.arange(len(eps_values))
    ax.set_xticks(xaxis)
    ax.set_yticks(yaxis)
    ax.set_xticklabels(min_samples_values)
    ax.set_yticklabels([round(eps, 2) for eps in eps_values])
    ax.set_xlabel("Min. Samples")
    ax.set_ylabel('$\epsilon_{DBSCAN}$')
    plt.savefig("grid_search_n_clusters.png")

    t2 = time.time()

    print("DBSCAN grid search took:", round(t2 - t1, 2), "seconds.")

        
    


if __name__ == "__main__":
    main()