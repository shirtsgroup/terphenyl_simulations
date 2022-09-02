import heteropolymer_simulations as hs
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn import preprocessing 
import time
from DBCV import DBCV
from scipy.spatial.distance import euclidean
import sys
from tqdm import tqdm
plt.rcParams.update({'font.size': 10})

def main():
    t1 = time.time()
    n_min_samples = 30
    n_eps = 30
    n_processes = 32


    # Load trajectory
    traj = md.load(["sim0/npt_new.whole.xtc", "sim1/npt_new.whole.xtc", "sim3/npt_new.whole.xtc"], top = "sim0/berendsen_npt.gro")
    
    # Remove solvent
    top = traj.topology
    selection = top.select("resname TET")
    oct_traj = traj.atom_slice(selection)

    # Construct RMSD matrix once
    rmsd_matrix = hs.clustering.construct_rmsd_matrix(oct_traj)

    # Grid search for DBSCAN hyperparameters
    max_rmsd = np.max(rmsd_matrix)
    eps_values = np.linspace(0.01 * max_rmsd, max_rmsd / 2, n_eps)
    total_frames = oct_traj.n_frames
    min_sample_values = np.linspace(0.01 * total_frames, 0.25 * total_frames, n_min_samples, dtype = int)
    min_samples_values = [int(a) for a in min_sample_values]
    ss = np.zeros((len(eps_values), len(min_sample_values)))
    n_clusters = np.zeros((len(eps_values), len(min_sample_values)))
    # dbcvs = np.zeros((len(eps_values), len(min_sample_values)))
    for i, eps in enumerate(tqdm(eps_values)):
        for j, ms in enumerate(min_samples_values):
            dbscan, labels = hs.clustering.DBSCAN_RMSD_clustering(rmsd_matrix, eps, ms, parallel = n_processes)
            n_clusters[i,j] = len(np.unique(labels))
            if len(np.unique(labels)) > 2:
                ss[i,j] = metrics.silhouette_score(rmsd_matrix, labels)
                # dbcvs[i,j] = DBCV(rmsd_matrix, labels, dist_function=euclidean)

            else:
                ss[i,j] = np.nan
                # dbcvs[i,j] = np.nan

    m1 = (ss - np.nanmin(ss))/(np.nanmax(ss) - np.nanmin(ss))
    m2 = (-n_clusters - np.min(-n_clusters))/(np.max(-n_clusters) - np.min(-n_clusters))

    combination = 1.4*np.power(m1, 2) + np.power(m2, 2)

    # SS figure
    fig = plt.figure(dpi = 300, figsize = [5,5])
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

    # n_clusters figure
    fig = plt.figure(dpi = 300, figsize = [5,5])
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

    # combination metric figure
    fig = plt.figure(dpi = 300, figsize = [5,5])
    ax = fig.add_subplot(111)
    cax = ax.matshow(combination)
    fig.colorbar(cax)
    xaxis = np.arange(len(min_samples_values))
    yaxis = np.arange(len(eps_values))
    ax.set_xticks(xaxis)
    ax.set_yticks(yaxis)
    ax.set_xticklabels(min_samples_values)
    ax.set_yticklabels([round(eps, 2) for eps in eps_values])
    ax.set_xlabel("Min. Samples")
    ax.set_ylabel('$\epsilon_{DBSCAN}$')
    plt.savefig("grid_search_combo.png")

    # DBCV metric figure
    # fig = plt.figure(dpi = 300)
    # ax = fig.add_subplot(111)
    # cax = ax.matshow(dbcvs)
    # fig.colorbar(cax)
    # xaxis = np.arange(len(min_samples_values))
    # yaxis = np.arange(len(eps_values))
    # ax.set_xticks(xaxis)
    # ax.set_yticks(yaxis)
    # ax.set_xticklabels(min_samples_values)
    # ax.set_yticklabels([round(eps, 2) for eps in eps_values])
    # ax.set_xlabel("Min. Samples")
    # ax.set_ylabel('$\epsilon_{DBSCAN}$')
    # plt.savefig("grid_search_dbcv.png")

    # Get max value where silhouette score is maximized and cluster number is minimized

    max_metric = np.nanmax(combination)

    max_indices = list(zip(*np.where(combination == max_metric)))

    print("The maximum metric value of:", combination[max_indices[0]])
    print("N Clusters:", int(n_clusters[max_indices[0]]), "Silhouette Score", ss[max_indices[0]])
    print("Eps:", eps_values[max_indices[0][0]], "Min_samples:", min_samples_values[max_indices[0][1]])

    dbscan, labels = hs.clustering.DBSCAN_RMSD_clustering(rmsd_matrix, eps_values[max_indices[0][0]], min_samples_values[max_indices[0][1]])
    hs.clustering.write_clusters_to_file(labels, oct_traj)

    t2 = time.time()

    print("DBSCAN grid search took:", round(t2 - t1, 2), "seconds.")

    


if __name__ == "__main__":
    main()