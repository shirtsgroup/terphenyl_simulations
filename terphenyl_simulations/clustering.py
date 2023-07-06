import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import shutil as sh
from sklearn.cluster import DBSCAN
from sklearn import metrics, preprocessing
from .utils import backoff_directory
from .plotting import plot_grid_search


def construct_rmsd_matrix(traj):
    print("Constructing RMSD Matrix...")
    rmsd_matrix = np.zeros((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        rmsd = md.rmsd(
            traj, traj, frame=i, precentered=False, parallel=False
        )  # RMSD between frame i and j
        rmsd_matrix[i, :] = rmsd

    return rmsd_matrix


def DBSCAN_RMSD_clustering(rmsd_matrix, eps, min_samples, parallel=True):

    # Create DBSCAN object
    # print("Running DBSCAN clustering...")
    dbscan_parallel = None if parallel is True else parallel
    dbscan = DBSCAN(
        eps=eps, min_samples=min_samples, metric="precomputed", n_jobs=dbscan_parallel
    )
    dbscan.fit(rmsd_matrix)

    labels = dbscan.labels_
    cluster_ids = np.unique(labels)
    # print("Identified", len(cluster_ids), "cluster(s)!")
    return dbscan, labels

def write_medoids_to_file(
    labels,
    sil_scores,
    traj_object,
    output_dir = "clustering_output",
    prefix = None,
    output_format="gro",
    backoff = False
):
    if os.path.isdir(output_dir):
        if backoff:
            backoff_directory(output_dir)
            os.mkdir(output_dir)
        else:
            pass
    else:
        os.mkdir(output_dir)

    if prefix is not None:
        prefix += "_"
    else:
        prefix = ""

    for label in np.unique(labels):
        cluster_traj = traj_object[np.where(labels == label)]
        ss_cluster = sil_scores[np.where(labels == label)]
        medoid = cluster_traj[np.argmax(ss_cluster)]
        medoid_filename = os.path.join(
            output_dir, prefix + "medoid_" + str(label) + "." + output_format
        )
        medoid.save(medoid_filename)
    

def write_clusters_to_file(
    labels,
    traj_object,
    output_dir="clustering_output",
    prefix=None,
    output_format="gro",
    backoff = True
):

    # Write to output directory
    # print("Writing output files...")
    if os.path.isdir(output_dir):
        if backoff:
            backoff_directory(output_dir)
            os.mkdir(output_dir)
        else:
            pass
    else:
        os.mkdir(output_dir)

    if prefix is not None:
        prefix += "_"
    else:
        prefix = ""

    cluster_ids = np.unique(labels)
    for cluster_id in cluster_ids:
        frame_index = np.where(labels == cluster_id)
        cluster_traj = traj_object[frame_index]
        cluster_traj.superpose(cluster_traj, frame=0)
        cluster_filename = os.path.join(
            output_dir, prefix + "cluster_" + str(cluster_id) + "." + output_format
        )
        cluster_traj.save(cluster_filename)


def silhouette_score_metric(rmsd_matrix, dbscan_list_of_lists):
    result = np.zeros((len(dbscan_list_of_lists), len(dbscan_list_of_lists[0])))
    print("Calculating Silhouette Scores...")
    for i in tqdm(range(len(dbscan_list_of_lists))):
        for j in range(len(dbscan_list_of_lists[i])):
            labels = dbscan_list_of_lists[i][j]
            if len(np.unique(labels)) > 2:
                result[i, j] = metrics.silhouette_score(rmsd_matrix, labels)
            else:
                result[i, j] = np.nan
    return result

def n_clusters_metric(rmsd_matrix, dbscan_list_of_lists):
    result = np.zeros((len(dbscan_list_of_lists), len(dbscan_list_of_lists[0])))
    print("Getting N Clusters...")
    for i in tqdm(range(len(dbscan_list_of_lists))):
        for j in range(len(dbscan_list_of_lists[i])):
            labels = dbscan_list_of_lists[i][j]
            result[i,j] = len(np.unique(labels))
    return result


def combination_metric(rmsd_matrix, dbscan_list_of_lists):
    ss = np.zeros((len(dbscan_list_of_lists), len(dbscan_list_of_lists[0])))
    n_clusters = np.zeros((len(dbscan_list_of_lists), len(dbscan_list_of_lists[0])))
    print("Calculating Combined Metrics...")
    for i in tqdm(range(len(dbscan_list_of_lists))):
        for j in range(len(dbscan_list_of_lists[i])):
            labels = dbscan_list_of_lists[i][j]
            n_clusters[i,j] = len(np.unique(labels))
            if len(np.unique(labels)) > 2:
                ss[i,j] = metrics.silhouette_score(rmsd_matrix, labels)
            else:
                ss[i,j] = np.nan

    m1 = (ss - np.nanmin(ss)) / (np.nanmax(ss) - np.nanmin(ss))
    m2 = (-n_clusters - np.min(-n_clusters)) / (
        np.max(-n_clusters) - np.min(-n_clusters)
    )

    result = np.power(m1, 2) + np.power(m2, 2)

    return result

def plot_RMSD_histogram(rmsd_matrix, prefix = "remd"):
    rmsd_values = rmsd_matrix.reshape(-1)
    plt.hist(rmsd_values[rmsd_values > 0.00001], bins = 100, density = True)
    plt.xlabel("RMSD (nm)")
    plt.ylabel("Density")
    plt.savefig(prefix + "_rmsd_hist.png")

def clustering_grid_search(
    file_list,
    top_file,
    selection,
    n_min_samples=40,
    n_eps=40,
    n_processes=16,
    eps_limits=[0.01, 0.5],
    min_sample_limits=[0.005, 0.25],
    prefix="grid_search",
    frame_start=0,
    frame_end=-1,
    frame_stride=1,
    plot_filename = "ss.png",
    output_dir = "clustering_output",
    overwrite = True
):

    if overwrite == False:
        if os.path.isdir(output_dir):
            print("The output directory", output_dir, "already exists. " + 
                  "Skipping clustering. If you want to overwrite the existsing " +
                  "clustering ouput, set `overwrite = True`."
            )
            return


    # Load trajectory
    if type(file_list) == list:
        traj = md.load(file_list[0], top=top_file)
        for i in range(1, len(file_list)):
            tmp_traj = md.load(file_list[i], top = top_file)
            tmp_traj = tmp_traj[frame_start:frame_end:frame_stride]
            traj = traj.join(tmp_traj)

    else:
        traj = md.load(file_list, top=top_file)
        traj = traj[frame_start:frame_end:frame_stride]

    # Remove solvent
    top = traj.topology
    selection = top.select(selection)
    traj_object = traj.atom_slice(selection)

    # Construct RMSD matrix once
    rmsd_matrix = construct_rmsd_matrix(traj_object)

    plot_RMSD_histogram(rmsd_matrix)

    # Grid search for DBSCAN hyperparameters
    max_rmsd = np.max(rmsd_matrix)
    eps_values = np.linspace(eps_limits[0] * max_rmsd, eps_limits[1] * max_rmsd, n_eps)
    total_frames = traj_object.n_frames
    min_sample_values = np.linspace(
        min_sample_limits[0] * total_frames,
        min_sample_limits[1] * total_frames,
        n_min_samples,
        dtype=int,
    )
    min_samples_values = [int(a) for a in min_sample_values]
    ss = np.zeros((len(eps_values), len(min_sample_values)))
    n_clusters = np.zeros((len(eps_values), len(min_sample_values)))
    print("Performing grid search...")
    for i, eps in enumerate(tqdm(eps_values)):
        for j, ms in enumerate(min_samples_values):
            dbscan, labels = DBSCAN_RMSD_clustering(
                rmsd_matrix, eps, ms, parallel=n_processes
            )
            n_clusters[i, j] = len(np.unique(labels))
            if len(np.unique(labels)) > 1:
                ss[i, j] = metrics.silhouette_score(rmsd_matrix, labels)
            else:
                ss[i, j] = np.nan

    plot_grid_search(
        ss,
        min_samples_values,
        [round(eps, 2) for eps in eps_values],
        "Min. Samples",
        "$\epsilon_{DBSCAN}$",
        prefix  + "_" + plot_filename,
        "Avg. Silhouette Score"
    )

    # Get max value of first metric
    max_ss = np.nanmax(ss)
    max_indices = list(zip(*np.where(ss == max_ss)))
    print(max_indices)

    print("The maximum avg. silhouette score:", max_ss)
    print(
        "N Clusters:",
        int(n_clusters[max_indices[0]]),
        "Silhouette Score",
        ss[max_indices[0]],
    )
    print(
        "Eps:",
        eps_values[max_indices[0][0]],
        "Min_samples:",
        min_samples_values[max_indices[0][1]],
    )

    dbscan, labels = DBSCAN_RMSD_clustering(
        rmsd_matrix,
        eps_values[max_indices[0][0]],
        min_samples_values[max_indices[0][1]],
    )

    # Identify cluster medoids
    sil_scores = metrics.silhouette_samples(rmsd_matrix, labels)
    write_clusters_to_file(labels, traj_object, output_dir = output_dir)
    write_medoids_to_file(labels, sil_scores, traj_object, output_dir = output_dir)


def main():
    if not os.path.exists("npt_new_preloaded.h5"):
        # Load in trajectory
        test_traj = md.load_xtc(
            "../simulations/terphenyl_mop/tetramer_remd/t_200_350/sim0/npt_new.whole.xtc",
            top="../simulations/terphenyl_mop/tetramer_remd/t_200_350/sim0/npt_new.gro",
        )

        for i in range(1, 6):
            print(
                "Loading",
                "../simulations/terphenyl_mop/tetramer_remd/t_200_350/sim"
                + str(i)
                + "/npt_new.whole.xtc",
            )
            test_traj = test_traj.join(
                md.load_xtc(
                    "../simulations/terphenyl_mop/tetramer_remd/t_200_350/sim"
                    + str(i)
                    + "/npt_new.whole.xtc",
                    top="../simulations/terphenyl_mop/tetramer_remd/t_200_350/sim0/npt_new.gro",
                )
            )

        # Remove solvent
        top = test_traj.topology
        selection = top.select("resname TET")
        traj_object = test_traj.atom_slice(selection)

        # Frame Stride
        traj_object = traj_object[::]
        traj_object.save_hdf5("npt_new_preloaded.h5")
    else:
        traj_object = md.load("npt_new_preloaded.h5")

    # Run clustering
    dbscan, labels = DBSCAN_RMSD_clustering(traj_object, 0.3, 6 * 50, parallel=16)


if __name__ == "__main__":
    main()
