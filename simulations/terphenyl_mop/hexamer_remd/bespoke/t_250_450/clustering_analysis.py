import heteropolymer_simulations as hs
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing 
import time
import sys
import os
from tqdm import tqdm
# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)
    
        
def main():
    t1 = time.time()
    # Clustering workflow
    # hs.clustering.clustering_grid_search(["sim0/npt_new.whole.xtc", "sim1/npt_new.whole.xtc", "sim2/npt_new.whole.xtc", "sim3/npt_new.whole.xtc", "sim4/npt_new.whole.xtc"], "sim0/berendsen_npt.gro", "resname HEX", n_min_samples = 40, n_eps = 40, n_processes = 32, prefix = "grid_search", min_sample_limits = [0.005, 0.25])
    
    # Read in cluster outputs and REMD trajs
    cluster_file_list = os.listdir("clustering_output")
    print("Loading Cluster trajectory files...")
    cluster_trajs = [md.load("clustering_output/" + gro_file) for gro_file in tqdm(cluster_file_list)]

    # A1 Torsion
    print("Working on Aromatic Torsions 1...")
    torsion_atom_names = [
        ["C133", "C66", "C61", "C65"],
        ["C69", "C70", "C33", "C32"],
        ["C111", "C7", "C99", "C25"],
        ["C123", "C50", "C44", "C45"],
        ["C24", "C106", "C100", "C105"],
        ["C87", "C86", "C83", "C82"],
    ]
    for i, traj in enumerate(tqdm(cluster_trajs)):
        hs.plotting.plot_torsions_distributions(traj,
                                torsion_atom_names,
                                'Aromatic Torsion 1 (radians)',
                                "cluster_torsions/a1_" + cluster_file_list[i].split(".")[0],
                                cluster_file_list[i].split(".")[0],
                                )


    # A2 Torsion
    print("Working on Aromatic Torsions 2...")
    torsion_atom_names = [
        ["C61", "C65", "C74", "C73"],
        ["C33", "C32", "C72", "C9"],
        ["C99", "C57", "C53", "C1"],
        ["C44", "C11", "C47", "C48"],
        ["C100", "C105", "C108", "C109"],
        ["C83", "C84", "C92", "C93"],
    ]
    for i, traj in enumerate(tqdm(cluster_trajs)):
        hs.plotting.plot_torsions_distributions(traj,
                                torsion_atom_names,
                                "Aromatic Torsion 2 (radians)",
                                "cluster_torsions/a2_" + cluster_file_list[i].split(".")[0],
                                cluster_file_list[i].split(".")[0],
                                )

    # P1 Torsion
    print("Working on Peptide Torsions 1...")
    torsion_atom_names = [
        ["C73", "C130", "C136", "N6"],
        ["C9", "C114", "C116", "N1"],
        ["C1", "C56", "C31", "N3"],
        ["C48", "C35", "C124", "N4"],
        ["C21", "C15", "C36", "N2"],
        ["C95", "C94", "C98", "N5"],
    ]
    for i, traj in enumerate(tqdm(cluster_trajs)):
        hs.plotting.plot_torsions_distributions(traj,
                                torsion_atom_names,
                                "Peptide Torsion 1",
                                "cluster_torsions/p1_" + cluster_file_list[i].split(".")[0],
                                cluster_file_list[i].split(".")[0],
                                )

    # P2 Torsion
    print("Working on Peptide Torsions 2...")
    torsion_atom_names = [
        ["C130", "C136", "N6", "C28"],
        ["C114", "C116", "N1", "C43"],
        ["C56", "C31", "N3", "C2"],
        ["C35", "C124", "N4", "C37"],
        ["C15", "C36", "N2", "C97"],
        ["C94", "C98", "N5", "C16"],
    ]

    for i, traj in enumerate(tqdm(cluster_trajs)):
        hs.plotting.plot_torsions_distributions(traj,
                                torsion_atom_names,
                                "Peptide Torsion 2",
                                "cluster_torsions/p2_" + cluster_file_list[i].split(".")[0],
                                cluster_file_list[i].split(".")[0],
                                )

    # P3 Torsion
    print("Working on Peptide Torsions 3...")
    torsion_atom_names = [
        ["N6", "C28", "C68", "C20"],
        ["N1", "C43", "C119", "C117"],
        ["N3", "C2", "C8", "C34"],
        ["N4", "C37", "C40", "C39"],
        ["N2", "C97", "C89", "C88"],
        ["N3", "C2", "C8", "C34"],
    ]

    for i, traj in enumerate(tqdm(cluster_trajs)):
        hs.plotting.plot_torsions_distributions(traj,
                                torsion_atom_names,
                                "Peptide Torsion 3",
                                "cluster_torsions/p3_" + cluster_file_list[i].split(".")[0],
                                cluster_file_list[i].split(".")[0],
                                )

    p2_names = [
        [["C130", "C136", "N6", "C28"]],
        [["C114", "C116", "N1", "C43"]],
        [["C56", "C31", "N3", "C2"]],
        [["C35", "C124", "N4", "C37"]],
        [["C15", "C36", "N2", "C97"]],
        [["C94", "C98", "N5", "C16"]],
    ]

    # Plot individual torsions against neighboring peptide torsions
    for i, traj in enumerate(tqdm(cluster_trajs)):
        hs.plotting.plot_neighbor_dependent_2D_histogram(traj, p2_names, ["p2_r1", "p2_r2", "p2_r3", "p2_r4", "p2_r5", "p2_r6"], "cluster_torsions/p2_neighbors_" + cluster_file_list[i].split(".")[0], n_bins = 50)


    t2 = time.time()

    print("Analysis took:", round(t2 - t1, 2), "seconds.")

    

if __name__ == "__main__":
    main()
