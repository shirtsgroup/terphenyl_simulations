import heteropolymer_simulations as hs
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing 
import time
import sys
from tqdm import tqdm
# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)
    
        
def main():
    t1 = time.time()
    # Clustering workflow
    hs.clustering.clustering_grid_search(["sim0/npt_new.whole.xtc"], "sim0/berendsen_npt.gro", "resname TET", n_min_samples = 40, n_eps = 40, n_processes =16, prefix = "grid_search")
    
    # Read in cluster outputs and REMD trajs
    cluster_file_list = ["clustering_output/cluster_0.gro", "clustering_output/cluster_1.gro", "clustering_output/cluster_-1.gro" ]
    print("Loading Cluster trajectory files...")
    cluster_trajs = [md.load(gro_file) for gro_file in tqdm(cluster_file_list)] 

    remd_file_list = ["sim" + str(i) + "/npt_new.xtc" for i in range(64)]
    print("Loading REMD trajectory files...")
    remd_trajs = [md.load(xtc_file, top = "sim0/berendsen_npt.gro") for xtc_file in tqdm(remd_file_list)]

    # A1 Torsion
    print("Working on Aromatic Torsions 1...")
    torsion_atom_names = [
        ["C7", "C86", "C80", "C85"],
        ["C72", "C67", "C64", "C65"],
        ["C37", "C26", "C20", "C25"],
        ["C48", "C47", "C44", "C45"],
    ]

    hs.plotting.plot_torsions_distributions(cluster_trajs,
                               torsion_atom_names,
                               'Aromatic Torsion 1 (radians)',
                               "a1",
                               legend = ["cluster_0.gro", "cluster_1.gro", "cluster_-1.gro" ],
                               mirror_sym = True
                            )
    
    hs.plotting.plot_torsions_distributions(remd_trajs,
                               torsion_atom_names,
                               "Aromatic Torsion 1 (radians)",
                               "remd_a1",
                               mirror_sym = True
                            )


    # A2 Torsion
    print("Working on Aromatic Torsions 2...")
    torsion_atom_names = [
        ["C80", "C85", "C88", "C89"],
        ["C64", "C65", "C73", "C74"],
        ["C20", "C25", "C28", "C29"],
        ["C44", "C45", "C53", "C54"],
    ]

    hs.plotting.plot_torsions_distributions(cluster_trajs,
                               torsion_atom_names,
                               "Aromatic Torsion 2 (radians)",
                               "a2",
                               legend = ["cluster_0.gro", "cluster_1.gro", "cluster_-1.gro" ],
                            )

    hs.plotting.plot_torsions_distributions(remd_trajs,
                               torsion_atom_names,
                               "Aromatic Torsion 2",
                               "remd_a2",
                               mirror_sym = True
                            )


    # P1 Torsion
    print("Working on Peptide Torsions 1...")
    torsion_atom_names = [
        ["C89", "C4", "C15", "N1"],
        ["C74", "C75", "C79", "N4"],
        ["C29", "C30", "C38", "N2"],
    ]

    hs.plotting.plot_torsions_distributions(cluster_trajs,
                               torsion_atom_names,
                               "Peptide Torsion 1",
                               "p1",
                               legend = ["cluster_0.gro", "cluster_1.gro", "cluster_-1.gro" ]
                            )

    hs.plotting.plot_torsions_distributions(remd_trajs,
                               torsion_atom_names,
                               "Peptide Torsion 1",
                               "remd_p1",
                            )

    # P2 Torsion
    print("Working on Peptide Torsions 2...")
    torsion_atom_names = [
        ["C4", "C15", "N1", "C78"],
        ["C75", "C79", "N4", "C39"],
        ["C30", "C38", "N2", "C58"],
    ]

    hs.plotting.plot_torsions_distributions(cluster_trajs,
                               torsion_atom_names,
                               "Peptide Torsion 2",
                               "p2",
                               legend = ["cluster_0.gro", "cluster_1.gro", "cluster_-1.gro" ]
                            )

    hs.plotting.plot_torsions_distributions(remd_trajs,
                               torsion_atom_names,
                               "Peptide Torsion 2",
                               "remd_p2",
                            )

    # P3 Torsion
    print("Working on Peptide Torsions 3...")
    torsion_atom_names = [
        ["N1", "C78", "C70", "C71"],
        ["N4", "C39", "C35", "C36"],
        ["N2", "C58", "C50", "C49"],
    ]

    hs.plotting.plot_torsions_distributions(cluster_trajs,
                               torsion_atom_names,
                               "Peptide Torsion 3",
                               "p3",
                               legend = ["cluster_0.gro", "cluster_1.gro", "cluster_-1.gro"],
                               mirror_sym = True
                            )

    hs.plotting.plot_torsions_distributions(remd_trajs,
                               torsion_atom_names,
                               "Peptide Torsion 3",
                               "remd_p3",
                               mirror_sym = True
                            )

    p2_names = [
        [["C4", "C15", "N1", "C78"]],
        [["C75", "C79", "N4", "C39"]],
        [["C30", "C38", "N2", "C58"]],
        [["C55", "C59", "N3", "C5"]]
    ]

    hs.plotting.plot_neighbor_dependent_2D_histogram(cluster_trajs[0], p2_names, ["p2_r1", "p2_r2", "p2_r3", "p2_r4"], "p2_neighbors_cluster_0", n_bins = 50)
    hs.plotting.plot_neighbor_dependent_2D_histogram(remd_trajs[0], p2_names, ["p2_r1", "p2_r2", "p2_r3", "p2_r4"], "p2_neighbors_remd_0", n_bins = 50)


    t2 = time.time()

    print("Analysis took:", round(t2 - t1, 2), "seconds.")

    

if __name__ == "__main__":
    main()