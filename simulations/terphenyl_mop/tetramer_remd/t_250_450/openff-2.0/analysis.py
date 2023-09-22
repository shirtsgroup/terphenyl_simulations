import terphenyl_simulations as hs
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing
import time
import sys
import os
import csv
import glob
from tqdm import tqdm
import MDAnalysis as mda

# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)


def main():
    t1 = time.time()
    # Clustering workflow
    hs.clustering.clustering_grid_search(
        [
            "sim0/npt_new.whole.xtc",
            "sim1/npt_new.whole.xtc",
            "sim2/npt_new.whole.xtc",
            "sim3/npt_new.whole.xtc",
            "sim4/npt_new.whole.xtc",
            "sim5/npt_new.whole.xtc",
            "sim6/npt_new.whole.xtc",
            "sim7/npt_new.whole.xtc",
            "sim8/npt_new.whole.xtc",
            "sim9/npt_new.whole.xtc",
        ],
        "sim0/berendsen_npt.gro",
        "resname TET or resname CAP",
        n_min_samples=40,
        n_eps=40,
        n_processes=32,
        prefix="grid_search",
        eps_limits=[0.01, 0.5],
        min_sample_limits=[0.005, 0.5],
        plot_filename = "ss.png",
        frame_stride = 10,
        overwrite = False
    )

    # Read in cluster outputs and REMD trajs
    cluster_file_list = glob.glob("clustering_output/cluster*")
    print("Loading Cluster trajectory files...")
    cluster_trajs = [
        md.load(gro_file) for gro_file in tqdm(cluster_file_list)
    ]

    remd_file_list = ["sim" + str(i) + "/npt_new.whole.xtc" for i in range(64)]
    print("Loading REMD trajectory files...")
    remd_trajs = [
        md.load(xtc_file, top="sim0/berendsen_npt.gro")
        for xtc_file in tqdm(remd_file_list)
    ]

    hexamer_u = mda.Universe("sim0/npt_new.tpr", "sim0/npt_new.gro")

    # Torsion definitions for first residue
    hexamer_r1_torsions = {
        "a1": ["C4", "C5", "C6", "C7"],
        "a2": ["C6", "C7", "C13", "C18"],
        "p1": ["C18", "C17", "C19", "N1"],
        "p2": ["C17", "C19", "N1", "H14"],
        "p3": ["O1", "C1", "C2", "C3"],
    }

    # Torsion Analysis
    cluster_entropy = np.zeros(len(cluster_trajs))

    for torsion_type in hexamer_r1_torsions.keys():
        print("Working on", torsion_type, "torsion...")
        torsion_atom_names = hs.utils.get_torsion_ids(
            hexamer_u, "TET", hexamer_r1_torsions[torsion_type], template_residue_i = 0
        )

        hs.plotting.plot_torsions_distributions(
            remd_trajs,
            torsion_atom_names,
            torsion_type + " Torsion (Degrees)",
            torsion_type + "_remd",
            torsion_type + " Torsion Plot",
            figsize = [5,5],
            cbar_params = [250, 450, "Temperature (K)"]
        )
        for i, traj in enumerate(cluster_trajs):

            entropy = hs.observables.calculate_torsion_entropy(traj, torsion_atom_names)
            cluster_entropy[i] += entropy

            hs.plotting.plot_torsions_distributions(
                traj,
                torsion_atom_names,
                torsion_type.upper() + " Torsion (Degrees)",
                "torsion_plots/" + torsion_type + "_" + cluster_file_list[i].split("/")[-1].split(".")[0],
                torsion_type + " Torsion Plot " + "Entropy: " + str(entropy),
                figsize=[5,5],
            )

    # Write entropy to file

    csv_header = [name.split("/")[-1].split(".")[0] for name in cluster_file_list]

    print(cluster_entropy)

    with open("cluster_torsion_entropy.csv", "w", newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(csv_header)
        csv_writer.writerow(cluster_entropy)
        

    t2 = time.time()
    print("Analysis took:", round(t2 - t1, 2), "seconds.")


if __name__ == "__main__":
    main()
