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
        "resname HEX or resname CAP",
        n_min_samples=40,
        n_eps=40,
        n_processes=32,
        prefix="grid_search",
        eps_limits=[0.01, 0.2],
        min_sample_limits=[0.005, 0.1],
        plot_filename = "ss.png",
        frame_stride = 2
    )

    # Read in cluster outputs and REMD trajs
    cluster_file_list = os.listdir("clustering_output")
    print("Loading Cluster trajectory files...")
    cluster_trajs = [
        md.load("clustering_output/" + gro_file) for gro_file in tqdm(cluster_file_list)
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
        "a1": ["C51", "C50", "C49", "C58"],
        "a2": ["C58", "C59", "C60", "C61"],
        "p1": ["C62", "C63", "C64", "N66"],
        "p2": ["C63", "C64", "N66", "H85"],
        "p3": ["O55", "C54", "C53", "C52"],
    }

    # Torsion Analysis
    hs.utils.make_path("torsion_plots")
    for torsion_type in hexamer_r1_torsions.keys():
        print("Working on", torsion_type, "torsion...")
        torsion_atom_names = hs.utils.get_torsion_ids(
            hexamer_u, "HEX", hexamer_r1_torsions[torsion_type], template_residue_i = 1
        )

        hs.plotting.plot_torsions_distributions(
            remd_trajs,
            torsion_atom_names,
            torsion_type + "Torsion (radians)",
            torsion_type + "_remd",
            torsion_type + " Torsion Plot",
            figsize = [5,5],
            cbar_limits = [250, 450]
        )
        
        for i, traj in enumerate(cluster_trajs):
            hs.plotting.plot_torsions_distributions(
                traj,
                torsion_atom_names,
                torsion_type.upper() + " Torsion (radians)",
                "torsion_plots/" + torsion_type + "_" + cluster_file_list[i].split(".")[0],
                torsion_type + " Torsion Plot",
                figsize=[5,5]
            )

    t2 = time.time()
    print("Analysis took:", round(t2 - t1, 2), "seconds.")


if __name__ == "__main__":
    main()
