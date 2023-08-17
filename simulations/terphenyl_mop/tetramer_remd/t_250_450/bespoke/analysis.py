import terphenyl_simulations as ts
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
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)


def main():
    t1 = time.time()

    commmon_factor = (450/250)**(1/63)
    temperatures = [250 * (commmon_factor)**i for i in range(64)]

    # Clustering workflow
    ts.clustering.clustering_grid_search(
        [
            "sim0/npt_new.whole.xtc",
            "sim1/npt_new.whole.xtc",
            "sim2/npt_new.whole.xtc",
            "sim3/npt_new.whole.xtc",
            "sim4/npt_new.whole.xtc",
            "sim5/npt_new.whole.xtc",
            "sim6/npt_new.whole.xtc",
        ],
        "sim0/berendsen_npt.gro",
        "resname TET or resname CAP",
        n_min_samples=40,
        n_eps=40,
        n_processes=16,
        prefix="grid_search",
        min_sample_limits=[0.01, 0.15],
        eps_limits=[0.05, 0.3],
        frame_stride = 1,
        overwrite = False
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
        "a1": ["C5", "C6", "C7", "C8"],
        "a2": ["C7", "C8", "C15", "C20"],
        "p1": ["C20", "C19", "C21", "N22"],
        "p2": ["C19", "C21", "N22", "H39"],
        "p3": ["O2", "C1", "C3", "C4"],
    }

    print(remd_trajs)

    # Torsion Analysis
    ts.utils.make_path("torsion_plots")
    for torsion_type in hexamer_r1_torsions.keys():
        if torsion_type == -1:
            continue
        print("Working on", torsion_type, "torsion...")
        torsion_atom_names = ts.utils.get_torsion_ids(
            hexamer_u, "TET", hexamer_r1_torsions[torsion_type], template_residue_i = 0
        )

        ts.plotting.plot_torsions_distributions(
            remd_trajs,
            torsion_atom_names,
            torsion_type + " Torsion (radians)",
            "torsion_plots/" + torsion_type + "_remd",
            torsion_type + " Torsion Plot",
            cbar_limits = [250, 450]
        )
        
        for i, traj in enumerate(cluster_trajs):
            ts.plotting.plot_torsions_distributions(
                traj,
                torsion_atom_names,
                torsion_type.upper() + " Torsion (radians)",
                "torsion_plots/" + torsion_type + "_" + cluster_file_list[i].split(".")[0],
                torsion_type + " Torsion Plot",
                figsize = [5,5]
            )

        

    t2 = time.time()
    print("Analysis took:", round(t2 - t1, 2), "seconds.")


if __name__ == "__main__":
    main()
