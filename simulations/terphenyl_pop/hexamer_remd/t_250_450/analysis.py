import terphenyl_simulations as hs
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing
from subprocess import Popen, PIPE
import time
import sys
import os
from tqdm import tqdm
import shutil
from glob import glob
import MDAnalysis as mda

# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)


def main():
    t1 = time.time()

    # Check for gmx executable and if npt_new.whole.xtc files exist
    # And apply -pbc whole to all files
    if shutil.which("gmx") and len(glob("sim*/npt_new.whole.xtc")) != len(glob("sim*")):
        for xtc_file in glob("sim*/npt_new.xtc"):
            p = Popen(["gmx", "trjconv", "-f", xtc_file, "-s", "sim0/npt_new.tpr", "-pbc", "whole", "-o", xtc_file.split(".")[0] + ".whole.xtc"], stdin=PIPE, stdout=PIPE, bufsize=1)
            p.communicate(input = b'0\n')
    if shutil.which("gmx_mpi") and len(glob("sim*/npt_new.whole.xtc")) != len(glob("sim*")):
        for xtc_file in glob("sim*/npt_new.xtc"):
            p = Popen(["mpirun", "-np", "1", "gmx_mpi", "trjconv", "-f", xtc_file, "-s", "sim0/npt_new.tpr", "-pbc", "whole", "-o", xtc_file.split(".")[0] + ".whole.xtc"], stdin=PIPE, stdout=PIPE, bufsize=1)
            p.communicate(input = b'0\n')


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
        eps_limits=[0.01, 0.3],
        min_sample_limits=[0.005, 0.3],
        plot_filename = "ss.png",
        frame_stride = 3
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
        "a1": ["C49", "C50", "C51", "C52"],
        "a2": ["C51", "C52", "C59", "C60"],
        "p1": ["C61", "C62", "C65", "N66"],
        "p2": ["C62", "C65", "N66", "H83"],
        "p3": ["O46", "C45", "C47", "C48"],
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
            cbar_params = [250, 450, "Temperature (K)"]
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
