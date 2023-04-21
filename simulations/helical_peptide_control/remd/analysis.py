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
import MDAnalysis as mda

# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)


def main():
    t1 = time.time()
    # Clustering workflow
    hs.clustering.clustering_grid_search(
        [ "sim" + str(i) + "/npt_new.whole.xtc" for i in range(64)],
        "sim0/berendsen_npt.gro",
        "resname GLU or resname LYS or resname NTC or resname CTC",
        n_min_samples=40,
        n_eps=40,
        n_processes=32,
        prefix="grid_search",
        eps_limits=[0.1, 0.5],
        min_sample_limits=[0.005, 0.15],
        plot_filename = "ss.png",
        frame_stride = 15
    )

    t2 = time.time()
    print("Analysis took:", round(t2 - t1, 2), "seconds.")


if __name__ == "__main__":
    main()
