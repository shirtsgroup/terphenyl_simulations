import terphenyl_simulations as ts
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing
import time
import sys
import os
import glob
from tqdm import tqdm
import MDAnalysis as mda

# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)


def main():
    t1 = time.time()
    # Clustering workflow
    ts.clustering.clustering_grid_search(
        [ "sim" + str(i) + "/npt_new.whole.xtc" for i in range(10)],
        "sim0/berendsen_npt.gro",
        "resname GLU or resname LYS or resname NTC or resname CTC",
        n_min_samples=40,
        n_eps=40,
        n_processes=32,
        prefix="grid_search",
        eps_limits=[0.1, 0.5],
        min_sample_limits=[0.005, 0.2],
        plot_filename = "ss.png",
        frame_stride = 2,
        overwrite = False
    )


    # Plot unbiased Native H-bonds for each simulation
    # Should turn this in to an an analysis workflow script
    cur_dir = os.path.abspath("")
    for sim_dir in tqdm(glob.glob("_sim*")):
        os.chdir(sim_dir)
        print(os.path.abspath(""))
        subprocess.Popen("plumed driver --plumed ../plumed_hbond_dist.dat --mf_xtc npt_new.xtc --igro npt_new.gro", shell=True).wait()
        data_file = "HBOND_SUMS"
        with open(data_file) as f:
            headers = f.readline().strip()
        headers = headers.split(" ")[2:]
        hb_data = pd.read_csv(data_file, sep = " ", header = None, comment="#", names = headers)
        plt.figure()
        plt.plot(hb_data.time/5, hb_data.n_hbonds)
        plt.xlabel("Time (ns)")
        plt.ylabel("N Native H-Bonds")
        half_index = int(len(hb_data.n_hbonds) / 2)
        print(half_index)
        plt.title("Mean: " + str(round(np.mean(hb_data.n_hbonds.values[half_index:]), 5)) +  " STD: " +  str(round(np.std(hb_data.n_hbonds.values[half_index:]), 5)))
        plt.savefig("remd_hbonds.png", dpi = 300)
        plt.close()
        os.chdir(cur_dir)
    
    ts.utils.make_path("ramachandran_plots")
    for sim_dir in tqdm(glob.glob("sim*")):        
        # Pull temperature from mdp file
        with open(os.path.join(sim_dir, "npt_new.mdp")) as f:
            title = None
            for line in f.readlines():
                if "gen-temp" in line:
                    title = "Temperature: " + str(round(float(line.split()[-1]), 1)) + " K"

        config_points = [os.path.join(sim_dir, fn) for fn in ["em_solvated.gro", "berendsen_nvt.gro", "berendsen_npt.gro"]]

        ts.plotting.plot_ramachandran_plot(
            os.path.join(sim_dir, "npt_new.whole.xtc"),
            os.path.join(sim_dir, "npt_new.gro"),
            prefix = "ramachandran_plots/"+sim_dir,
            bins = 100,
            title = title,
            scatter_points_files = config_points
        )

    t2 = time.time()
    print("Analysis took:", round(t2 - t1, 2), "seconds.")


if __name__ == "__main__":
    main()
