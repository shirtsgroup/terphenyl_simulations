import heteropolymer_simulations as hs
from tqdm import tqdm
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 10})

def main():

    # Load in REMD trajectory
    remd_file_list = ["sim" + str(i) + "/npt_new.whole.xtc" for i in range(64)]
    print("Loading REMD trajectory files...")
    remd_trajs = [
        mda.Universe("sim0/berendsen_npt.gro", xtc_file)
        for xtc_file in tqdm(remd_file_list)
    ]

    # Load in helix medoid structure
    helix_medoid = mda.Universe("clustering_output/medoid_0.gro")
    folded_percentage = []
    error = []

    for i, traj in enumerate(remd_trajs):
        print("Replica", str(i))
        traj_rmsd = mda.analysis.rms.RMSD(
                        traj,
                        helix_medoid, 
                        select = "resname TET or resname CAP",
                        ref_frame=0
                    )
        traj_rmsd.run()
        plt.figure()
        plt.hist(traj_rmsd.rmsd[500:,2], density = True, bins = 30)
        plt.xlabel("RMSD to Helix Medoid (A)")
        plt.ylabel("Density")
        plt.savefig("test_plots/replica_" + str(i) + "_rsmd_hist.png")
        plt.close()

        folded = [int(item) for item in list(traj_rmsd.rmsd[500:,2] < 3)]
        folded_percentage.append(np.mean(folded))
        error.append(1.96 * np.std(folded)/np.sqrt(len(folded)))
    plt.figure(figsize = [5,5])
    plt.errorbar(list(range(len(remd_trajs))), folded_percentage, yerr = error)
    plt.xlabel("Replica Index")
    plt.ylabel("Percentage Folded")
    plt.savefig("folded_replica_percentage.png", dpi = 300, transparent=True)
    plt.close()

    common_ratio = np.power(450/250, 1/(64-1))
    t_range = [250 * (common_ratio ** i) for i in range(len(remd_trajs))]

    plt.figure(figsize = [5,5])
    plt.errorbar(t_range, folded_percentage, yerr = error)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Fraction Folded")
    plt.savefig("folded_T_percentage.png", dpi = 300, transparent=True)
    plt.close()



    

    # Per residue helix definition
        



if __name__ == "__main__":
    main()