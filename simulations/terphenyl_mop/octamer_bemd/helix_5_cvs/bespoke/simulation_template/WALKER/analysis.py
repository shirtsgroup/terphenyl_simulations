import heteropolymer_simulations as hs
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing 
import time
import sys
import os
from tqdm import tqdm
import yaml
# sns.set_style('whitegrid')
# sns.color_palette("bwr", as_cmap=True)


# Going to move this into heteropolyerm_simulations package
def get_torsion_ids(universe, resname, torsion_id, template_residue_i = 0):
    """
    Using an MDAnalysis universe file with proper residue definitions, this function
    will extract the torsion ids of all torsions propagated along the chain. Specifically
    torsions should try to be fully defined within a single residue.
    """

    # Get index in residue
    atoms_in_residue = [a.name for a in universe.residues[template_residue_i].atoms]
    residue_atom_index  = [atoms_in_residue.index(a) for a in torsion_id if a in atoms_in_residue ]    
    dihedral_ids = []
    for residue in universe.residues:
        if residue.resname == resname:
            torsion_atoms = [residue.atoms[i] for i in residue_atom_index]
            dihedral_ids.append([ta.name for ta in torsion_atoms])
    
    return(dihedral_ids)



def main():
    t1 = time.time()

    # Torsion analysis

    # Read in MD trajectory (both as an MDAnlysis and MDTraj object)
    print("Loading trajectories...")
    traj = md.load("npt_new.whole.xtc", top = "berendsen_npt.gro")
    octamer_u = mda.Universe("npt_new.tpr", "npt_new.whole.xtc")
    
    # Load helix cluster trajectory and get torsion IDs
    with open("old_hexamer_torsion_ids.yml", 'r') as yaml_file:
        old_hexamer_torsions = yaml.safe_load(yaml_file)
    helix_cluster = md.load("build_system/clustering_output/cluster_9.gro")



    octamer_r1_torsions = {
        "a1" : ["C4", "C5", "C6", "C7"],
        "a2" : ["C6", "C7", "C13", "C14"],
        "p1" : ["C18", "C17", "C19", "N1"],
        "p2" : ["C17", "C19", "N1", "H14"],
        "p3" : ["O1", "C1", "C2", "C3"],
    }

    octamer_torsions = {}
    for torsion_type in octamer_r1_torsions.keys():
        t_ids = hs.utils.get_torsion_ids(octamer_u, "OCT", octamer_r1_torsions[torsion_type], template_residue_i=0)
        octamer_torsions[torsion_type] = t_ids
        
        # Need to adjust some of the torsion distributions
        # based on which atoms were used to measure the torsions
        # In the newer scripts we try to define torsions fully within
        # the residues, so torsion definitions may be off by 180 degrees

        offsets = None
        if torsion_type == "a2":
            offsets = [-np.pi, 0]
        if torsion_type == "p2":
            offsets = [np.pi, 0]

        # Get torsions from simulation
        hs.plotting.plot_torsions_distributions(
            [traj, helix_cluster],
            [octamer_torsions[torsion_type], old_hexamer_torsions[torsion_type]],
            torsion_type + " Torsion (Radians)",
            "torsion_plots/" + torsion_type,
            torsion_type + " Torsion Plot",
            legend = ["Helix MD", "Helix Cluster"],
            mirror_sym = False,
            offsets = offsets
        )

        hs.plotting.plot_torsion_timeseries(
            traj,
            octamer_torsions[torsion_type],
            [torsion_type + "_res_" + str(i+1) + ".png" for i in range(len(octamer_torsions[torsion_type]))],
            titles = [torsion_type.upper() + " torsion for residue " + str(i+1) for i in range(len(octamer_torsions[torsion_type]))],
            degrees = True
        )

    # Hydogen Bond analysis
    hbond_dir = "hbonds"
    hs.utils.make_path(hbond_dir)
    hbonds = HydrogenBondAnalysis(
        octamer_u,
        donors_sel = None,
        hydrogens_sel = "name H14 H33 H52 H71 H90 H109 H128 H147",
        acceptors_sel = "element O",
        d_a_cutoff = 3.5,
        d_h_a_angle_cutoff = 150,
        update_selections = False
    )

    hbonds.run(verbose=True)

    # Number of hydrogen bonds
    plt.figure(dpi = 300)
    plt.plot(hbonds.times, hbonds.count_by_time())
    plt.title("Number of hydrogen bonds over time", weight="bold")
    plt.xlabel("Time (ps)")
    plt.ylabel(r"$N_{HB}$")
    # plt.figtext("Mean: " + str(np.mean(hbonds.count_by_time())))
    # plt.figtext("STD: " + str(np.std(hbonds.count_by_time())))
    plt.savefig(hbond_dir + "/n_hydrogen_bonds.png")

    print("Average number of H-bonds:", np.mean(hbonds.count_by_time()))
    print("STD of H-bonds:", np.std(hbonds.count_by_time()))

    plt.close()

    # Distribution of hydrogen bond distances

    distances = hbonds.results.hbonds[:, 4]
    angles = hbonds.results.hbonds[:, 5]
    
    plt.figure(dpi = 300)
    plt.hist(distances, density = True, bins = 40)
    plt.title("Distribution of hydrogen bond distances", weight="bold")
    plt.xlabel("Distance (A)")
    plt.ylabel("Density")
    plt.savefig(hbond_dir + "/h_bond_distances.png")
    plt.close()

    plt.figure(dpi = 300)
    plt.hist(angles, density = True, bins = 40)
    plt.title("Distribution of hydrogen bond angles", weight="bold")
    plt.xlabel("Distance (A)")
    plt.ylabel("Density")
    plt.savefig(hbond_dir + "/h_bond_angles.png")
    plt.close()

    hb_types = list(hbonds.count_by_type())
    hb_types.sort(key = lambda x: int(x[2]), reverse = True)
    counts = np.array([int(hb_entry[2]) for hb_entry in hb_types])
    percentage = counts / octamer_u.trajectory.n_frames
    names = [hb_entry[0].split(":")[1] + " " + hb_entry[1].split(":")[1] for hb_entry in hb_types]

    plt.figure(dpi = 300)
    plt.bar(range(len(percentage)), percentage, tick_label = names)
    plt.xticks(rotation=90)
    plt.xlabel("Hydrogen Bond Types")
    plt.ylabel("Percentage of Total Simulation Time")
    plt.tight_layout()
    plt.savefig(hbond_dir + "/h_bond_types.png",  bbox_inches = "tight")
    plt.close()

    # Clustering first 100 ns
    # This will take sometime

    hs.clustering.clustering_grid_search(
        "npt_new.whole.xtc",
        "berendsen_npt.gro", 
        "resname OCT or resname CAP",
        n_min_samples = 40,
        n_eps = 40,
        n_processes = 16,
        prefix = "grid_search",
        min_sample_limits = [0.01, 0.1],
        eps_limits = [0.01, 0.4],
        frame_start = 0,
        frame_end = 500, # first 100 ns
        frame_stride = 1
    )
    

if __name__ == "__main__":
    main()