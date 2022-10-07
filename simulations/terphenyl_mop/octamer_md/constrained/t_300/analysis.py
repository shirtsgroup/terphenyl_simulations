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
    hexamer_u = mda.Universe("npt_new.tpr", "npt_new.whole.xtc")
    
    # Load helix cluster trajectory and get torsion IDs
    with open("old_hexamer_torsion_ids.yml", 'r') as yaml_file:
        old_hexamer_torsions = yaml.safe_load(yaml_file)
    helix_cluster = md.load("build_system/clustering_output/cluster_9.gro")



    hexamer_r1_torsions = {
        "a1" : ["C4", "C5", "C6", "C7"],
        "a2" : ["C6", "C7", "C13", "C14"],
        "p1" : ["C18", "C17", "C19", "N1"],
        "p2" : ["C17", "C19", "N1", "H14"],
        "p3" : ["O1", "C1", "C2", "C3"],
    }

    octamer_torsions = {}
    for torsion_type in hexamer_r1_torsions.keys():
        t_ids = get_torsion_ids(hexamer_u, "OCT", hexamer_r1_torsions[torsion_type])
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
            torsion_type,
            torsion_type + " Torsion Plot",
            legend = ["Helix MD", "Helix Cluster"],
            mirror_sym = False,
            offsets = offsets
        )

    # Hydogen Bond analysis

    hbonds = HydrogenBondAnalysis(
        hexamer_u,
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
    plt.savefig("n_hydrogen_bonds.png")
    plt.close()

    # Distribution of hydrogen bond distances

    distances = hbonds.results.hbonds[:, 4]
    angles = hbonds.results.hbonds[:, 5]
    
    plt.figure(dpi = 300)
    plt.hist(distances, density = True, bins = 40)
    plt.title("Distribution of hydrogen bond distances", weight="bold")
    plt.xlabel("Distance (A)")
    plt.ylabel("Density")
    plt.savefig("h_bond_distances.png")
    plt.close()

    plt.figure(dpi = 300)
    plt.hist(angles, density = True, bins = 40)
    plt.title("Distribution of hydrogen bond angles", weight="bold")
    plt.xlabel("Distance (A)")
    plt.ylabel("Density")
    plt.savefig("h_bond_angles.png")
    plt.close()

    print(hbonds.count_by_type())

    



if __name__ == "__main__":
    main()