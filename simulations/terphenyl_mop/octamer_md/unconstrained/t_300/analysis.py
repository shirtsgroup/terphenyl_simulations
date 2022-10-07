import heteropolymer_simulations as hs
import mdtraj as md
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn import preprocessing 
import time
import sys
import os
from tqdm import tqdm
# sns.set_style('whitegrid')
sns.color_palette("bwr", as_cmap=True)


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

    # Read in MD trajectory files
    traj = md.load("npt_new.whole.xtc", top = "berendsen_npt.gro")

    # Get hexamer torsion ids
    hexamer_u = mda.Universe("terphenyl_mop_hexamer.itp", "build_system/em_adjusted.gro")

    hexamer_r1_torsions = {
        "a1" : ["C4", "C5", "C6", "C7"],
        "a2" : ["C6", "C7", "C13", "C14"],
        "p1" : ["C18", "C17", "C19", "N1"],
        "p2" : ["C17", "C19", "N1", "H14"],
        "p3" : ["O1", "C1", "C2", "C3"],
    }

    octamer_torsions = {}
    for torsion_type in hexamer_r1_torsions.keys():
        t_ids = get_torsion_ids(hexamer_u, "HEX", hexamer_r1_torsions[torsion_type])
        octamer_torsions[torsion_type] = t_ids
    
        # Get torsions from simulation
        hs.plotting.plot_torsions_distributions(
            traj,
            octamer_torsions[torsion_type],
            torsion_type,
            torsion_type,
            torsion_type,
            legend = None,
            mirror_sym = False,
        )

    # Hydogen Bond analysis

    

    



if __name__ == "__main__":
    main()