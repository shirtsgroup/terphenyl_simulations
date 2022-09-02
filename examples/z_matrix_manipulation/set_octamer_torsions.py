import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import copy
from heteropolymer_simulations.edit_conf import InternalCoordinateEditor
import numpy as np
import sys

def main():
    # Here we want to identify the torsions accross the peptide bond and change it by
    # X degrees.

    angle = -np.pi/8

    # These are the carbon atoms present in these torsions
    peptide_bond_c = ["C5", "C58", "C39", "C78", "C16", "C146", "C127", "C166"]

    # Using the InternalCoordinateEditor we can interface with the MDAnalysis bond-angle-torsion analysis object
    # to change torsions of our structure
    ice = InternalCoordinateEditor("mop_octamer.gro", "terphenyl_mop_octamer.itp")
    
    # We use ICE to identify torsions with N elements and each of the atoms in the
    # center two positions
    torsion_ids = []
    torsions = []
    for p_carbon in peptide_bond_c:
        nc_t_ids, ts = ice.find_torsions(["N", p_carbon], set_op = "and")
        for t_id, t in zip(nc_t_ids, ts):
            # Here we filter torsions that have N and carbonyl carbons no in the middle of the torsion
            if "N" not in t_id.split(" ")[0] and "N" not in t_id.split(" ")[3]:
                if p_carbon not in t_id.split(" ")[0] and p_carbon not in t_id.split(" ")[3]:
                    torsion_ids.append(t_id)
                    torsions.append(t)

    # Here we want to identify which move propagates the chain the most
    # Some torsions with N and C* will propagate a single atom, rather than
    # the whole structure. So we use an RMSD between old and new structures
    # to identify the torsion which propagates a larger change in coordinates

    old_coordinates = copy.deepcopy(ice.universe.atoms.positions)
    old_torsions = copy.deepcopy(ice.torsions)
    old_ic_list = copy.deepcopy(ice.ic_list)

    final_torsion_ids = []
    for p_carbon in peptide_bond_c:
        # Get torsions with specific carbon
        t_ids = [t_id for t_id in torsion_ids if p_carbon in t_id.split(" ")]
        
        # If there are multiple torsions defined across the bond we specifed
        # We verify which torsion propagates all atoms on the torsion
        if len(t_ids) > 1:
            # Get list of all atoms in all identified 
            t_ids_list = [t_id.split(" ") for t_id in t_ids]

            # Flatten list of lists
            flat_list = np.array([item for sublist in t_ids_list for item in sublist])

            # Get instances of unique atoms
            atoms_unique, counts = np.unique(flat_list, return_counts=True)
            rotated_atoms = atoms_unique[counts < len(t_ids)]

            d_atoms = []

            for t_id in t_ids:
                # Trial torsions by rotating them and evaluating with RMSD
                ice.set_torsion(t_id, torsions[torsion_ids.index(t_id)] + angle)
                ice.update_internal_coordinates()

                # Calculate distance moved by both rotated atoms
                # t_id that gives the largest change is the correct torsion
                d_tot = 0
                for atom in rotated_atoms:
                    atom_index = ice.universe.select_atoms("name " + atom)[0].index
                    d_tot += np.linalg.norm(ice.universe.atoms.positions[atom_index] - old_coordinates[atom_index])                                
                d_atoms.append(d_tot)

                ice.write_structure("test_"+t_id.replace(" ", "_") + ".gro")

                # Reset to original coordinates
                # ice.universe.atoms.positions = copy.deepcopy(old_coordinates)
                ice.torsions = copy.deepcopy(old_torsions)
                ice.ic_list = copy.deepcopy(old_ic_list)
                ice.update_internal_coordinates()

            max_index = np.argmax(d_atoms)
            final_torsion_ids.append(t_ids[max_index])
        else:
            final_torsion_ids.append(t_ids[0])
    
    print("The final torsions we will change are:")
    for t_id in final_torsion_ids:
        print(t_id)
        ice.write_structure("test_1.gro")    
        ice.set_torsion(t_id, torsions[torsion_ids.index(t_id)] + angle)
        ice.update_internal_coordinates()
        ice.write_structure("test_2.gro")

    ice.write_structure("mop_octamer_cis.gro")   


if __name__ == "__main__":
    main()