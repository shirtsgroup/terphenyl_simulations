import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import copy
from heteropolymer_simulations.edit_conf import InternalCoordinateEditor
import numpy as np
import sys

def main():
    # Here we want to identify the torsions accross the peptide bond and change it by
    # X degrees.

    angle = np.pi

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
            if "N" not in t_id.split(" ")[0] and "N" not in t_id.split(" ")[3]:
                if p_carbon not in t_id.split(" ")[0] and p_carbon not in t_id.split(" ")[3]:
                    torsion_ids.append(t_id)
                    torsions.append(t)

    # Here we want to identify which move propagates the chain the most
    # Some torsions with N and C* will propagate a single atom, rather than
    # the whole structure. So we use an RMSD between old and new structures
    # to identify the torsion which propagates a larger change in coordinates

    old_u = ice.universe.copy()
    old_torsions = copy.deepcopy(ice.torsions)
    old_ic_list = copy.deepcopy(ice.ic_list)

    final_torsion_ids = []
    for p_carbon in peptide_bond_c:
        # Get torsions with specific carbon
        t_ids = [t_id for t_id in torsion_ids if p_carbon in t_id.split(" ")]
        rmsds = []
        for t_id in t_ids:
            # Trial torsions by rotating them and evaluating with RMSD
            ice.set_torsion(t_id, torsions[torsion_ids.index(t_id)] + np.pi)
            ice.update_internal_coordinates()
            prop_rmsd = RMSD(old_u, ice.universe)
            prop_rmsd.run()

            # Reset to original coordinates
            ice.universe = old_u
            ice.torsions = old_torsions
            ice.ic_list = old_ic_list

            rmsds.append(prop_rmsd.results.rmsd[0][2])
        max_index = np.argmax(rmsds)
        final_torsion_ids.append(t_ids[max_index])
    
    print("The final torsions we will change are:")
    for t_id in final_torsion_ids:
        print(t_id)
    for t_id in final_torsion_ids:
        ice.set_torsion(t_id, torsions[torsion_ids.index(t_id)] + angle)
    ice.update_internal_coordinates()
    ice.write_structure("mop_octamer_cis.gro")    


if __name__ == "__main__":
    main()