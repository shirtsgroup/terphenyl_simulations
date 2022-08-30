import MDAnalysis as mda
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.tests.datafiles import PSF, DCD
import numpy as np
import sys


def main():
    # Load .gro file
    u = mda.Universe("terphenyl_mop_octamer.itp", "mop_octamer.gro")

    # atom selection
    selected_atoms = u.select_atoms("all")

    # Run bond-angle-torsion analysis on structure
    R = BAT(selected_atoms)
    R.run()

    # Remove header from BAT array
    bat = R.results.bat[0, :]
    z_matrix_header = bat[0:9]
    z_matrix_body = bat[9:]

    # Each bond, angle and torsion entry listed sequentially in z_matrix_body
    # Here I extract out all the internal DOFs after the definition of the
    # first 3 atoms
    n_entries = int(len(z_matrix_body) / 3)
    bond_indices = [0, n_entries]
    bond_angle_indices = [n_entries, 2*n_entries]
    torsion_indices = [2*n_entries, 3*n_entries]

    bonds = z_matrix_body[0:n_entries]
    bond_angles = z_matrix_body[n_entries:2*n_entries]
    torsions = z_matrix_body[2*n_entries:3*n_entries]

    assert(len(bonds) == len(bond_angles))
    assert(len(bonds) == len(torsions))
    print(len(R._torsions))
    print(len(torsions))

    count = 0
    for i, torsion in enumerate(R._torsions):
        torsion_id = " ".join(torsion.names)
        element_id = ''.join([i for i in torsion_id if not i.isdigit()]).replace(" ", "")
        
        print("Torsion ID:", torsion_id)
        print("Angle:", torsions[i] * 180 / np.pi)
        # if torsion is a peptide bond, flip to cis bond
        if "HNCO" == element_id or "OCNH" == element_id or \
           "HNCC" == element_id or "CCNH" == element_id or \
           "CNCO" == element_id or "OCNC" == element_id or \
           "CNCC" == element_id or "CCNC" == element_id:
            print("Rotating", torsion_id, "by 180 degrees...")
            bat[9 + torsion_indices[0] + i] += np.pi
            count += 1
    
    print("Modified", count, "torsions.")
    
    # Reconstruct structure with modified torsions 
    XYZ = R.Cartesian(bat)
    u.atoms.positions = XYZ

    selected_atoms = u.select_atoms("all")
    selected_atoms.write("mop_octamer_cis.gro")
            


    


if __name__ == "__main__":
    main()