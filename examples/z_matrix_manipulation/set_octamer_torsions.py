import MDAnalysis as mda
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis.rms import RMSD
import copy
import numpy as np
import sys

class InternalCoordinateEditor:
    def __init__(self, structure_file, topology_file):
        """
        Initialize an InternalCoordinateEditor object. This object can extract
        all internal coordinates of a system and has tools to identify and modify
        specific bond lengths, bond-angles and torsions.

        Parameters
        ----------
        structure_file : string
            String of structure file to run MDAnalysis BAT method on
        topology_file : string
            String of the topology file used to identify internal coordinates. Compatible with all 
            MDAnalysis compatible topology files.
        """
        self._structure_file = structure_file
        self._top_file = topology_file
        self.universe = mda.Universe(self._top_file, self._structure_file)
        self.get_IC_lists()



    def get_IC_lists(self, frame = 0):
        """
        Function used to extract all torsions from a given structure file.
        This function assumes a single structure and will default to reading
        the first frame of a given trajectory.

        Parameters
        ----------
        frame : int
            Optional parameter to specify which frame of a trajectory to extract torsions.
        """
        # Load .gro file

        # atom selection
        selected_atoms = self.universe.select_atoms("all")

        # Run bond-angle-torsion analysis on structure
        self.bat = BAT(selected_atoms)
        self.bat.run()

        # Remove header from BAT array
        self.ic_list = self.bat.results.bat[frame, :]
        z_matrix_header = self.ic_list[0:9]
        z_matrix_body = self.ic_list[9:]
        
        # Extract root atoms and information
        self.root_atoms = self.bat._root.names
        self.root_atom_com = z_matrix_header[0:3]
        self.root_ext_angles = z_matrix_header[3:6]
        self.root_ics = z_matrix_header[6:9]

        # Extract remaining atoms and internal coorinates
        n_entries = int(len(z_matrix_body) / 3)
        self.bond_indices = [0, n_entries]
        self.bond_angle_indices = [n_entries, 2*n_entries]
        self.torsion_indices = [2*n_entries, 3*n_entries]

        self.bonds = z_matrix_body[0:n_entries]
        self.bond_angles = z_matrix_body[n_entries:2*n_entries]
        self.torsions = z_matrix_body[2*n_entries:3*n_entries]
        self.torsion_ids = [" ".join(torsion.names) for torsion in self.bat._torsions]

        assert(len(self.bonds) == len(self.bond_angles))
        assert(len(self.bonds) == len(self.torsions))

    def find_torsions(self, atom_list, set_op = "or"):
        """
        Return a list of torsions ids from with a specific atoms
        
        Parameters
        ----------
        atom_list : list
            List of atom names used to identify torsions including those atoms. Here you
            can use specific atom names or elements names.
        set_op : string
            Set operator defining how selection should be dealt with multiple atom identifiers.
        
        Returns
        -------
        result : list
            List of torsions ids including atoms in atom_list
        torsions : list
            List of torsions angles reported in radians
        """

        result = []
        torsions = []

        for i, torsion_id in enumerate(self.torsion_ids):
            include = 0
            for atom_id in atom_list:
                if any(i.isdigit() for i in atom_id):
                    if atom_id in torsion_id.split(" "):
                        include += 1
                else:
                    if atom_id in torsion_id:
                        include += 1
            if set_op == "or":
                if include > 0:
                    result.append(torsion_id)
                    torsions.append(self.torsions[i])
            if set_op == "and":
                if include == len(atom_list):
                    result.append(torsion_id)
                    torsions.append(self.torsions[i])
        
        return result, torsions
    
    def set_torsion(self, torsion_id, new_torsion):
        """
        Set a torsion, specified by its torsion id, to a new value

        Parameters
        ----------
        torsion_id : string
            Torsion ID for torsion to change
        new_torsion : float
            New torsion value in radians
        """

        if torsion_id in self.torsion_ids:
            torsion_index = self.torsion_ids.index(torsion_id)
            self.torsions[torsion_index] = new_torsion
            self.ic_list[self.torsion_indices[0] + 9 + torsion_index] = new_torsion

        else:
            print(torsion_id, "is not a valid torsion ID.")

    def update_internal_coordinates(self):
        """
        Update internal universe coordinates from BAT object
        """
        XYZ = self.bat.Cartesian(self.ic_list)
        self.universe.atoms.positions = XYZ

    def write_structure(self, structure_file):
        """
        Function to rebuild structure with updated internal coordiante

        Parameters
        ----------
        structure_file : string
            String of filename to output new structure
        """
        selected_atoms = self.universe.select_atoms("all")
        selected_atoms.write(structure_file)

def main():
    # Here we want to identify the torsions accross the peptide bond and change it by
    # X degrees.

    angle = np.pi

    # These are the carbon atoms present in these torsions
    peptide_bond_c = ["C5", "C58", "C39", "C78", "C16", "C146", "C127", "C166"]

    # Using the InternalCoordinateEditor we can interface with the MDAnalysis bond-angle-torsion analysis object
    # to change torsions of our structure
    ice = InternalCoordinateEditor("mop_octamer.gro", "terphenyl_mop_octamer.itp")
    
    # We use ICE to identify torsions with N elemets and each of the atoms in the
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