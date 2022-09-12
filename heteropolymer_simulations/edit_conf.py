import MDAnalysis as mda
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.rms import RMSD
import copy
import numpy as np


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
                if any(c.isdigit() for c in atom_id):
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

    def filter_torsions(self, torsion_id, torsions, exclusion_atoms, exlusion_positions):
        """
        This function is to filter torsion IDs based on 
        """

        # Should rework above function to include position  specification
        pass

    def identify_chain_prop_torsion(self, torsion_id_list, ):
        """
        This function identifies the torsion id that propagates downstream changes to
        the overall structure. Provided torsion IDs must contain the same internal atoms.

        Parameters
        ----------
        torsion_id_list : list
            List of strings of torsion ids

        Returns
        -------
        torsion_id : string
            The torsion id that propagates movements downstream
        """
        central_atoms = torsion_id_list[0].split(" ")[1:3]
        
        for c_atoms in central_atoms:
            pass

    
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