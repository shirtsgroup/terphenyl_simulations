import MDAnalysis as mda
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.rms import RMSD
import sys
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

    def get_IC_lists(self, frame=0):
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
        self.bond_angle_indices = [n_entries, 2 * n_entries]
        self.torsion_indices = [2 * n_entries, 3 * n_entries]

        self.bonds = z_matrix_body[0:n_entries]
        self.bond_angles = z_matrix_body[n_entries : 2 * n_entries]
        self.torsions = z_matrix_body[2 * n_entries : 3 * n_entries]
        self.torsion_ids = [list(torsion.names) for torsion in self.bat._torsions]

        assert len(self.bonds) == len(self.bond_angles)
        assert len(self.bonds) == len(self.torsions)

    def find_bonds(self, atom_list):
        """
        Return a list of bond ids with a given atoms

        Parameters
        ----------
        atom_list : list
            List of strings of atom names used to identify bond

        Returns
        -------
        result : list
            List of bond ids including atoms in atom_list
        bond_length : list
            List of bond lengths corresponding to bond IDs
        """

        result = []
        bond_lengths = []
        for i, torsion_id in enumerate(self.torsion_ids):
            bond_id = torsion_id[0:2]
            include = 0
            for atom_id in atom_list:
                if any(c.isdigit() for c in atom_id):
                    if any(atom_id in b_atom for b_atom in bond_id):
                        include += 1
                else:
                    if atom_id in bond_id:
                        include += 1
            if include == len(atom_list):
                result.append(bond_id)
                bond_lengths.append(self.bonds[i])

        return result, bond_lengths

    def find_torsions(self, atom_list, positions=None):
        """
        Return a list of torsions ids from with a specific atoms

        Parameters
        ----------
        atom_list : list
            List of atom names used to identify torsions including those atoms. Here you
            can use specific atom names or elements names.
        positions : list
            List of integers indicating which positions to search for provided atom names

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
            torsion_id_list = torsion_id
            # Filter torsion_id positions based on given position selection
            if positions is not None:
                torsion_id_list = [
                    i for i in torsion_id_list if torsion_id_list.index(i) in positions
                ]
            for atom_id in atom_list:
                if any(c.isdigit() for c in atom_id):
                    if any(atom_id == t_atom for t_atom in torsion_id_list):
                        include += 1
                else:
                    if atom_id in torsion_id:
                        include += 1
            if include == len(atom_list):
                if include == len(atom_list):
                    result.append(torsion_id)
                    torsions.append(self.torsions[i])
        return result, torsions

    def identify_chain_prop_torsion(self, torsion_id_list):
        """
        This function identifies the torsion id that propagates downstream changes to
        the overall structure. Provided torsion IDs must contain the same internal atoms.

        Parameters
        ----------
        torsion_id_list : list
            List of lists of atom ids

        Returns
        -------
        torsion_id : list
            The torsion id that propagates movements downstream
        """
        # Check for the earliest definition of a given torsion
        # Earlier torsions propagate chains
        t_inds = [self.torsion_ids.index(t_id) for t_id in torsion_id_list]
        min_index = np.argmin(t_inds)
        return torsion_id_list[min_index]

    def get_torsion(self, torsion_id):
        """
        Get torsion value for a specified torsion id string

        Parameters
        ----------
        torsion_id : list
            Torsion ID for torsion value to return
        """
        if torsion_id in self.torsion_ids:
            torsion_index = self.torsion_ids.index(torsion_id)
            return self.torsions[torsion_index]
        else:
            print(
                " ".join(torsion_id),
                "is not a valid torsion ID. Available torsion \
                can be found in InternalCoordinateEditor.torsion_ids",
            )

    def get_bond_angle(self, bond_angle_id):
        """
        Get bond angle value for a specified angle id

        Parameters
        angle_id : string
            Angle ID for desired angle
        """

        ice_angle_ids = [t_id[1:] for t_id in self.torsion_ids]
        if bond_angle_id in ice_angle_ids:
            angle_index = ice_angle_ids.index(bond_angle_id)
            return self.bond_angles[angle_index]
        elif bond_angle_id == self.root_atoms or bond_angle_id == self.root_atoms[::-1]:
            return self.root_ics[2]
        else:
            print(bond_angle_id, "is not a valid angle ID.")

    def get_bond(self, bond_id):
        """
        Get bond length value for a specified bond id

        Parameters
        bond_id : list
            Bond ID for desired bond length
        """

        ice_bond_ids = [t_id[2:] for t_id in self.torsion_ids]
        if bond_id in ice_bond_ids:
            bond_index = ice_bond_ids.index(bond_id)
            return self.bond_angles[bond_index]
        elif bond_id in self.root_atoms[:2] or bond_id in self.root_atoms[:2][::-1]:
            return self.root_ics[0]
        elif bond_id in self.root_atoms[2:] or bond_id in self.root_atoms[2:][::-1]:
            return self.root_ics[1]
        else:
            print(bond_id, "is not a valid angle ID.")

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
            # print("Setting", torsion_id, "to", new_torsion * 180 / np.pi)
        else:
            print(torsion_id, "is not a valid torsion ID.")

    def set_angle(self, angle_id, new_angle):
        """
        Set a new angle within the InternalCoordinateEditor specified by its atom ids

        Parameters
        ----------
        angle_id : string
            Angle ID with atom names in angle
        new_angle : float
            New angle value in radians
        """

        # print("Angle ID:", angle_id)
        # print("fwd:", angle_id in [t_id[:-1] for t_id in self.torsion_ids])
        # print("bwd:", angle_id[::-1] in [t_id[:-1] for t_id in self.torsion_ids])
        angle_ids = [t_id[:-1] for t_id in self.torsion_ids]

        if angle_id in angle_ids:
            # print("1) Setting angle", angle_id, "to", new_angle * 180 / np.pi)
            angle_index = angle_ids.index(angle_id)
            self.bond_angles[angle_index] = new_angle
            self.ic_list[self.bond_angle_indices[0] + 9 + angle_index] = new_angle
        if angle_id[::-1] in angle_ids:
            # print("2) Setting angle", angle_id[::-1], "to", new_angle * 180 / np.pi)
            angle_index = angle_ids.index(angle_id[::-1])
            self.bond_angles[angle_index] = new_angle
            self.ic_list[self.bond_angle_indices[0] + 9 + angle_index] = new_angle
        if angle_id == list(self.root_atoms):
            # print("root atoms")
            self.root_ics[2] = new_angle
            self.ic_list[8] = new_angle
        if angle_id[::-1] == list(self.root_atoms):
            # print("root atoms")
            self.root_ics[2] = new_angle
            self.ic_list[8] = new_angle
        else:
            print(angle_id, "is not a valid angle ID.")

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
