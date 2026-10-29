import MDAnalysis as mda
import numpy as np
import sys
from terphenyl_simulations.utils import get_torsion_ids

class InternalCoordinateEditor:
    def __init__(self, structure_file, topology_file, initial_atom = None):
        """
        Initialize an InternalCoordinateEditor object. This object can identify and modify
        specific bond lengths, bond-angles and torsions in a molecular structure.

        Parameters
        ----------
        structure_file : string
            String of structure file to edit
        topology_file : string
            String of the topology file used to identify internal coordinates. Compatible with all
            MDAnalysis compatible topology files.
        initial_atom : string, optional
            Name of the initial atom for internal coordinate definition (used by some analyses)
        """
        self._structure_file = structure_file
        self._top_file = topology_file
        self.universe = mda.Universe(self._top_file, self._structure_file, context = 'default', to_guess = ['types', 'masses', 'angles', 'bonds', 'dihedrals'])
        self.initial_atom = initial_atom

    @staticmethod
    def _get_downstream_group(universe, start_atom, blocked_atom):
        """Get all atoms reachable from start_atom without passing through blocked_atom"""
        # Build adjacency list from bonds
        bonds = universe.bonds
        adjacency = {}
        for bond in bonds:
            a1, a2 = bond.atoms[0], bond.atoms[1]
            if a1 not in adjacency:
                adjacency[a1] = []
            if a2 not in adjacency:
                adjacency[a2] = []
            adjacency[a1].append(a2)
            adjacency[a2].append(a1)

        # BFS from start_atom, avoiding blocked_atom
        visited = set()
        queue = [start_atom]
        visited.add(start_atom)

        while queue:
            current = queue.pop(0)
            if current in adjacency:
                for neighbor in adjacency[current]:
                    if neighbor == blocked_atom:
                        continue
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

        visited_indices = [atom.index for atom in visited]
        return universe.atoms[visited_indices]  # Return as AtomGroup

    @staticmethod
    def _rotate_group_around_axis(group, point_on_axis, axis_vector, angle):
        """Rotate group of atoms around axis defined by point and vector using Rodrigues' formula"""
        # Convert to numpy arrays
        point = np.array(point_on_axis)
        axis = np.array(axis_vector)
        # Normalize axis
        axis_norm = np.linalg.norm(axis)
        if axis_norm < 1e-12:
            return  # Avoid division by zero
        axis = axis / axis_norm

        # Get positions
        positions = group.positions

        # Translate so point_on_axis is at origin
        translated = positions - point

        # Rodrigues' formula: v_rot = v*cosθ + (k×v)*sinθ + k*(k·v)*(1-cosθ)
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)

        # Cross product k × v
        cross = np.cross(axis, translated)

        # Dot product k · v
        dot = np.dot(translated, axis)  # shape (N,)

        # Apply rotation
        rotated = translated * cos_theta + cross * sin_theta + axis * dot[:, np.newaxis] * (1 - cos_theta)

        # Translate back
        group.positions = rotated + point

    def set_torsions_from_template(self, template_universe, residue_name,
                                 template_residue_index=0):
        """
        Set torsions in the current universe to match those in a template universe
        by copying torsion values residue-by-residue with periodic boundary conditions.

        Parameters
        ----------
        template_universe : MDAnalysis.Universe
            Universe containing the template structure (e.g., helix)
        residue_name : str
            Name of the residue to match torsions for
        template_residue_index : int, default 0
            Index of the residue in template_universe to use as the pattern template

        Notes
        -----
        This method:
        1. Identifies a torsion pattern from the specified residue in the template
        2. Finds all matching torsions in both template and current universe
        3. Extracts torsion values per residue from the template
        4. Applies these values to the current universe with periodic boundary conditions
        5. Uses direct coordinate rotation (not BAT/Z-matrix) for better compatibility
           with aromatic systems and complex chemistry
        """
        # Get residues of interest
        template_residues = [r for r in template_universe.residues if r.resname == residue_name]
        current_residues = [r for r in self.universe.residues if r.resname == residue_name]

        if not template_residues:
            print(f"ERROR: No residues named {residue_name} found in template")
            return False
        if not current_residues:
            print(f"ERROR: No residues named {residue_name} found in current structure")
            return False
        if template_residue_index >= len(template_residues):
            print(f"ERROR: template_residue_index {template_residue_index} out of range")
            return False

        ref_residue = template_residues[template_residue_index]
        ref_resid = ref_residue.resid

        # Collect dihedral patterns fully contained in the reference residue
        patterns = []
        for dih in template_universe.dihedrals:
            if all(a.resid == ref_resid for a in dih):
                pattern = tuple(a.name for a in dih)
                if pattern not in patterns:
                    patterns.append(pattern)

        if not patterns:
            print(f"ERROR: No dihedral patterns found in reference residue")
            return False

        print(f"Found {len(patterns)} dihedral patterns in reference residue")

        torsions_set = 0
        total_torsions = 0

        for pattern in patterns:
            # Get torsion IDs for all residues using the helper
            template_ids = get_torsion_ids(template_universe, residue_name, list(pattern),
                                          template_residue_i=template_residue_index)
            current_ids = get_torsion_ids(self.universe, residue_name, list(pattern),
                                          template_residue_i=template_residue_index)

            if not template_ids or not current_ids:
                continue

            n_template = len(template_ids)
            n_current = len(current_ids)

            for i in range(n_current):
                t_idx = i % n_template
                t_names = template_ids[t_idx]
                c_names = current_ids[i]

                t_res = template_residues[t_idx]
                c_res = current_residues[i]

                try:
                    t_atoms = [template_universe.select_atoms(f"resid {t_res.resid} and name {n}")[0] for n in t_names]
                    c_atoms = [self.universe.select_atoms(f"resid {c_res.resid} and name {n}")[0] for n in c_names]
                except IndexError:
                    continue

                total_torsions += 1
                t_angle = mda.lib.distances.calc_dihedrals(*[a.position for a in t_atoms])
                c_angle = mda.lib.distances.calc_dihedrals(*[a.position for a in c_atoms])

                delta = t_angle - c_angle
                delta = (delta + np.pi) % (2 * np.pi) - np.pi

                if abs(delta) < 1e-6:
                    continue

                # Rotate downstream group around B-C bond
                point_on_axis = c_atoms[1].position
                axis_vector = c_atoms[2].position - c_atoms[1].position
                downstream_group = self._get_downstream_group(self.universe, c_atoms[3], c_atoms[2])
                if len(downstream_group) == 0:
                    continue
                self._rotate_group_around_axis(downstream_group, point_on_axis, axis_vector, delta)
                torsions_set += 1

        print(f'Successfully set {torsions_set} out of {total_torsions} torsions')
        return torsions_set > 0

    def write_structure(self, structure_file):
        """
        Function to rebuild structure with updated internal coordinates

        Parameters
        ----------
        structure_file : string
            String of filename to output new structure
        """
        selected_atoms = self.universe.select_atoms("all")
        selected_atoms.write(structure_file)