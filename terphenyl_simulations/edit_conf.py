import MDAnalysis as mda
import numpy as np
import sys
from MDAnalysis.lib.distances import calc_dihedrals, calc_angles
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
        self.universe = mda.Universe(self._top_file, self._structure_file, context = 'default', to_guess = ['types', 'masses', 'angles', 'bonds', 'dihedrals', 'impropers'])
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

    @staticmethod
    def _is_heavy(atom):
        if hasattr(atom, 'element') and atom.element is not None:
            return atom.element != 'H'
        return not atom.name.startswith('H')

    def set_torsions_from_template(self, template_universe, residue_name,
                                 template_residue_index=0, window_size=2,
                                 verbose=False, fit_angles=False,
                                 fit_impropers=False, heavy_only=False):
        """
        Set torsions in the current universe to match those in a template universe
        by copying internal coordinates with a residue window.

        Mirrors the logic in examples/build_helix/transfer_helix_torsions.py:
        collects dihedral/angle/improper patterns within a window of
        window_size residues, then applies them along the target chain
        with periodic repetition over template windows.

        Parameters
        ----------
        template_universe : MDAnalysis.Universe
            Universe containing the template structure (e.g., helix)
        residue_name : str
            Name of the residue to match torsions for
        template_residue_index : int, default 0
            Index of the starting residue in template_universe for the pattern window
        window_size : int, default 2
            Number of consecutive residues defining the pattern window
        verbose : bool, default False
            Print each internal coordinate change
        fit_angles : bool, default False
            If True, also transfer bond angles
        fit_impropers : bool, default False
            If True, also transfer improper dihedrals
        heavy_only : bool, default False
            If True, only transfer coordinates where all atoms are heavy atoms

        Notes
        -----
        This method:
        1. Identifies patterns from the windowed residues in the template
        2. Maps them via resid_offsets so cross-residue motifs are captured
        3. Applies values to the current universe with periodic boundary i % n_template_windows
        4. Uses direct coordinate rotation for angles/dihedrals/impropers
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

        if template_residue_index + window_size > len(template_residues):
            print(f"ERROR: window_size {window_size} exceeds template residues from index {template_residue_index}")
            return False

        ref_residues = template_residues[template_residue_index:template_residue_index+window_size]
        # Lookup for atom-name offset handling (mirrors standalone script)
        atom_per_residue = {}
        for atom in ref_residues[0].atoms:
            atom_type = ''.join([c for c in atom.name if not c.isdigit()])
            atom_per_residue[atom_type] = atom_per_residue.get(atom_type, 0) + 1

        allowed_resids = {r.resid for r in ref_residues}
        base_resid = ref_residues[0].resid

        dihedral_patterns = []
        for dih in template_universe.dihedrals:
            if all(a.resid in allowed_resids for a in dih):
                if heavy_only and not all(self._is_heavy(a) for a in dih):
                    continue
                pattern = tuple(a.name for a in dih)
                resid_offsets = tuple(a.resid - base_resid for a in dih)
                key = (pattern, resid_offsets)
                if key not in dihedral_patterns:
                    dihedral_patterns.append(key)

        angle_patterns = []
        for ang in template_universe.angles:
            if all(a.resid in allowed_resids for a in ang):
                if heavy_only and not all(self._is_heavy(a) for a in ang):
                    continue
                pattern = tuple(a.name for a in ang)
                resid_offsets = tuple(a.resid - base_resid for a in ang)
                key = (pattern, resid_offsets)
                if key not in angle_patterns:
                    angle_patterns.append(key)

        improper_patterns = []
        try:
            for imp in template_universe.impropers:
                if all(a.resid in allowed_resids for a in imp):
                    if heavy_only and not all(self._is_heavy(a) for a in imp):
                        continue
                    pattern = tuple(a.name for a in imp)
                    resid_offsets = tuple(a.resid - base_resid for a in imp)
                    key = (pattern, resid_offsets)
                    if key not in improper_patterns:
                        improper_patterns.append(key)
        except Exception:
            improper_patterns = []

        print(f"Found {len(dihedral_patterns)} dihedral patterns, {len(angle_patterns)} angle patterns and {len(improper_patterns)} improper patterns in residue window")

        if not fit_angles:
            angle_patterns = []
        if not fit_impropers:
            improper_patterns = []

        torsions_set = 0
        total_torsions = 0
        angles_set = 0
        total_angles = 0
        impropers_set = 0
        total_impropers = 0

        def get_atoms(universe, base_resid, pattern, resid_offsets, ref_resid=None):
            atoms = []
            # ref_resid defaults to base_resid so template_offset is 0 for
            # per-residue atom naming; if atom names are globally incremented
            # the offset arithmetic will remap them.
            if ref_resid is None:
                ref_resid = base_resid
            template_offset = base_resid - ref_resid
            for name, off in zip(pattern, resid_offsets):
                resid = base_resid + off
                atom_type = ''.join([c for c in name if not c.isdigit()])
                # Only apply numeric remapping if the base name contains digits
                # and we have a count for that element type
                if atom_type in atom_per_residue and any(c.isdigit() for c in name):
                    try:
                        atom_num = int(name.replace(atom_type, '')) + template_offset * atom_per_residue[atom_type]
                        atom_name = atom_type + str(atom_num)
                    except ValueError:
                        atom_name = name
                else:
                    atom_name = name
                sel = universe.select_atoms(f"resid {resid} and name {atom_name}")
                if len(sel) == 0:
                    # Fallback to original pattern name for per-residue naming
                    if atom_name != name:
                        sel = universe.select_atoms(f"resid {resid} and name {name}")
                    if len(sel) == 0:
                        raise IndexError(f"Atom {atom_name} resid {resid} not found")
                atoms.append(sel[0])
            return atoms

        # Angles — fixed template lookup so the same angle is applied to
        # every target window; atom_per_residue remaps globally-incrementing names
        for pattern, resid_offsets in angle_patterns:
            # Get atoms from template region
            t_atoms = get_atoms(template_universe, base_resid, pattern, resid_offsets)
            max_offset = max(resid_offsets)
            for i in range(len(current_residues) - max_offset):
                target_base_resid = current_residues[i].resid
                # Get atoms that match in the target structure
                try:
                    c_atoms = get_atoms(self.universe, target_base_resid, pattern, resid_offsets, ref_resid=base_resid)
                except IndexError:
                    continue
                total_angles += 1
                t_angle = calc_angles(*[a.position for a in t_atoms])
                c_angle = calc_angles(*[a.position for a in c_atoms])
                if np.isnan(t_angle) or np.isnan(c_angle):
                    continue
                delta = (t_angle - c_angle + np.pi) % (2 * np.pi) - np.pi
                if not np.isfinite(delta):
                    continue
                if abs(delta) < 1e-6:
                    continue
                if verbose:
                    print(f"Angle {pattern} target base resid {target_base_resid} offsets {resid_offsets} delta={delta:.4f} rad")
                v1 = c_atoms[0].position - c_atoms[1].position
                v2 = c_atoms[2].position - c_atoms[1].position
                axis = np.cross(v1, v2)
                axis_norm = np.linalg.norm(axis)
                if axis_norm < 1e-12:
                    continue
                axis = axis / axis_norm
                downstream = self._get_downstream_group(self.universe, c_atoms[0], c_atoms[1])
                if len(downstream) == 0:
                    continue
                self._rotate_group_around_axis(downstream, c_atoms[1].position, axis, delta)
                angles_set += 1

        # Impropers — fixed template lookup (mirrors dihedral logic)
        for pattern, resid_offsets in improper_patterns:
            # Get atoms from template region
            t_atoms = get_atoms(template_universe, base_resid, pattern, resid_offsets)
            max_offset = max(resid_offsets)
            for i in range(len(current_residues) - max_offset):
                target_base_resid = current_residues[i].resid
                # Get atoms that match in the target structure
                try:
                    c_atoms = get_atoms(self.universe, target_base_resid, pattern, resid_offsets, ref_resid=base_resid)
                except IndexError:
                    continue
                total_impropers += 1
                t_angle = calc_dihedrals(*[a.position for a in t_atoms])
                c_angle = calc_dihedrals(*[a.position for a in c_atoms])
                if np.isnan(t_angle) or np.isnan(c_angle):
                    continue
                delta = (t_angle - c_angle + np.pi) % (2 * np.pi) - np.pi
                if not np.isfinite(delta):
                    continue
                if abs(delta) < 1e-6:
                    continue
                if verbose:
                    print(f"Improper {pattern} target base resid {target_base_resid} offsets {resid_offsets} delta={delta:.4f} rad")
                central = c_atoms[1]
                v1 = c_atoms[2].position - central.position
                v2 = c_atoms[3].position - central.position
                axis = np.cross(v1, v2)
                axis_norm = np.linalg.norm(axis)
                if axis_norm < 1e-12:
                    continue
                axis = axis / axis_norm
                downstream = self._get_downstream_group(self.universe, c_atoms[0], central)
                if len(downstream) == 0:
                    continue
                self._rotate_group_around_axis(downstream, central.position, axis, delta)
                impropers_set += 1

        # Dihedrals — fixed template lookup; same torsion applied along chain
        for pattern, resid_offsets in dihedral_patterns:
            # Get atoms from template region
            t_atoms = get_atoms(template_universe, base_resid, pattern, resid_offsets)
            max_offset = max(resid_offsets)
            for i in range(len(current_residues) - max_offset):
                target_base_resid = current_residues[i].resid
                # Get atoms that match in the target structure
                try:
                    c_atoms = get_atoms(self.universe, target_base_resid, pattern, resid_offsets, ref_resid=base_resid)
                except IndexError:
                    continue
                total_torsions += 1
                t_angle = calc_dihedrals(*[a.position for a in t_atoms])
                c_angle = calc_dihedrals(*[a.position for a in c_atoms])
                if np.isnan(t_angle) or np.isnan(c_angle):
                    continue
                delta = (t_angle - c_angle + np.pi) % (2 * np.pi) - np.pi
                if not np.isfinite(delta):
                    continue
                if abs(delta) < 1e-6:
                    continue
                if verbose:
                    print(f"Dihedral {pattern} target base resid {target_base_resid} offsets {resid_offsets} delta={delta:.4f} rad")
                point_on_axis = c_atoms[1].position
                axis_vector = c_atoms[2].position - c_atoms[1].position
                downstream = self._get_downstream_group(self.universe, c_atoms[3], c_atoms[2])
                if len(downstream) == 0:
                    continue
                self._rotate_group_around_axis(downstream, point_on_axis, axis_vector, delta)
                torsions_set += 1

        print(f"Successfully set {torsions_set} out of {total_torsions} torsions, {angles_set} out of {total_angles} angles and {impropers_set} out of {total_impropers} impropers")
        return (torsions_set + angles_set + impropers_set) > 0

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
