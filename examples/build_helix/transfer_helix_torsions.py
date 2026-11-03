#!/usr/bin/env python3
"""
Standalone script to transfer helical torsions from a template structure
to a target oligomer.

Usage:
    python transfer_helix_torsions.py --template helix.gro --target target.gro \
        --residue-name MOP --output target_helix.gro
"""

import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals, calc_angles


def _get_downstream_group(universe, start_atom, blocked_atom):
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
    return universe.atoms[visited_indices]


def _rotate_group_around_axis(group, point_on_axis, axis_vector, angle):
    point = np.array(point_on_axis)
    axis = np.array(axis_vector)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-12:
        return
    axis = axis / axis_norm
    positions = group.positions
    translated = positions - point
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    cross = np.cross(axis, translated)
    dot = np.dot(translated, axis)
    rotated = translated * cos_theta + cross * sin_theta + axis * dot[:, np.newaxis] * (1 - cos_theta)
    group.positions = rotated + point


def _is_heavy(atom):
    if hasattr(atom, 'element') and atom.element is not None:
        return atom.element != 'H'
    return not atom.name.startswith('H')


def transfer_torsions(template_path, target_path, residue_name, output_path,
                      template_residue_index=0, window_size=2, verbose=False,
                      fit_angles=False, fit_impropers=False, heavy_only=False):
    template_u = mda.Universe(template_path, template_path, context='default', to_guess=['bonds','angles','dihedrals','impropers'])
    target_u = mda.Universe(target_path, target_path, context='default', to_guess=['bonds','angles','dihedrals','impropers'])

    template_residues = [r for r in template_u.residues if r.resname == residue_name]
    target_residues = [r for r in target_u.residues if r.resname == residue_name]

    if not template_residues:
        raise ValueError(f"No residues named {residue_name} found in template")
    if not target_residues:
        raise ValueError(f"No residues named {residue_name} found in target")
    if template_residue_index >= len(template_residues):
        raise ValueError(f"template_residue_index {template_residue_index} out of range")

    if template_residue_index + window_size > len(template_residues):
        raise ValueError(f"window_size {window_size} exceeds template residues from index {template_residue_index}")

    ref_residues = template_residues[template_residue_index:template_residue_index+window_size]
    # Build look up table for atom specific offset numbers
    atom_per_residue = {}
    for atom in ref_residues[0].atoms:
        atom_type = ''.join([c for c in atom.name if not c.isdigit()])
        atom_per_residue[atom_type] = atom_per_residue.get(atom_type, 0) + 1


    allowed_resids = {r.resid for r in ref_residues}
    base_resid = ref_residues[0].resid

    dihedral_patterns = []
    for dih in template_u.dihedrals:
        if all(a.resid in allowed_resids for a in dih):
            if heavy_only and not all(_is_heavy(a) for a in dih):
                continue
            pattern = tuple(a.name for a in dih)
            resid_offsets = tuple(a.resid - base_resid for a in dih)
            key = (pattern, resid_offsets)
            if key not in dihedral_patterns:
                dihedral_patterns.append(key)

    angle_patterns = []
    for ang in template_u.angles:
        if all(a.resid in allowed_resids for a in ang):
            if heavy_only and not all(_is_heavy(a) for a in ang):
                continue
            pattern = tuple(a.name for a in ang)
            resid_offsets = tuple(a.resid - base_resid for a in ang)
            key = (pattern, resid_offsets)
            if key not in angle_patterns:
                angle_patterns.append(key)

    improper_patterns = []
    try:
        for imp in template_u.impropers:
            if all(a.resid in allowed_resids for a in imp):
                if heavy_only and not all(_is_heavy(a) for a in imp):
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
        # ref_resid is the template's base_resid; default keeps offset 0.
        # Passing ref_resid=base_resid (template) when querying the target
        # enables atom_per_residue remapping for globally-incrementing names.
        if ref_resid is None:
            ref_resid = base_resid
        template_offset = base_resid - ref_resid
        for name, off in zip(pattern, resid_offsets):
            resid = base_resid + off
            atom_type = ''.join([c for c in name if not c.isdigit()])
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
        t_atoms = get_atoms(template_u, base_resid, pattern, resid_offsets)
        max_offset = max(resid_offsets)
        for i in range(len(target_residues) - max_offset):
            target_base_resid = target_residues[i].resid
            # Get atoms that match in the target structure
            try:
                c_atoms = get_atoms(target_u, target_base_resid, pattern, resid_offsets, ref_resid=base_resid)
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
            downstream = _get_downstream_group(target_u, c_atoms[0], c_atoms[1])
            if len(downstream) == 0:
                continue
            _rotate_group_around_axis(downstream, c_atoms[1].position, axis, delta)
            angles_set += 1

    # Impropers — fixed template lookup (mirrors dihedral logic)
    for pattern, resid_offsets in improper_patterns:
        # Get atoms from template region
        t_atoms = get_atoms(template_u, base_resid, pattern, resid_offsets)
        max_offset = max(resid_offsets)
        for i in range(len(target_residues) - max_offset):
            target_base_resid = target_residues[i].resid
            # Get atoms that match in the target structure
            try:
                c_atoms = get_atoms(target_u, target_base_resid, pattern, resid_offsets, ref_resid=base_resid)
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
            downstream = _get_downstream_group(target_u, c_atoms[0], central)
            if len(downstream) == 0:
                continue
            _rotate_group_around_axis(downstream, central.position, axis, delta)
            impropers_set += 1

    # Dihedrals
    for pattern, resid_offsets in dihedral_patterns:
        # Get atoms from template region
        t_atoms = get_atoms(template_u, base_resid, pattern, resid_offsets)
        max_offset = max(resid_offsets)
        for i in range(len(target_residues) - max_offset):
            target_base_resid = target_residues[i].resid
            # Get atoms that match in the target structure
            try:
                c_atoms = get_atoms(target_u, target_base_resid, pattern, resid_offsets, ref_resid=base_resid)
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
            downstream = _get_downstream_group(target_u, c_atoms[3], c_atoms[2])
            if len(downstream) == 0:
                continue
            _rotate_group_around_axis(downstream, point_on_axis, axis_vector, delta)
            torsions_set += 1

    print(f"Successfully set {torsions_set} out of {total_torsions} torsions, {angles_set} out of {total_angles} angles and {impropers_set} out of {total_impropers} impropers")
    target_u.atoms.write(output_path)
    print(f"Wrote output to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transfer helical torsions")
    parser.add_argument("--template", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--residue-name", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--template-residue-index", type=int, default=0)
    parser.add_argument("--window-size", type=int, default=2)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--fit-angles", action="store_true", help="Fit bond angles from template")
    parser.add_argument("--fit-impropers", action="store_true", help="Fit improper dihedrals from template")
    parser.add_argument("--heavy-only", action="store_true", help="Only transfer internal coordinates involving heavy atoms")
    args = parser.parse_args()
    transfer_torsions(args.template, args.target, args.residue_name, args.output,
                      args.template_residue_index, args.window_size, args.verbose,
                      args.fit_angles, args.fit_impropers, args.heavy_only)
