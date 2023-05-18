#!/usr/bin/env python
import getpass
import os
from datetime import datetime
import rdkit
from rdkit import Chem
import platform
import numpy as np
import re
import shutil as sh
import sys


class TopFileObject:
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename, "r") as f:
            self.top_file = f.readlines()
        self.parse_file()

    def parse_file(self):
        section_name = None
        section_data = {}
        for i, line in enumerate(self.top_file):
            # Skip commented lines
            if line[0] == ";":
                continue
            # New section
            if "[" in line and "]" in line:
                section_name = line.split()[1]
                section_data[section_name] = {"data": []}
                # Use legends to store values a dict
                if self.top_file[i + 1][0] == ";":
                    labels = self.top_file[i + 1]
                    section_data[section_name]["labels"] = labels
                else:
                    labels = None
                continue
            if section_name is not None:
                entries = line.split()
                if len(entries) > 0:
                    section_data[section_name]["data"].append(line)
                else:
                    continue
        self.__dict__.update(**section_data)


def write_itp_file(top_object, filename, itp_sections=None):

    if itp_sections is None:
        itp_sections = [
            "moleculetype",
            "atoms",
            "bonds",
            "pairs",
            "angles",
            "dihedrals",
        ]

    if ".itp" not in filename:
        filename += ".itp"
    with open(filename, "w") as f:
        # write header
        f.write("; itp topology file\n")
        f.write("; Generated using heteropolymer_simulations.util.write_itp_file()\n")
        f.write("; Original file: " + top_object.filename + "\n")
        f.write("; Author: " + getpass.getuser() + "\n")
        f.write("; Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
        f.write("; Time: " + datetime.now().strftime("%I:%M%p") + "\n")
        f.write("; System: " + platform.platform() + "\n\n")

        for itp_s in itp_sections:
            if itp_s in dir(top_object):
                f.write("[ " + itp_s + " ]\n")
                if getattr(top_object, itp_s)["labels"] is not None:
                    f.write(getattr(top_object, itp_s)["labels"])
                for line in getattr(top_object, itp_s)["data"]:
                    f.write(line)
                f.write("\n\n")


def renumber_pdb_atoms(pdb_file, out_pdb):
    rdmol = Chem.rdmolfiles.MolFromPDBFile(pdb_file, removeHs=False)

    for atom in rdmol.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        new_name = "{0:<4}".format(atom.GetSymbol() + str(atom.GetIdx() + 1))
        ri.SetName(new_name)
        ri.SetIsHeteroAtom(False)

    Chem.rdmolfiles.MolToPDBFile(rdmol, out_pdb)


def make_path(prefix):
    """
    Utility function to make a directory for a given prefix if it doesnt
    exist yet.

    Parameters
    ----------
    prefix : string
        String of a directory path
    """

    prefix_list = prefix.split("/")
    if "." in prefix_list[-1]:
        prefix_list.pop()
    prefix_wo_file = "/".join(prefix_list)
    if len(prefix_wo_file) == 0:
        return
    if not os.path.exists(prefix_wo_file):
        os.makedirs(prefix_wo_file)
    else:
        backoff_directory(prefix)
        os.makedirs(prefix_wo_file)


def backoff_directory(dir_name):
    """
    Utility function to backoff a directory if it exists already

    Parameters
    ----------
    dir_name : string
        Directory name
    """
    i_backoff = 1
    old_path = "#" + dir_name + "." + str(i_backoff) + "#"
    while os.path.isdir(old_path):
        i_backoff += 1
        old_path = "#" + dir_name + "." + str(i_backoff) + "#"
    print("Backoff! Moving the old " + dir_name + " to " + old_path)
    sh.move(dir_name, old_path)

def get_torsion_ids(universe, resname, torsion_id, template_residue_i = 1):
    """
    Using an MDAnalysis universe file with proper residue definitions, this function
    will extract the torsion ids of all torsions propagated along the chain. Specifically
    torsions should try to be fully defined within a single residue.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        MDAnalysis universe of trajectory files used to get torsion ids
    resname : string
        Residue name from which the torsions are defined
    torsion_id : list of strings
        Atom names that make up the torsion in the provided residue
    template_residue_i : int, default 0
        Residue in which the atom names are defined

    Returns
    -------
    dihedral_ids : list of lists of strings
        List of atom names defining the torsion provided by torsion_id
        propaged in for each residue with name resname.
    """

    # Get index in residue
    atoms_in_residue = [a.name for a in universe.select_atoms("resid " + str(template_residue_i + 1)).atoms]

    # Case where all atoms in dihedral are defined the given residue
    if all([atom_id in atoms_in_residue for atom_id in torsion_id]):
        residue_atom_index  = [atoms_in_residue.index(a) for a in torsion_id if a in atoms_in_residue ]    
        dihedral_ids = []
        for residue in universe.residues:
            if residue.resname == resname:
                torsion_atoms = [residue.atoms[i] for i in residue_atom_index]
                dihedral_ids.append([ta.name for ta in torsion_atoms])
    # Case where 1 or more atoms in dihedral are not in given residue
    elif any([atom_id in atoms_in_residue for atom_id in torsion_id]):
        if all([universe.select_atoms("name " + atom_id).resnames[0] == resname for atom_id in torsion_id]):
            dihedral_ids = []
            resid_in_torsion = np.array([universe.select_atoms("name " + atom_id).resids[0] for atom_id in torsion_id])
            resid_diff = resid_in_torsion - template_residue_i - 1
            atom_resid_res_inds = []
            for i, atom in enumerate(torsion_id):
                atoms_in_res_i = [a.name for a in universe.select_atoms("resid " + str(resid_in_torsion[i])).atoms]
                atom_res_index = atoms_in_res_i.index(atom)
                atom_resid_index = (resid_in_torsion[i] - template_residue_i - 1, atom_res_index)
                atom_resid_res_inds.append(atom_resid_index)
            for i in range(len(universe.residues)):
                if np.max(np.array(resid_diff) + i) == len(universe.residues):
                    continue
                if all([universe.residues[i].resname == universe.residues[i + diff].resname for diff in resid_diff]) and universe.residues[i].resname == resname:
                    torsion_atoms = [universe.residues[i + diff].atoms[j].name for diff, j in atom_resid_res_inds]
                    dihedral_ids.append(torsion_atoms)
    else:
        print("All atoms in the dihedral are not present in the template residue. " + \
            "Please change `template_residue_i` to be consistent with the residue  " + \
            "where the atom names are originally from.")
        return(None)

        
    return(dihedral_ids)

def get_angle_ids(universe, resname, angle_id, template_residue_i = 0):
    """
    Using an MDAnalysis universe with proper residue definitions, this function
    will extract the angle atom ids of all angles propagated along the chain. Specifically
    bond-angles should try to be fully defined within a single residue.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        MDAnalysis universe of trajectory files used to get torsion ids
    resname : string
        Residue name from which the torsions are defined
    torsion_id : list of strings
        Atom names that make up the torsion in the provided residue
    template_residue_i : int, default 0
        Residue in which the atom names are defined

    Returns
    -------
    angle_ids : list of lists of strings
        List of atom names defining the torsion provided by torsion_id
        propaged in for each residue with name resname.
    """

    # Get index in residue
    atoms_in_residue = [a.name for a in universe.select_atoms("resid " + str(template_residue_i + 1)).atoms]
    
    if all([atom_id in atoms_in_residue for atom_id in angle_id]):
        residue_atom_index  = [atoms_in_residue.index(a) for a in angle_id if a in atoms_in_residue ]    
        angle_ids = []
        for residue in universe.residues:
            if residue.resname == resname:
                angle_atoms = [residue.atoms[i] for i in residue_atom_index]
                angle_ids.append([ta.name for ta in angle_atoms])
        return(angle_ids)

    # Case where 1 or more atoms in dihedral are not in given residue
    elif any([atom_id in atoms_in_residue for atom_id in angle_id]):
        if all([universe.select_atoms("name " + atom_id).resnames[0] == resname for atom_id in angle_id]):
            angle_ids = []
            resid_in_torsion = np.array([universe.select_atoms("name " + atom_id).resids[0] for atom_id in angle_id])
            resid_diff = resid_in_torsion - template_residue_i - 1
            atom_resid_res_inds = []
            for i, atom in enumerate(angle_id):
                atoms_in_res_i = [a.name for a in universe.select_atoms("resid " + str(resid_in_torsion[i])).atoms]
                atom_res_index = atoms_in_res_i.index(atom)
                atom_resid_index = (resid_in_torsion[i] - template_residue_i - 1, atom_res_index)
                atom_resid_res_inds.append(atom_resid_index)
            for i in range(len(universe.residues)):
                if np.max(np.array(resid_diff) + i) == len(universe.residues):
                    continue
                if all([universe.residues[i].resname == universe.residues[i + diff].resname for diff in resid_diff]) and universe.residues[i].resname == resname:
                    torsion_atoms = [universe.residues[i + diff].atoms[j].name for diff, j in atom_resid_res_inds]
                    angle_ids.append(torsion_atoms)
    else:
        pass
    return(angle_ids)


def replace_all_pattern(pattern, replace, file):
    """
    Function to replace patterns present in with a particular string

    Parameters
    ----------
    pattern : string
        string pattern to replace in file
    replace : string
        string that will replace the pattern found in file
    file : string
        path to file
    """
    with open(file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = re.sub(pattern, replace, lines[i])
    with open(file, "w") as f:
        f.writelines(lines)

def main():
    top = TopFileObject("test.top")
    write_itp_file(top, "test")


if __name__ == "__main__":
    main()
