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

ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) # This is your Project Root

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
                else:
                    labels = None
                section_data[section_name]["labels"] = labels
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


class GromacsLogFile:
    def __init__(self, filename):
        with open(filename, "r") as fn:
            self.log_file_lines = fn.readlines()
        self._extract_info()

    def _find_expression(self, expression):
        lines = list(filter(lambda x: expression in x, self.log_file_lines))
        lines = [line.strip() for line in lines]
        return lines
    
    def _find_text_block(self, expression, lines_after):
        line_idx= [i for i in range(len(self.log_file_lines)) if expression in self.log_file_lines[i]]
        block = []
        for i in line_idx:
            lines = [line.strip() for line in self.log_file_lines[i:i+lines_after+1]]
            block.append(lines)
        return block
    
    def _extract_info(self):
        # Basic information about simulation
        print("Extracting simulation information...")
        self.gmx_version = self._find_expression("GROMACS version:")[0].split()[-1]
        self.cmd_line_call = self._find_text_block("Command line:", lines_after = 1)[0][-1].strip()
        self.hostname = self._find_expression("Hardware detected on host")[0].split()[4]
        self.working_dir = self._find_expression("Working dir:")[0].split()[-1]
        self.input_mdp = self._find_text_block("Input Parameters:", 174)[0]
        

        print("Extracting Observables Output...")
        # Time block
        time_blocks = self._find_text_block("Step           Time", 1)
        obs_blocks = self._find_text_block("Energies (kJ/mol)", 8)
        observables = {
            "step" : [],
            "time" : [],
            "bond" : [],
            "angle" : [],
            "proper_dih" : [],
            "per_imp_dih" : [],
            "lj_14" : [],
            "coulomb_14" : [],
            "LJ_sr" : [],
            "disper_corr" : [],
            "coulomb_sr" : [],
            "coulomb_recip" : [],
            "potential" : [],
            "kinetic" : [],
            "total" : [],
            "conserved" : [],
            "temperature" : [],
            "pres_dc" : [],
            "pressure" : [],
            "constr_rmsd" : []
        }

        for time_block in time_blocks:
            observables["step"].append(int(time_block[1].split()[0]))
            observables["time"].append(float(time_block[1].split()[1]))

        self.n_steps = max(observables["step"])

        # Exclude the last observable block because this it is the averages over
        # The entire simulation
        for obs_block in obs_blocks[:-1]:
            observables["bond"].append(float(obs_block[2].split()[0]))
            observables["angle"].append(float(obs_block[2].split()[1]))
            observables["proper_dih"].append(float(obs_block[2].split()[2]))
            observables["per_imp_dih"].append(float(obs_block[2].split()[3]))
            observables["lj_14"].append(float(obs_block[2].split()[4]))
            observables["coulomb_14"].append(float(obs_block[4].split()[0]))
            observables["LJ_sr"].append(float(obs_block[4].split()[1]))
            observables["disper_corr"].append(float(obs_block[4].split()[2]))
            observables["coulomb_sr"].append(float(obs_block[4].split()[3]))
            observables["coulomb_recip"].append(float(obs_block[4].split()[4]))
            observables["potential"].append(float(obs_block[6].split()[0]))
            observables["kinetic"].append(float(obs_block[6].split()[1]))
            observables["total"].append(float(obs_block[6].split()[2]))
            observables["conserved"].append(float(obs_block[6].split()[3]))
            observables["temperature"].append(float(obs_block[6].split()[4]))
            observables["pres_dc"].append(float(obs_block[8].split()[0]))
            observables["pressure"].append(float(obs_block[8].split()[1]))
            observables["constr_rmsd"].append(float(obs_block[8].split()[2]))

        self.observables = observables
        if self._find_expression("Order After Exchange:"):
            # Extract replica exchange state trajectory
            print("Extracting state trajectories...")
            sp_lines = self._find_expression("Order After Exchange:")
            sp_lines = [sl.replace("Order After Exchange:", "") for sl in sp_lines]
            states = [sl.replace("x","").split() for sl in sp_lines]
            states = [[int(s_i) for s_i in state] for state in states]
            self.states = states
            self.n_states = max(states[0]) + 1

            # Collect Empirical Exchange Transition Matrix
            ex_matrix = self._find_text_block("Empirical Transition Matrix", self.n_states+1)

            if len(ex_matrix) != 0:
                ex_probs = []
                for i in range(self.n_states):
                    ex_prob_str = ex_matrix[i+2].split()[1:self.n_states+1]
                    ex_probs.append([float(p) for p in ex_prob_str])
                
                self.transition_matrix = np.array(ex_probs)


def main():
    pass


if __name__ == "__main__":
    main()
