#!/usr/bin/env python
import getpass
import os
from datetime import datetime
import rdkit
from rdkit import Chem
import platform

class TopFileObject:
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename, 'r') as f:
            self.top_file = f.readlines()
        self.parse_file()
    
    def parse_file(self):
        section_name = None
        section_data = {}
        for i, line in enumerate(self.top_file):
            # Skip commented lines
            if  line[0] == ";":
                continue
            # New section
            if "[" in line  and "]" in line:
                section_name = line.split()[1]
                section_data[section_name] = {"data" : []}
                # Use legends to store values a dict
                if self.top_file[i+1][0] == ";":
                    labels = self.top_file[i+1]
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

def write_itp_file(top_object, filename, itp_sections = None):

    if itp_sections is None:
        itp_sections = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals"]

    if ".itp" not in filename:
        filename += ".itp"
    with open(filename, "w") as f:
        # write header
        f.write("; itp topology file\n")
        f.write("; Generated using heteropolymer_simulations.util.write_itp_file()\n")
        f.write("; Original file: " + top_object.filename + "\n")
        f.write("; Author: " + getpass.getuser() + "\n")
        f.write("; Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
        f.write("; Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
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
        new_name = '{0:<4}'.format(atom.GetSymbol()+str(atom.GetIdx()+1))
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
    prefix_wo_file = "/".join(prefix.split("/")[:-1])
    if not os.path.exists(prefix_wo_file):
        os.makedirs(prefix_wo_file)


def main():
    top = TopFileObject("test.top")
    write_itp_file(top, "test")

if __name__ == "__main__":
    main()



                




                

            