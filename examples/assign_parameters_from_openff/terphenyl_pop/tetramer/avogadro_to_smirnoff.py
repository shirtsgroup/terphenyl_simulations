# quick script to add different atom names per entry in a PDB
import argparse
from rdkit import Chem
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description = "A quick script to rename atoms output \
                       from Avogadro/ChemDraw to list them as \
                       unique atoms"
                    )

    parser.add_argument(
        "-f", "--file",
        type = str,
        help = "file name of original file to change"
    )

    parser.add_argument(
        "-o", "--output",
        type = str,
        help = "output file name to write new pdb file"
    )

    return parser.parse_args()


def main():
    
    args = parse_args()

    rdmol = Chem.rdmolfiles.MolFromPDBFile(args.file, removeHs=False)
    
    for atom in rdmol.GetAtoms():
        print(f"{atom.GetSymbol()}{atom.GetIdx()}")
        ri = atom.GetPDBResidueInfo()
        print(dir(ri))
        new_name = '{0:<4}'.format(atom.GetSymbol()+str(atom.GetIdx()+1))
        ri.SetName(new_name)
        ri.SetResidueName("TET")
        print(ri.GetResidueName())
        ri.SetIsHeteroAtom(False)

    Chem.rdmolfiles.MolToPDBFile(rdmol, args.output)

    


if __name__ == "__main__":
    main()