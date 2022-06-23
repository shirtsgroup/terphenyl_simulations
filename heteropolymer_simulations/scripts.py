import argparse
import  heteropolymer_simulations as hs

def renumber_pdb_atoms():
    
    def parse_args():
        parser = argparse.ArgumentParser(
            description = "A quick script to rename atoms in a pdb \
                        file. This can be used to give unique atom \
                        to all atoms in the pdb file."
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

        parser.add_argument(
            "--resname",
            type = str,
            help = "output file name to write new pdb file",
            default = "TET"
        )


        return parser.parse_args()
    
    args = parse_args()
    hs.utils.renumber_pdb_atoms(args.file, args.output, args.resname)

def top_to_itp():
    def parse_args():
        parser = argparse.ArgumentParser(
            description = "A quick script to convert system .top \
                        files output from the openff-toolkit to  \
                        molecule specific .itp files."
                        )

        parser.add_argument(
            "-t", "--top",
            type = str,
            help = "file name of original .top file to convert"
        )

        parser.add_argument(
            "-o", "--output",
            type = str,
            help = "output file name to write new .itp file"
        )

        return parser.parse_args()

    args = parse_args()
    top = hs.utils.TopFileObject(args.top)
    hs.utils.write_itp_file(top, args.output)
    




