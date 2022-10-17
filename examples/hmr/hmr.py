import parmed as pmd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description = "A quick script to run HMR on  .top \
                    files"
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

    parser.add_argument(
        "--hmr_ratio",
        type = float,
        help = "ratio of the the HMR hydrogen mass to the original mass.\
            An hmr_ratio of 3 will scale the hydrogen masses by 3."
    )

    return parser.parse_args()    

def main():
    args = parse_args()
    gmx_top = pmd.load_file(args.top)

    for atom in gmx_top.atoms:
        if atom.element_name == "H":
            new_mass = atom.mass * args.hmr_ratio
            d_mass = new_mass - atom.mass
            for b_atom in atom.bond_partners:
                b_atom.mass -= d_mass
            atom.mass = new_mass

    gmx_top.write(args.output)

if __name__ == "__main__":
    main()