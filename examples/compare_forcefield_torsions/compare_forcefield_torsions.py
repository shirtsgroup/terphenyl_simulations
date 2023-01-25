import argparse
import parmed
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

def parse_args():
    parser = argparse.ArgumentParser(
        description = "A script to compare torsion parameters \
            between two gromacs topology files"
    )
    parser.add_argument(
        "--top_files", 
        type = str, 
        help = "list of topology files to compare torsions",
        nargs = "+"
    )
    parser.add_argument(
        "--legend_labels", 
        type = str, 
        help = "list of topology files to compare torsions",
        nargs = "+"
    )
    parser.add_argument(
        "--output_dir",
        type = str,
        help = "name of directory to ouput comparsion plots",
        default = "output"
    )
    
    return parser.parse_args()

def extract_dihedral_parameters(dihe_list):
    dihe_dict = {}
    for dihe in dihe_list:
        id = " ".join([dihe.atom1.name, dihe.atom2.name, dihe.atom3.name, dihe.atom4.name])
        if id in dihe_dict.keys():
            dihe_dict[id].append(dihe.type)
        else:
            dihe_dict[id] = [dihe.type]
    return dihe_dict

def plot_torsion_potentials(dihe_potential_list, torsion_id, file_name, legend_labels):
    x = np.linspace(-np.pi, np.pi, 100)
    plt.figure()
    for dihe_potential in dihe_potential_list:
        y = 0
        for dihe_type in dihe_potential[torsion_id]:
            y += dihe_type.phi_k * (1 + np.cos(dihe_type.per * x - dihe_type.phase * np.pi / 180))
        plt.plot(x, y)
    plt.title(torsion_id)
    plt.xlabel("Torsion Angle (radians)")
    plt.ylabel("Potential Energy")
    plt.legend(legend_labels)
    plt.savefig(file_name)
    plt.close()


def main():
    args = parse_args()

    # Make output directory
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    else:
        backoff_dir = args.output_dir
        backoff_count = 1
        while os.path.exists(backoff_dir):
            backoff_dir = "#" + args.output_dir + "." + str(backoff_count) + "#"
            backoff_count += 1
        shutil.move(args.output_dir, backoff_dir)
        os.mkdir(args.output_dir)

    top_files = []
    for fn in args.top_files:
        top_files.append(parmed.load_file(fn))
    
    dihedral_potentials = []
    for parameters in top_files:
        dihes = parameters.dihedrals
        dihe_dict = extract_dihedral_parameters(dihes)
        dihedral_potentials.append(dihe_dict)

    for i, torsion_id in enumerate(dihedral_potentials[0].keys()):
        if len(dihedral_potentials[0][torsion_id]) != len(dihedral_potentials[1][torsion_id]):
            plot_torsion_potentials(dihedral_potentials, torsion_id, args.output_dir + "/torsion_" + str(i) + ".png", args.legend_labels)
        

if __name__ == "__main__":
    main()