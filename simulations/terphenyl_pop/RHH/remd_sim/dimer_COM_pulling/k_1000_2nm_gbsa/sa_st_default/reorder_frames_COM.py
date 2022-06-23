import mdtraj as md
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--traj_file_list",
        required = True,
        type = str,
        nargs = "+"
    )

    parser.add_argument(
        "--top",
        required = True,
        type = str,
    )
    parser.add_argument(
        "--out_name",
        required = True,
        type = str
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    print(args.top)
    traj = md.load(args.traj_file_list, top = args.top)
    
    # Select both octomers
    topology = traj.topology
    selection = topology.select("resname OCT")

    # Atom Slice atoms of interest
    traj_dimer = traj.atom_slice(selection)
    top_dimer = traj_dimer.topology

    # Reorder trajectory based on distance between two chains
    res_0 = top_dimer.select("residue 1")
    res_1 = top_dimer.select("residue 2")

    res_0 = traj_dimer.atom_slice(res_0)
    res_1 = traj_dimer.atom_slice(res_1)

    com_distances = []
    indices = []
    for i in range(traj_dimer.n_frames):
        com_0 = md.compute_center_of_mass(res_0[i])
        com_1 = md.compute_center_of_mass(res_1[i])
        d = com_0 - com_1
        d = np.sqrt(np.dot(d,d.T))[0][0]
        
        com_distances.append(d)
        indices.append(i)

    sorted_indices = [x for _, x in sorted(zip(com_distances, indices))]

    sorted_traj = traj_dimer[sorted_indices]

    sorted_traj[0].save(args.out_name + ".pdb")
    sorted_traj.save_xtc(args.out_name + ".xtc")





if __name__ == "__main__":
    main()