#! ~/anaconda3/envs/python36/bin/python

import mdtraj as md
import numpy as np
import panedr
import os
from pymbar import timeseries


# Load in files

def read_in_trajectory(structure, trajectories):
    traj = md.load(trajectories, top=structure)
    return(traj)


def read_in_energy(edr_file):
    energy = panedr.edr_to_df(edr_file)
    return(energy)


def get_equilibrium_frames(energy):
    t0, g, Neff_max = timeseries.detectEquilibration(
        energy['Potential'].values)
    return(t0)

# Dihe of interest : C3, C4, C6, C7
# Dihe of interest : C4, C6, C7, C12
# Dihe of interest : C6, C7, C12, C13
# Dihe of interest : C6, C7, C12, C13


def construct_selector(atom_list, res_dict, n_residues):
    """
    Since our molecule is all in one residue, we need some interesting ways to extract residues, instead of a simple
    select string

    atom_list : list
        String name of base atoms in the first residue
    res_dict : dict
        Dictionary telling how many atoms of a given element are in each residue
    """
    atom_info = [[atom[0], int(atom[1:])] for atom in atom_list]

    selector_list = []
    for i in range(n_residues):
        atomnos = [atom[0]+str(atom[1] + i*res_dict[atom[0]])
                   for atom in atom_info]
        # print(atomnos)
        for j in range(len(atomnos)):
            selector_list.append(atomnos[j])

    selector = ' '.join(selector_list)
    return(selector)


def get_dihedrals(selector, traj):
    print(selector)
    top = traj.topology
    i_dihes = []
    for atom in selector.split(' ')[1:]:
        i_dihe = top.select("name "+atom)
        i_dihes.append(i_dihe)
    i_dihes = np.array(i_dihes)
    i_dihes = i_dihes.reshape(int(len(i_dihes)/4), 4)
    dihes = md.compute_dihedrals(traj, i_dihes, periodic=True).flatten()
    print(i_dihe)
    return(dihes)


def main():
    import matplotlib.pyplot as plt

    traj_files = ['../octamer/initial_configs/1/PR.trr',
                  '../octamer/initial_configs/3/PR.trr',
                  '../octamer/initial_configs/4/PR.trr',
                  '../tetramer/PR.trr']
    struct_files = ['../heteropolymer_inputs/octamer/initial_configs/frame_433.gro',
                    '../heteropolymer_inputs/octamer/initial_configs/frame_618.gro',
                    '../heteropolymer_inputs/octamer/initial_configs/solvated.gro',
                    '../heteropolymer_inputs/tetramer/solvated.gro']
    energy_files = ['../octamer/initial_configs/1/PR.edr',
                    '../octamer/initial_configs/3/PR.edr',
                    '../octamer/initial_configs/4/PR.edr',
                    '../tetramer/PR.edr']
    fig_names = ['outputs/octamer_config_1_dihe_',
                 'outputs/octamer_config_3_dihe_',
                 'outputs/octamer_config_4_dihe_',
                 'outputs/tetramer_dihe_']
    selectors = [['C23', 'C24', 'C26', 'C27'], ['C24', 'C26', 'C27', 'C32'], ['C26', 'C27', 'C32', 'C33'], [
        'C14', 'C15', 'C19', 'N1'], ['C15', 'C19', 'N1', 'C38'], ['C19', 'N1', 'C38', 'C21'], ['N1', 'C38', 'C21', 'C20']]

    chain_sizes = [6, 6, 6, 2]

    dihe_nums = ['1', '2', '3', '4', '5', '6', '7']

    for traj_file, struct_file, energy_file, fig_name, chain_size in zip(traj_files, struct_files, energy_files, fig_names, chain_sizes):

        if not np.all([os.path.exists(f) for f in [fig_name+dihe_num+'.npy' for dihe_num in dihe_nums]]):
            traj = read_in_trajectory(struct_file, traj_file)  # expensive
            energy = read_in_energy(energy_file)
            t0 = get_equilibrium_frames(energy)

        # Dihe of interest : C3, C4, C6, C7
        # Dihe of interest : C4, C6, C7, C12
        # Dihe of interest : C6, C7, C12, C13
        for sel, dihe_num in zip(selectors, dihe_nums):
            selector = construct_selector(
                sel, {'C': 20, 'N': 1, 'O': 1}, chain_size)
            if os.path.exists(fig_name+dihe_num+'.npy'):
                dihes = np.load(fig_name+dihe_num+'.npy')
            else:
                dihes = get_dihedrals(
                    'name '+selector, traj[t0:])*180/np.pi  # expensive
                np.save(fig_name+dihe_num+'.npy', dihes)
            plt.figure()
            plt.hist(dihes, bins=30, density=True)
            plt.ylabel('Probability Density')
            plt.xlabel('Dihedral Angle')
            plt.savefig(fig_name+dihe_num+'.pdf')
            plt.savefig(fig_name+dihe_num+'.png')
            plt.close()

            # plt.show()


if __name__ == '__main__':
    main()
