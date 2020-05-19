#! ~/anaconda3/envs/python36/bin/python

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import extract_1D_dihedrals
import argparse
from utils import ReadItpFile

matplotlib.rcParams.update({'font.size': 30})
matplotlib.rcParams.update({'font.family': "serif"})
plt.margins(x=0)


# dihedrals of interest
# Crystal monomer:
# Dihedral 1 : 143.962
# Dihedral 2: -6.792
# Dihedral 3: 124.974

dihe_nums = ['1', '2', '3', '4', '5', '6', '7']
dihe_num_to_annotation = {1:"A", 2:"B", 3:"CC", 4:"D", 5:"E", 6:"F", 7:"G"}


# to investigate: Time steps after 3500


# def get_md_dihedrals(dihe_id, output_file):
#     extract_1D_dihedrals_PR.read_in


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file_id",
                        help="ID used to label figure and npy file output",
                        required=True,
                        type=str
                        )
    parser.add_argument("--dihedral_ids",
                        help="Ids of which dihedral you would like to generate",
                        required=False,
                        nargs='+',
                        type=int,
                        )
    parser.add_argument("--struct_file",
                        help="location of structure file to load into MDtraj",
                        required=True,
                        type=str,
                        )
    parser.add_argument("--traj_files",
                        help="path to trajectory files to load into MDtraj",
                        required=True,
                        nargs='+',
                        type=str,
                        )
    parser.add_argument("--itp_path",
                        help="path to the itp file for a trajectory",
                        required=True,
                        type=str,
    )

    parser.add_argument("--scatter",
                        help="datasets to scatter on plots",
                        default= None,
                        nargs="+",
                        required=False,
    )

    args = parser.parse_args()

    # Generate/Load dihedrals from minimized structure


    color = cm.get_cmap("Set1", lut=len(args.scatter))
    extra_points = []
    for points in args.scatter:
        extra_points.append(np.load(points, allow_pickle = True))
    mm2_gellman_dihes = np.load("mm2_torsions_gellman.npy", allow_pickle = True)
    mmf94_gellman_dihes = np.load("mmf94_torsions_gellman.npy", allow_pickle = True)
    # crystal_monomer = [143.962, -6.792, 124.974]

    if "outputs" not in args.file_id:
        args.file_id = "outputs/" + args.file_id

    dihe_nums = args.dihedral_ids
    # Generate dataset
    fig_names = []
    traj = extract_1D_dihedrals.read_in_trajectory(
        args.struct_file, args.traj_files)
    for dihe_num in dihe_nums:
        figname = args.file_id+"_"+str(dihe_num)+".npy"
        fig_names.append(figname)
        if os.path.exists(figname):
            continue
        else:
            annote_list = [dihe_num_to_annotation[dihe_num]+str(i) for i in range(1,5)]
            itp =  ReadItpFile(args.itp_path)
            selector = itp.construct_dihe_selection(annote_list)
            dihes = extract_1D_dihedrals.get_dihedrals(
                'name ' + selector, traj)*180/np.pi
            print(selector)
            np.save(args.file_id+"_"+str(dihe_num)+".npy", dihes)

    index_to_axis = {1: r'Dihedral Angle 1',
                     2: r'Dihedral Angle 2',
                     3: r'Dihedral Angle 3',
                     4: r'Dihedral Angle $\theta$',
                     5: r'Dihedral Angle $\phi$',
                     6: r'Dihedral Angle $\omega$',
                     7: r'Dihedral Angle $\psi$',
                     }

    fig, ax = plt.subplots(nrows=len(args.dihedral_ids),
                           ncols=len(args.dihedral_ids), figsize=[40, 40])
    # use LaTeX fonts in the plot
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    
    
    for i in range(len(args.dihedral_ids)):
        for j in range(len(args.dihedral_ids)):
            if i == j:
                # 1D distribution
                dist = np.load(fig_names[i])
                ax[i, j].hist(dist, density=True,
                                  bins=60, range=[-180, 180])
                ax[i,j].set_xlim([-180, 180])
                ax[i, j].set_xlabel(index_to_axis[args.dihedral_ids[i]])
                ax[i, j].set_ylabel('Probability Density')
                c_i = 0
                for extra_point in extra_points:
                    for each_point in extra_point[args.dihedral_ids[j] - 1]:
                        ax[i, j].axvline(each_point, ls = "--", color=color(c_i))
                    c_i += 1

            elif i > j:
                fig.delaxes(ax[i, j])
            else:
                # 2D distribution

                dist1 = np.load(fig_names[i])
                dist2 = np.load(fig_names[j])
                print(dist1)
                print(dist2)
                print("dist 1:", len(dist1), "dist 2:", len(dist2))
                print("dist 1:", dist1.dtype, "dist 2:", dist2.dtype)
                print(fig_names[i])
                print(fig_names[j])
                ax[i, j].hist2d(dist2, dist1, bins=60, density=True, range=[[-180, 180], [-180, 180]])
                ax[i,j].set_xlim([-180, 180])
                ax[i,j].set_ylim([-180, 180])
                # ax[i, j].set_xlabel(index_to_axis[args.dihedral_ids[j]])
                # ax[i, j].set_ylabel(index_to_axis[args.dihedral_ids[i]])
                # ax[i-1, j-1].scatter(crystal_monomer[i-1], crystal_monomer[j-1], color='red')
                c_i = 0
                for extra_point in extra_points:
                    print(extra_point)
                    print(args.dihedral_ids[j] - 1)
                    ax[i, j].scatter(extra_point[args.dihedral_ids[j] - 1],
                                        extra_point[args.dihedral_ids[i] - 1],  s=50, color=color(c_i))
                    c_i += 1

    handles, labels = ax[1,1].get_legend_handles_labels()
    print(labels)
    fig.legend(handles, labels, loc='lower left')
    fig.savefig(args.file_id+'_torsions.pdf')
    fig.savefig(args.file_id+'_torsions.png')


if __name__ == '__main__':
    main()
