import mdtraj as md
import numpy as np
import argparse
import panedr
import pandas
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import os

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--directories",
        help = "directories to search for specified file names. Example sim{0..63}",
        nargs='+'
    )
    parser.add_argument(
        "--top",
        help = "topology file to use when loading in trajectory files"
    )

    parser.add_argument(
        "--base_file_name",
        help = "basename of .edr, .trr, .log, files to search for in each directory",
        nargs='+'
    )

    parser.add_argument(
        "--traj_format",
        help = "format of file to search for. Default is xtc",
        default = "xtc"
    )

    return(parser.parse_args())


def get_COM_distance(traj_object, selection1, selection2, time_per_frame = 0.1):
    # Get COM1 Trajectory
    top = traj_object.topology
    sel_1 = top.select(selection1)
    sel_1_traj = traj_object.atom_slice(sel_1)
    com_1 = md.compute_center_of_mass(sel_1_traj)

    # Center coordinates around COM1
    print("nsteps:", com_1.shape[0])
    for i_frame in range(com_1.shape[0]):
        traj_object.xyz[i_frame, :, :] = traj_object.xyz[i_frame, :, :] - com_1[i_frame]
    # Center coordinates about a single molecule
    # com_1 = md.compute_center_of_mass()
    sel_2 = top.select(selection2)
    sel_2_traj = traj_object.atom_slice(sel_2)
    com_2 = md.compute_center_of_mass(sel_2_traj)
    box_lenghts = traj_object.unitcell_lengths

    print(com_2.shape[0])
    d_com = np.zeros(com_2.shape[0])
    for i_frame in range(com_2.shape[0]):
        # apply PBCs to COM_2
        com_2[i_frame] = com_2[i_frame] - box_lenghts[i_frame] * np.floor(com_2[i_frame]/(box_lenghts[i_frame]) + 0.5)
        d_com[i_frame] = np.sqrt(np.dot(com_2[i_frame], com_2[i_frame]))

    time = np.array(range(com_1.shape[0]))*time_per_frame
    return(time, d_com)

def get_edr_obs(edr_df, key):
    time = edr_df["Time"]*0.002
    observable = edr_df[key]
    return(time, observable)


def main():
    args = parse_arguments()

    sim_obs = {}
    for sim_dir in args.directories:
        all_files  = os.listdir(sim_dir)
        for file_id in args.base_file_name:
            sim_files = [fn for fn in all_files if file_id in fn]
            
            # Load trajectory file
            traj_files = [fn for fn in sim_files if args.traj_format in fn]
            traj = None
            traj_files.sort() # Need partXXXX on files to ensure order is correct here
            for traj_file in traj_files:
                print("Loading:", os.path.join(sim_dir,traj_file))
                if traj is None:
                    traj = md.load(os.path.join(sim_dir,traj_file), top = os.path.join(sim_dir, args.top))
                else:
                    temp_traj =  md.load(os.path.join(sim_dir,traj_file), top = os.path.join(sim_dir, args.top))
                    traj = traj.join(temp_traj)

            # Load edr file
            edr_files = [fn for fn in sim_files if ".edr" in fn]
            edr_files.sort()
            energy_df = None
            for edr_file in edr_files:
                print("Loading:", os.path.join(sim_dir,edr_file))
                if energy_df is None:
                    energy_df = panedr.edr_to_df(os.path.join(sim_dir, edr_file))
                else:
                    energy_df = pd.concat([energy_df, panedr.edr_to_df(os.path.join(sim_dir, edr_file))])
            energy_df.sort_values(by="Time")


            # observables of interest
            # 'Time', 'Bond', 'Angle', 'Proper Dih.', 'LJ-14', 'Coulomb-14',
            # 'LJ (SR)', 'Coulomb (SR)', 'COM Pull En.', 'Potential', 'Kinetic En.',
            # 'Total Energy', 'Temperature', 'Pressure', 'Constr. rmsd', 'Vir-XX',
            #'Vir-XY', 'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', 'Vir-ZY',
            #'Vir-ZZ', 'Pres-XX', 'Pres-XY', 'Pres-XZ', 'Pres-YX', 'Pres-YY',
            #'Pres-YZ', 'Pres-ZX', 'Pres-ZY', 'Pres-ZZ', '#Surf*SurfTen',
            #'T-System'

            time, com = get_COM_distance(traj, "residue 1", "residue 2")
            te_time, tot_energy = get_edr_obs(energy_df, "Total Energy")
            pe_time, pull_energy = get_edr_obs(energy_df, "COM Pull En.")
            sim_obs[sim_dir] = {"com"    : [time, com],
                                "tot_energy" : [te_time, tot_energy],
                                "pull_energy" : [pe_time, pull_energy]
                                }
    
    # Plot all simulation observables on one plot
    # com plot
    fig, ax = plt.subplots(figsize = [5,5], dpi = 300)
    for key in sim_obs.keys():
        ax.plot(sim_obs[key]["com"][0], sim_obs[key]["com"][1], lw = 0.5, alpha = 0.5)
    ax.set_xlabel(r'Time (ns)')
    ax.set_ylabel(r'$d_{COM}$ (nm)')
    plt.tight_layout()
    fig.savefig("COM_plot.jpg")

    fig, ax = plt.subplots(figsize = [5,5], dpi = 300)
    for key in sim_obs.keys():
        ax.plot(sim_obs[key]["tot_energy"][0], sim_obs[key]["tot_energy"][1], lw = 0.5, alpha = 0.5)
    ax.set_xlabel(r'Time (ns)')
    ax.set_ylabel(r'Total Energy (kJ mol$^{-1}$)')
    plt.tight_layout()
    fig.savefig("energy_plot.jpg")

    fig, ax = plt.subplots(figsize = [5,5], dpi = 300)
    for key in sim_obs.keys():
        ax.hist(sim_obs[key]["tot_energy"][1], bins = 50)
    ax.set_xlabel(r'Time (ns)')
    ax.set_ylabel(r'Total Energy (kJ mol$^{-1}$)')
    plt.tight_layout()
    fig.savefig("energy_hist.jpg")

    fig, ax = plt.subplots(figsize = [5,5], dpi = 300)
    for key in sim_obs.keys():
        ax.plot(sim_obs[key]["pull_energy"][0], sim_obs[key]["pull_energy"][1], lw = 0.5, alpha = 0.5)
    ax.set_xlabel(r'Time (ns)')
    ax.set_ylabel(r'Total Energy (kJ mol$^{-1}$)')
    plt.tight_layout()
    fig.savefig("pull_energy.jpg")


if __name__ == "__main__":
    main()
