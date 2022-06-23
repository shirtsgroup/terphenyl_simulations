import mdtraj as md
import numpy as np
import argparse
import panedr
import matplotlib.pyplot as plt
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

    return parser.parse_args()

def plot_energies(edr_df):
    pass

def plot_COM_distance(traj_object, selection1, selection2):
    
    # Get COM1 Trajectory
    top = traj_object.topology
    sel_1 = top.select(selection1)
    sel_1_traj = traj_object.atom_slice(sel_1)
    com_1 = md.compute_center_of_mass(sel_1_traj)
    print(com_1)

    # Center coordinates around COM1 and apply PBCs
    print(traj_object.xyz[0])
    traj_object.xyz = traj_object.xyz - com_1
    print(traj_object.xyz[0])

    traj_object.save_xtc("test_output.xtc")

    # Center coordinates about a single molecule
    com_1 = md.compute_center_of_mass()

    com_distances = []

    

def plot_density(traj_object):
    pass

def plot_pressure(edr_df):
    pass



def main():
    args = parse_arguments()

    for sim_dir in args.directories:
        all_files  = os.listdir(sim_dir)
        for id in args.base_file_name:
            sim_files = [fn for fn in all_files if id in fn]
            
            # Load trajectory file
            print(args.traj_format)
            traj_files = [fn for fn in sim_files if args.traj_format in fn]
            traj = md.load(os.path.join(sim_dir, args.top))
            for traj_file in traj_files:
                traj.join(md.load(traj_file, top = args.top))
            
            # Load edr file
            edr_files = [fn for fn in sim_files if ".edr" in fn]
            df = pd.DataFrame()
            for edr_file in traj_files:
                df.append(panedr.edr_to_df(edr_file))
            
            # Plotting
            plot_COM_distance(traj, "residue 1", "residue 2")




            

            


if __name__ == "__main__":
    main()