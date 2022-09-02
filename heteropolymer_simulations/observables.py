import mdtraj as md
import numpy as np
import argparse
import panedr
import pandas
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import os

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

    d_com = np.zeros(com_2.shape[0])
    for i_frame in range(com_2.shape[0]):
        # apply PBCs to COM_2
        com_2[i_frame] = com_2[i_frame] - box_lenghts[i_frame] * np.floor(com_2[i_frame]/(box_lenghts[i_frame]) + 0.5)
        d_com[i_frame] = np.sqrt(np.dot(com_2[i_frame], com_2[i_frame]))

    time = np.array(range(com_1.shape[0]))*time_per_frame
    return(time, d_com)

def get_ROG_distance(traj_object, selection, time_per_frame = 0.1):
    # Get COM1 trajectory
    top = traj_object.topology
    sel = top.select(selection)
    sel_traj = traj_object.atom_slice(sel)
    com_1 = md.compute_center_of_mass(sel_traj)

    # Center coordinates around COM1
    print("nsteps:", com_1.shape[0])
    for i_frame in range(com_1.shape[0]):
        sel_traj.xyz[i_frame, :, :] = sel_traj.xyz[i_frame, :, :] - com_1[i_frame]
    
    rg = md.compute_rg(sel_traj)
    time = np.array(range(com_1.shape[0]))*time_per_frame
    return(time, rg)
                

def get_edr_obs(edr_df, key):
    time = edr_df["Time"]*0.002
    observable = edr_df[key]
    return(time, observable)


def get_torsions(traj_obj, torsion_atom_names, mirror_sym = False):
    top = traj_obj.topology
    torsions_inds = []
    for torsion_atoms in torsion_atom_names:
        torsion_i =[top.select("name " + atom)[0] for atom in torsion_atoms]
        torsions_inds.append(torsion_i)

    # torsions_inds = np.array(torsions_inds)
    torsions = md.compute_dihedrals(traj_obj, torsions_inds, periodic=False)
    torsions = torsions.reshape(-1)
    if mirror_sym:
        for i in range(len(torsions)):
            if torsions[i] <= 0:
                torsions[i] = torsions[i] + np.pi
    return torsions