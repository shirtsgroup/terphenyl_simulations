import tqdm
import os
import itertools
import mdtraj
import numpy as np
from multiprocessing import Pool


class REMD_trajectories:
    def __init__(self, path, traj_id, traj_file_type, sim_id, top, np = 4):
        self.path = path
        self.top = top
        self.sim_id = sim_id
        self.n_replicas = len([a for a in os.listdir(self.path) if a.startswith(sim_id)])
        self.trajs = []
        self.temps = []
        self.traj_id = traj_id
        self.traj_file_type = traj_file_type

        pool = Pool(np)
        # read-in replica data
        self.trajs, self.temps = zip(*pool.map(self.read_replica, range(self.n_replicas)))


    def read_replica(self, i):
        print("Replica", i)
        with open(os.path.join(self.path, self.sim_id + str(i), "npt.mdp"), "r")as f:
            for line in f.readlines():
                if "ref-t" in line:
                    T = float(line.split()[-1])
        traj_files = []
        all_files = os.listdir(os.path.join(self.path, self.sim_id + str(i),))
        all_files.sort()
        all_files.remove(self.traj_id + "." + self.traj_file_type)
        all_files.insert(0, self.traj_id + "." + self.traj_file_type)
        
        for file in all_files:
            if file.startswith(self.traj_id) and file.endswith(self.traj_file_type):
                print(file)
                traj_files.append(os.path.join(self.path, self.sim_id + str(i),file))
        traj = mdtraj.load(traj_files, top=self.top)
        return(traj, T)

class HydrogenBondFinder:
    def __init__(self, ref_frame, topology, h_acceptor_element = "O",
                 h_donor_element = "N", residue_selector = "and resname OCT",
                 h_donor_bond_length = 0.105, h_bond_length = 0.3, h_bond_angle = 120):
        self.ref_frame = ref_frame
        self.top = topology
        self.h_acceptor_element = h_acceptor_element
        self.h_donor_element = h_donor_element
        self.h_donor_bond_length = h_donor_bond_length
        self.selection_type = "element"
        self.residue_selector = residue_selector
        self.h_bond_length = h_bond_length
        self.h_bond_angle = h_bond_angle
        self.donors = []
        self.acceptors = []
    
    def get_donors(self):
        heavy_atoms = self.top.select(" ".join([self.selection_type, self.h_donor_element, self.residue_selector]))
        hydrogens = self.top.select(" ".join(["element H", self.residue_selector]))
        for pair in itertools.product(heavy_atoms, hydrogens):
            r1 = self.ref_frame.xyz[0, pair[0], :] - self.ref_frame.xyz[0, pair[1], :]
            if np.dot(r1, r1) < self.h_donor_bond_length**2:
                self.donors.append(pair)
    
    def get_acceptors(self):
        self.acceptors = self.top.select(" ".join([self.selection_type, self.h_acceptor_element, self.residue_selector]))

    def get_hydrogen_bonds(self, trajectory):
        n_h_bonds = []
        h_bonds = []
        for i in tqdm.tqdm(range(trajectory.xyz.shape[0])):
            frame_n_h_bond = 0
            frame_h_bond = []
            for pair in itertools.product(self.donors, self.acceptors):
                r_hbond = trajectory.xyz[i, pair[0][1], :] - trajectory.xyz[i, pair[1], :]
                r_hb_2 = np.dot(r_hbond, r_hbond)
                if r_hb_2 < self.h_bond_length**2:
                    r1 = trajectory.xyz[i, pair[0][0], :] - trajectory.xyz[i, pair[0][1], :]
                    theta = np.arccos(np.dot(r1, r_hbond) / np.sqrt(r_hb_2 * np.dot(r1, r1)))
                    if np.abs(theta - np.pi) > (np.pi / 180 * self.h_bond_angle):
                        frame_n_h_bond += 1
                        frame_h_bond.append([pair[0], pair[1]])
            n_h_bonds.append(frame_n_h_bond)
            h_bonds.append(frame_h_bond)
        return(n_h_bonds, h_bonds)
            
def main():
    import argparse
    import plotting
    import matplotlib.pyplot as plt
    from multiprocessing import Pool
    import heat_capacity
    import time
    import pymbar

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sim_dir_name",
        required = True,
        help = "base name of replica directories",
        type = str
    )
    parser.add_argument(
        "--path",
        required = True,
        help = "path to RE exchange simulation directories"
    )
    parser.add_argument(
        "--n_replicas",
        required = True,
        help = "number of replicas run in replica exchange simulation",
        type = int
    )
    parser.add_argument(
        "--sim_id",
        help = "simulation name ID used in replica exchange simulation",
        type = str,
        default = "npt"
    )
    parser.add_argument(
        "--figure_id",
        required = True,
        help = "base name for figure generated in this script",
        type = str
    )
    parser.add_argument(
        "--n_estimated_points",
        help = "number of points to estimate between maximum and minimum temperature",
        type = int,
        default = 200
    )

    args = parser.parse_args()

    sim_dir_name = args.sim_dir_name
    path = args.path
    n_replicas = args.n_replicas
    sim_id = args.sim_id

    print("Reading Energy Files...")
    energies, temps = heat_capacity.get_energies(sim_dir_name, path, n_replicas, sim_name = sim_id)

    # Plot energy distributions of RE simualtion
    plotting.plot_RE_energy_distributions(energies)
    plt.savefig(args.figure_id + "energy_dist.png")

    # Plot energy trajectory of RE simulation
    plotting.plot_RE_energy_trajectory(energies)
    plt.savefig(args.figure_id + "energy_traj.png")

    u_kln, n_samples, t_list, betas = heat_capacity.construct_u_kln_matrix(temps, energies, add_temps = np.linspace(min(temps), max(temps), args.n_estimated_points))

    remd_trajs = REMD_trajectories("/mnt/summit/simulations/octamer_Rchiral/RHH/remd_sim/200K_to_350K", "npt", "whole.xtc","sim", "/mnt/summit/simulations/octamer_Rchiral/RHH/remd_sim/200K_to_350K/sim0/berendsen.gro", np = 1)

    h_bond_finder = h_bonds.HydrsogenBondFinder(remd_trajs.trajs[0][0], remd_trajs.trajs[0][0].top)
    h_bond_finder.get_donors()
    h_bond_finder.get_acceptors()
    n_h_bonds, h_bond_ids = h_bond_finder.get_hydrogen_bonds(remd_trajs.trajs[0])

    pool = Pool(4)

    n_h_bonds_remd, h_bonds_remd = zip(*pool.map(h_bond_finder.get_hydrogen_bonds, remd_trajs.trajs))
    n_h_bonds_remd = np.array(n_h_bonds_remd)

    h_bonds_kln = np.zeros(len(betas), len(betas), n_h_bonds_remd[1])
    for k in range(n_h_bonds_remd.shape[0]):
        for l in range(h_bonds_kln.shape[0]):
            h_bonds_kln[k,l,:] = n_h_bonds_remd[k]

    mbar_h_bonds = mbar_h_bonds = pymbar.MBAR(u_kln[:,:,:3608], n_samples_h_bonds, verbose = True, relative_tolerance = 1e-10, initial_f_k= None, maximum_iterations=1000)