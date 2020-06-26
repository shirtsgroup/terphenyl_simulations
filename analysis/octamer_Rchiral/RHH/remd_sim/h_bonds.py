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
            


