import terphenyl_simulations as ts
import mdtraj as md
import numpy as np
from tqdm import tqdm
from glob import glob
import sys
import csv
from subprocess import Popen, PIPE
import shutil
import matplotlib.pyplot as plt


def main():

    ts.utils.make_path("demux")

    traj_output = 50000
    state_output = 100
    freq = int(traj_output/state_output)

    gromacs_log = ts.utils.GromacsLogFile("sim0/npt_new.log")
    # Add a state for frame 1, where replicas start at their index
    replica_ids = list(range(max(gromacs_log.states[0]) + 1))
    gromacs_log.states.insert(0, replica_ids)

    exchange_states = len(gromacs_log.states)
    print("The simulation states written", exchange_states, "times.")

    demux_indices = []
    for state_i in np.unique(gromacs_log.states):
        state_traj = []
        for state in gromacs_log.states:
            state_traj.append(state.index(state_i))
        demux_indices.append(state_traj)
    
    print("Reading REMD trajectories...")
    remd_trajs = [md.load("sim" + str(i) + "/npt_new.xtc", top = "sim0/berendsen_npt.gro") for i in tqdm(range(len(glob("sim*"))))]
    n_frames = len(remd_trajs[0])
    tet_sele = remd_trajs[0].topology.select("resn TET or resn CAP")
    print("The simulation has", n_frames, "frames.")

    plt.plot(demux_indices[0])
    plt.savefig("test_state_traj.png")

    # Build demuxed trajectories for each replica
    # for rep_k, demux_i in enumerate(demux_indices):
    #     print("Working on Replica:", rep_k)
    #     # Makes an empty trajectory that can have individual frames joined to it
    #     new_traj = md.load("sim" + str(rep_k) + "/berendsen_npt.gro", top = "sim" + str(rep_k) + "/berendsen_npt.gro")[1:]
    #     
    #     for frame_i, replica_j in enumerate(demux_i[::freq]):
    #         print("Frame", frame_i)
    #         new_traj = md.join([new_traj, remd_trajs[replica_j][frame_i]])
    # 
    # 
    #     # print("Demux trajectory has", len(new_traj), "frames.")
    #     new_traj.save_xtc("demux/replica_" + str(rep_k) + ".xtc")

    # Using RMSDs to determine closest frame
    # This is really slow, every frame it has
    # to check RMSD to all other simulations
    # Could be improved by calculating RMSD matrices
    rmsd_demux = [ [] for _ in range(n_frames)]
    for rep_i in replica_ids:
        print("Working on replica", rep_i, "...")
        new_traj = md.load("sim" + str(rep_i) + "/berendsen_npt.gro", top = "sim" + str(rep_i) + "/berendsen_npt.gro")

        for frame_i in tqdm(range(n_frames)):
            remd_frames_i = md.join([remd_traj[frame_i] for remd_traj in remd_trajs])
            rmsds = list(md.rmsd(remd_frames_i, new_traj[-1], precentered = True))
            # print(rmsds)
            min_rmsd_i = rmsds.index(min(rmsds))
            new_traj = md.join([new_traj, remd_trajs[min_rmsd_i][frame_i]])
            rmsd_demux[frame_i].append(min_rmsd_i)

        filename = "demux/replica_" + str(rep_i) + ".xtc"
        new_traj.save_xtc(filename)
        if shutil.which("gmx"):
            p = Popen(["gmx", "trjconv", "-f", filename, "-s", "sim0/npt_new.tpr", "-pbc", "whole", "-o", "demux/replica_" + str(rep_i) + ".whole.xtc"], stdin = PIPE, stdout = PIPE)
            p.communicate(input = b'0\n')
            
    with open("rmsd_demux.csv", "wb") as f:
        writer = csv.writer(f)
        writer.writerows(rmsd_demux)

    
    


if __name__ == "__main__":
    main()