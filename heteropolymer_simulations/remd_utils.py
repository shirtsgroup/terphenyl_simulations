import mdtraj as md
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import heteropolymer_simulations as hs

plt.rcParams.update({"font.size": 7})


class REMDTrajectory:
    def __init__(self, directories, traj_files, top_file):
        self.directories = directories
        self.traj_files = traj_files
        self.top_file = top_file

        # get trajectories into MDTraj Trajectory object
        self.remd_trajs = []
        for d in self.directories:
            print("Loading files in", d)
            traj = md.load(
                [os.path.join(d, tf) for tf in traj_files],
                top=os.path.join(d, self.top_file),
            )
            self.remd_trajs.append(traj)

    def get_continuous_trajectories(self, remd_log, nstout, output_dir="continuous"):

        # Create output directory
        if os.path.isdir(output_dir):
            hs.clustering.backoff_directory(output_dir)
        os.mkdir(output_dir)

        # Get the states where frames were written
        steps = np.array(remd_log.steps)
        nst_index = np.where(np.mod(steps, nstout) == 0)
        out_states = [
            np.array(state_traj)[nst_index] for state_traj in remd_log.state_trajs
        ]

        # Get frame from correct state and create new continuous traj with
        # indices from specific replica
        const_traj = []
        for i in range(remd_log.n_states):
            xyz_cont = np.zeros(self.remd_trajs[i].xyz.shape)
            print("remd_traj[i].xyz.shape:", self.remd_trajs[i].xyz.shape)
            print(len(out_states[i]))
            for j, state in enumerate(out_states[i]):
                if j >= self.remd_trajs[i].xyz.shape[0]:
                    break
                xyz_cont[j, :, :] = self.remd_trajs[state].xyz[j, :, :]
            traj = md.Trajectory(xyz_cont, self.remd_trajs[i].topology)
            const_traj.append(traj)

        print("Writing output files...")

        for i in range(len(const_traj)):
            filename = os.path.join(output_dir, "replica_" + str(i) + ".xtc")
            const_traj[i].save(filename)


class REMDLogFile:
    def __init__(self, log_file):
        self.log_file = log_file
        self.get_state_trajectory()

    def get_state_trajectory(self):
        order_traj = []
        steps = []
        times = []
        with open(self.log_file, "r") as lf:
            for line in lf.readlines():
                if "Order After Exchange:" in line:
                    order = [int(a) for a in line.split()[3:]]
                    order_traj.append(order)
                if "Replica exchange at step" in line:
                    step = int(line.split()[4])
                    time = float(line.split()[6])
                    steps.append(step)
                    times.append(time)

        replicas = np.unique(order_traj[0])
        self.n_states = len(replicas)
        state_trajs = []
        for i in range(len(replicas)):
            state_traj = []
            state_traj = [o.index(i) for o in order_traj]
            state_trajs.append(state_traj)

        self.state_trajs = state_trajs
        self.steps = np.array(steps)
        self.times = np.array(times) / 1000

    def plot_state_trajectory(self):
        plot_dimensions = int(np.ceil(np.sqrt(self.n_states)))
        fig, axs = plt.subplots(
            plot_dimensions,
            plot_dimensions,
            figsize=[plot_dimensions * 2, plot_dimensions * 2],
        )
        for i, ax in enumerate(fig.axes):
            ax.plot(self.times[::], self.state_trajs[i], linewidth=0.5)
            ax.set_xlabel("Simulation Time (ns)")
            ax.set_ylabel("State Index")
            ax.set_title("Replica " + str(i))
        fig.tight_layout()
        fig.savefig("test.png", dpi=300)


def main():
    t1 = time.time()
    log_object = REMDLogFile(
        "/ocean/projects/cts160011p/tfobe/heteropolymer_simulations/simulations/terphenyl_mop/tetramer_remd/t_200_350/sim0/npt_new.log"
    )
    prefix = "/ocean/projects/cts160011p/tfobe/heteropolymer_simulations/simulations/terphenyl_mop/tetramer_remd/t_200_350/"
    directories = [prefix + "sim" + str(i) for i in range(64)]
    traj_files = ["npt_new.whole.xtc"]
    traj_object = REMDTrajectory(directories, traj_files, "berendsen_npt.gro")
    traj_object.get_continuous_trajectories(log_object, 50000)

    t2 = time.time()

    print("This analysis took", round(t2 - t1, 2), "second(s).")


if __name__ == "__main__":
    main()
