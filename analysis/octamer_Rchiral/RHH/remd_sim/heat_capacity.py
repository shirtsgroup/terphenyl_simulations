import os
import numpy as np
from scipy.constants import physical_constants
import mdtraj
import panedr
from multiprocessing import Pool
import time
import pymbar

kb = physical_constants["Boltzmann constant"][0] *  physical_constants["Avogadro constant"][0] / 1000 # J (molK)^-1

def get_energies(sim_dir_name, path, n_replicas, sim_name = "npt", energy_type = "Enthalpy", return_T = True):
    energies = []
    temps = []
    for i in range(n_replicas):
        sim_dir = os.path.join(path, sim_dir_name + str(i))
        T = None
        with open(os.path.join(sim_dir, sim_name + ".mdp"), "r")as f:
            for line in f.readlines():
                if "ref-t" in line:
                    T = float(line.split()[-1])
                    temps.append(T)
        sim_energies = []
        files = os.listdir(sim_dir)
        files.sort()
        for file in files:
            if file.startswith(sim_name) and file.endswith(".edr"):
                print(file)

                edr_file = os.path.join(sim_dir, file)
                energy_df = panedr.edr_to_df(edr_file)
                pot_energy = energy_df[energy_type].values
                sim_energies += list(pot_energy)
        sim_energies = np.array(sim_energies)
        energies.append(sim_energies)
    temps = np.array(temps)
    if return_T:
        return(energies, temps)
    else:
        return(energies)


def construct_u_kln_matrix(t_list, energies, add_temps = np.linspace(230, 325, 200)):
    # Initialize MBAR inputs
    K_samples = len(t_list) # number of states
    n_samples = [len(energies[i]) for i in range(len(energies))] # number of samples from each state

    # Additional temperatures
    if add_temps is not None:
        t_list = list(t_list)
        for T in add_temps:
            t_list.append(T)
            n_samples.append(0)
        t_list = np.array(t_list)

    betas = 1 / kb / np.array(t_list)
    K_all = len(t_list)
    n_samples = np.array(n_samples)



    u_kln = np.zeros([K_all, K_all, np.max(n_samples)], np.float64)

    for k in range(K_samples):
        for l in range(K_all):
            u_kln[k,l,:len(energies[k])] = betas[l] * energies[k]

    return u_kln, n_samples, t_list, betas


def calculate_cp(mbar, E_kln, t_list):
    
    # mean and variance of energy
    results = mbar.computeExpectations(E_kln, state_dependent=True, return_dict=True)
    E_expected = results['mu']
    dE_exoected = results['sigma']
    
    # mean and variance of difference of energies
    results = mbar.computeExpectations(E_kln, output='differences', state_dependent=True, return_dict=True)
    DeltaE_expected = results['mu']
    dDeltaE_expected = results['sigma']
    
    # mean and variance of energies squared
    results = mbar.computeExpectations(E_kln**2, state_dependent=True, return_dict=True)
    E2_expected = results['mu']
    dE2_expected = results['sigma']

    kb = physical_constants["Boltzmann constant"][0] *  physical_constants["Avogadro constant"][0] / 1000 # J (molK)^-1

    Cp_expect = (E2_expected - (E_expected*E_expected)) / (kb * t_list**2)

    return Cp_expect

def cp_bootstrap(j, temps, energies, n_estimated_points):
    print("Working on bootstrap", j,"...")
    t1 = time.time()
    energies_boot = []
    for i in range(len(energies)):
        energies_boot.append(np.random.choice(energies[i], size=len(energies[i])))
    # Solve MBAR equations
    u_kln_boot, n_samples_boot, t_list_boot, betas_boot = construct_u_kln_matrix(temps, energies_boot, add_temps = np.linspace(min(temps), max(temps), n_estimated_points))
    mbar_boot = pymbar.MBAR(u_kln_boot, n_samples_boot, verbose = False, relative_tolerance = 1e-10, initial_f_k= None, maximum_iterations=1000)
    # Compute expectations and variance of relevant terms
    E_kln_boot = u_kln_boot
    for k in range(u_kln_boot.shape[1]):
        E_kln_boot[:,k,:] *= betas_boot[k]**(-1)
    Cp_boot = calculate_cp(mbar_boot, E_kln_boot, t_list_boot)
    # compute boot strap cp
    t2 = time.time()
    print("bootstraj", j, "took", t2-t1, "seconds")
    return(Cp_boot)

def main():
    import argparse
    import plotting
    import matplotlib.pyplot as plt
    from multiprocessing import Pool
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
    energies, temps = get_energies(sim_dir_name, path, n_replicas, sim_name = sim_id)

    # Plot energy distributions of RE simualtion
    plotting.plot_RE_energy_distributions(energies)
    plt.savefig(args.figure_id + "energy_dist.png")

    # Plot energy trajectory of RE simulation
    plotting.plot_RE_energy_trajectory(energies)
    plt.savefig(args.figure_id + "energy_traj.png")

    u_kln, n_samples, t_list, betas = construct_u_kln_matrix(temps, energies, add_temps = np.linspace(min(temps), max(temps), args.n_estimated_points))

    N_boots = 100
    pool = Pool(2)

    cp_boot = pool.starmap(cp_bootstrap, ( (a, temps, energies, args.n_estimated_points) for a in range(N_boots) ))
    cp_boot = np.array(cp_boot)

    # Get bootstrapped means and variances of Cp
    mu_cp_boot = np.mean(cp_boot, axis = 0)
    sigma_cp_boot = np.std(cp_boot, axis = 0)

    plotting.plot_bootstrapped_heat_capacity(mu_cp_boot, sigma_cp_boot, t_list, args.n_replicas)
    plt.savefig(args.figure_id + "heat_capacity_bootstrap.png")

    print("Heat Capacity Analysis Complete!")
    exit()


if __name__ == "__main__":
    main()