import os
import numpy as np
from scipy.constants import physical_constants
import panedr

kb = physical_constants["Boltzmann constant"][0] *  physical_constants["Avogadro constant"][0] / 1000 # J (molK)^-1

def get_energies(sim_dir_name, path, n_replicas, energy_type = "Enthalpy", return_T = True):
    energies = []
    temps = []
    for i in range(n_replicas):
        sim_dir = os.path.join(path, sim_dir_name + str(i))
        T = None
        with open(os.path.join(sim_dir, "npt.mdp"), "r")as f:
            for line in f.readlines():
                if "ref-t" in line:
                    T = float(line.split()[-1])
                    temps.append(T)
        sim_energies = []
        for file in os.listdir(sim_dir):
            if file.endswith(".edr"):
                edr_file = os.path.join(sim_dir, file)
                energy_df = panedr.edr_to_df(edr_file)
                pot_energy = energy_df["Enthalpy"].values
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
            u_kln[k,l,:] = betas[l] * energies[k]

    return u_kln, n_samples, t_list, betas

