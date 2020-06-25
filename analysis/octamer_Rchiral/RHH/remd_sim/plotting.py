import matplotlib.pyplot as plt
import numpy as np


def plot_RE_energy_distributions(energies, energy_label = "Enthalpy ($kJ mol^{-1}$)"):
    plt.figure(figsize=[8,4], dpi=600)
    for i in range(len(energies)):
        # requires cutting off outlier (likely initial frame)
        plt.hist(energies[i][energies[i] > np.min(energies[i])], bins = 50, alpha = 0.6)
        plt.xlabel(energy_label)
        plt.ylabel("Counts")
    # plt.legend([str(np.round(t, 1)) for t in temps])

def plot_RE_energy_trajectory(energies, energy_label = "Enthalpy ($kJ mol^{-1}$)"):
    plt.figure(figsize=[8,4], dpi=600)
    for i in range(len(energies)):
        plt.plot(energies[i][energies[i] > np.min(energies[i])], lw=0.5)
        plt.ylabel(energy_label)
        plt.xlabel("Time Steps")


def plot_bootstrapped_heat_capacity(mu_cp_boot, sigma_cp_boot, t_list, n_replicas):
    plt.figure(figsize=[8,4], dpi=600)
    plt.plot(t_list[n_replicas:], mu_cp_boot[n_replicas:], color = "gray")
    plt.fill_between(t_list[n_replicas:], mu_cp_boot[n_replicas:]-2*sigma_cp_boot[n_replicas:],  mu_cp_boot[n_replicas:]+2*sigma_cp_boot[n_replicas:], color="gray", alpha = 0.3)
    plt.scatter(t_list[:n_replicas], mu_cp_boot[:n_replicas], color = "red", s = 5)
    plt.errorbar(t_list[:n_replicas], mu_cp_boot[:n_replicas], yerr=2*sigma_cp_boot[:n_replicas], fmt=".")
    plt.ylabel("Heat Capacity ($C_p$)")
    plt.xlabel("Temperature (K)")