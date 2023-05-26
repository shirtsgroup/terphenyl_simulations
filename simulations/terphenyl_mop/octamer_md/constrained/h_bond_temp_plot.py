import numpy as np
import matplotlib.pyplot as plt
from pymbar import timeseries
import glob
import sys

def main():
    dir_names = glob.glob("t_*")
    dir_names.sort()
    temperatures = []
    n_h_bonds = []
    h_bond_std = []
    for dir_name in dir_names:
        print(dir_name)
        temperatures.append(float(dir_name.split("_")[-1]))
        h_bonds = np.load(dir_name + "/h_bond_counts.npy")
        t0, g, Neff_max = timeseries.detect_equilibration(h_bonds)
        h_bonds_equil = h_bonds[t0:]
        indices = timeseries.subsample_correlated_data(h_bonds_equil, g = g)
        print("Using", len(indices), "samples out of", len(h_bonds), "H-bond counts.")
        hbonds_n = h_bonds_equil[indices].reshape(-1)
        n_h_bonds.append(np.mean(hbonds_n))
        h_bond_std.append(np.std(hbonds_n))
    print("H-bonds:", n_h_bonds)
    print("STDs:", h_bond_std)
    plt.figure(figsize=[5,5])
    plt.errorbar(temperatures, n_h_bonds, yerr=h_bond_std, capsize=3)
    plt.ylabel(r"$N_{HB}$")
    plt.xlabel("Temperature (K)")
    plt.savefig("h_bond_octamer.png", dpi = 150)
    plt.close()

if __name__ == "__main__":
    main()