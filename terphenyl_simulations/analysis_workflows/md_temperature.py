import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from pymbar import timeseries
import natsort
import glob
import sys
import subprocess
import os
import argparse
import terphenyl_simulations as ts

matplotlib.rc("axes", labelsize=15)

def parse_args():
    parser = argparse.ArgumentParser(
        description="A quick script to rename atoms in a pdb \
                    file. This can be used to give unique atom \
                    to all atoms in the pdb file."
    )

    parser.add_argument(
        "--shaded_region", type=int, help="limits of the shaded portion of the error plots \
                                           and violin plots input as a list with as \
                                           many pair entries as there are CVs.]",
        default=[5,7],
        nargs = "+"
    )

    parser.add_argument(
        "-p", "--prefix", type=str, help="prefix to put before file names"
    )

    parser.add_argument(
        "--gen_colvars", action='store_true', help = "flag to specifify whether or not to \
                                                      generate new COLVARS files using the \
                                                      plumed.dat file."
    )

    parser.add_argument(
        "--plumed", type = str, help = "path to the plumed.dat file used to generate COLVARS files.",
        default="../plumed_multi_cv.dat"
    )
    return parser.parse_args()


def read_plumed_data_file(filename):
    """
    Function for reading output from PRINT operations in plumed
    """
    with open(filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[1:]
    data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
    return data



def main():
    # Read cmd-line arguments
    args = parse_args()
    plt.rcParams.update({'font.size': 15})


    # Get temperatures and directory names
    # Will need to write this to be more general for other systems
    dir_names = glob.glob("t_*")
    dir_names = natsort.natsorted(dir_names)
    temperatures = [float(dir_name.split("_")[-1]) for dir_name in dir_names]
    cv_names = {"cv1" : "N Native Contacts", "cv2" : "$R_G$ (nm)", "cv3" : "$d_{end-to-end} (nm)$", "cv4" : "$d_{RMSD}$"}
    observables = {}
    ts.utils.make_path("timeseries")


    for i in range(len(dir_names)):
        os.chdir(dir_names[i])
        kt = 8.314 * 10 ** -3 * temperatures[i]
        subprocess.run(["plumed", "--no-mpi", "driver", "--plumed", args.plumed, "--kt", str(kt), "--mf_xtc", "npt_new.xtc", "--igro", "npt_new.gro"]) #, stdout=subprocess.DEVNULL)
        os.chdir("..")
        plumed_data = read_plumed_data_file(os.path.join(dir_names[i], "COLVARS"))
        for header in plumed_data.columns[2:]:
            print("Observable:", header, "Temperature:", temperatures[i])
            if header not in observables.keys():
                observables[header] = {"mean" : [], "std" : [], "data" : []}
            cv_data = plumed_data[header].values
            t0, g, Neff_max = timeseries.detect_equilibration(cv_data)
            data_equil = cv_data[t0:]
            indices = timeseries.subsample_correlated_data(data_equil, g = g)
            data_equil_ss = data_equil[indices]
            observables[header]["mean"].append(np.mean(cv_data[1000:]))
            observables[header]["std"].append(np.std(cv_data[1000:]))
            observables[header]["data"].append(cv_data[1000:])
            print("Using", len(data_equil_ss), "samples out of", len(cv_data), "H-bond counts.")

            # Plot timeseries
            plt.figure(figsize=[15,5])
            plt.plot(1/5 * plumed_data["time"].values, cv_data)
            plt.xlabel("Time (ns)")
            plt.ylabel(cv_names[header])
            plt.margins(x=0)
            plt.tight_layout()
            plt.savefig("timeseries/" + header + "_" + str(temperatures[i]) + ".png", dpi = 300)
            plt.close()




    # Plot observables
    for i, observable in enumerate(observables.keys()):
        if observable == "time":
            continue

        # Line plots
        plt.figure(figsize=[5,5])
        if len(args.shaded_region) / (2 * (i + 1) ) >= 1:
            plt.axhline(y = args.shaded_region[0], color = "r", linestyle = "--")
            plt.axhline(y = args.shaded_region[1], color = "orange", linestyle = "--")
            plt.axhspan(ymin = args.shaded_region[0], ymax = args.shaded_region[1], color = "k", alpha = 0.1)
        plt.errorbar(temperatures, observables[observable]["mean"], yerr = observables[observable]["std"])
        plt.ylabel(cv_names[observable])
        plt.xlabel("Temperature (K)")
        plt.tight_layout()

        plt.savefig(observable + "_temperature.png", dpi = 300)
        plt.close()

        # Violin plots
        fig = plt.figure(figsize=[5,5])
        if len(args.shaded_region) / (2 * (i + 1) ) >= 1:
            plt.axhline(y = args.shaded_region[0], color = "r", linestyle = "--")
            plt.axhline(y = args.shaded_region[1], color = "orange", linestyle = "--")
            plt.axhspan(ymin = args.shaded_region[0], ymax = args.shaded_region[1], color = "k", alpha = 0.1)
        plt.violinplot(observables[observable]["data"], showmedians=True)
        plt.xticks(ticks = [y + 1 for y in range(len(temperatures))], labels = [ str(int(t)) for t in  temperatures])
        plt.ylabel(cv_names[observable])
        plt.xlabel("Temperature (K)")
        plt.tight_layout()
        plt.savefig(observable + "_violin.png", dpi = 300)
        plt.close()
        

        # Write to CSV
        temperatures = np.array(temperatures)
        means = np.array(observables[observable]["mean"])
        stds = np.array(observables[observable]["std"])
        table = np.vstack((temperatures, means, stds)).T
        np.savetxt(observable + "_temp.csv", table, header="temperature,mean,std")
    

    





if __name__ == "__main__":
    main()