import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymbar import timeseries
import natsort
import glob
import sys
import subprocess
import os
import terphenyl_simulations as ts

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
    # Get temperatures and directory names
    dir_names = glob.glob("t_*")
    dir_names = natsort.natsorted(dir_names)
    temperatures = [float(dir_name.split("_")[-1]) for dir_name in dir_names]
    cv_names = {"cv1" : "N Native Contacts", "cv2" : "$R_G$ (nm)", "cv3" : "$d_{end-to-end} (nm)$", "cv4" : "$d_{RMSD}$"}
    observables = {}
    ts.utils.make_path("timeseries")


    for i in range(len(dir_names)):
        os.chdir(dir_names[i])
        kt = 8.314 * 10 ** -3 * temperatures[i]
        subprocess.run(["plumed", "--no-mpi", "driver", "--plumed", "../plumed_multi_cv.dat", "--kt", str(kt), "--mf_xtc", "npt_new.xtc", "--igro", "npt_new.gro"], stdout=subprocess.DEVNULL)
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
            observables[header]["mean"].append(np.mean(cv_data))
            observables[header]["std"].append(np.std(cv_data))
            observables[header]["data"].append(cv_data[500:])
            print("Using", len(data_equil_ss), "samples out of", len(cv_data), "H-bond counts.")

            # Plot timeseries
            plt.figure(figsize=[5,5])
            plt.plot(1/5 * plumed_data["time"].values, cv_data)
            plt.xlabel("Time (ns)")
            plt.ylabel(cv_names[header])
            plt.savefig("timeseries/" + header + "_" + str(temperatures[i]) + ".png")
            plt.close()




    # Plot observables
    for observable in observables.keys():
        if observable == "time":
            continue

        # Line plots
        plt.figure(figsize=[5,5])
        plt.errorbar(temperatures, observables[observable]["mean"], yerr = observables[observable]["std"])
        plt.ylabel(cv_names[observable])
        plt.xlabel("Temperature (K)")
        plt.savefig(observable + "_temperature.png")
        plt.close()

        # Violin plots
        fig = plt.figure(figsize=[5,5])
        plt.violinplot(observables[observable]["data"], showmedians=True)
        plt.xticks(ticks = [y + 1 for y in range(len(temperatures))], labels = [ str(int(t)) for t in  temperatures])
        plt.ylabel(cv_names[observable])
        plt.xlabel("Temperature (K)")
        plt.savefig(observable + "_violin.png")
        plt.close()
        

        # Write to CSV
        temperatures = np.array(temperatures)
        means = np.array(observables[observable]["mean"])
        stds = np.array(observables[observable]["std"])
        table = np.vstack((temperatures, means, stds)).T
        np.savetxt(observable + "_temp.csv", table, header="temperature,mean,std")
    

    





if __name__ == "__main__":
    main()