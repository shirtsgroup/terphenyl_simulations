import yaml
import os
import pandas as pd
from tqdm import tqdm
import glob
import subprocess
import sys


def read_yaml_parameter_file(yaml_file_name):
    with open(yaml_file_name, 'r') as file:
        yaml_data = yaml.safe_load(file)

    return yaml_data

def check_walker_file(job, filename, walker_dirs = ["WALKER0", "WALKER1", "WALKER2", "WALKER3"]):
    walkers = []
    for walker_dir in walker_dirs:
        walkers.append(job.isfile(os.path.join(walker_dir, filename)))
    return all(walkers)

def reweight_walker_trajectories(job, plumed_file, kt, gro_file = "npt_new.gro", xtc_file = "npt_new.xtc"):
    # Get original directory path
    current_dir = os.path.abspath("")

    # Change directory to job
    os.chdir(job.fn(""))
    
    print("Reweighting simulations...")
    for walker_dir in tqdm(glob.glob("WALKER*")):
        # navigate to specific WALKER dir
        os.chdir(walker_dir)

        # Run plumed driver to reweight biases of invidual simulations
        # We also use this to get H-bond measures from each frame from npt_new.xtc
        subprocess.run(["plumed", "--no-mpi", "driver", "--plumed", plumed_file, "--kt", str(kt), "--mf_xtc", xtc_file, "--igro", gro_file]) # Silences output from sum_hills
        os.chdir(job.fn(""))
    
    # Return to original directory
    os.chdir(current_dir)

def reweight_walker_colvars(job, plumed_file, kt):
    # Get original directory path
    current_dir = os.path.abspath("")

    # Change directory to job
    os.chdir(job.fn(""))
    
    print("Reweighting simulations...")
    for walker_dir in tqdm(glob.glob("WALKER*")):
        # navigate to specific WALKER dir
        os.chdir(walker_dir)

        # Run plumed driver to reweight biases of invidual simulations
        # We also use this to get H-bond measures from each frame from npt_new.xtc
        subprocess.run(["plumed", "--no-mpi", "driver", "--plumed", plumed_file, "--kt", str(kt), "--noatoms"]) # Silences output from sum_hills
        os.chdir(job.fn(""))
    
    # Return to original directory
    os.chdir(current_dir)

def read_plumed_data_file(filename):
    """
    Function for reading output from PRINT operations in plumed
    """
    with open(filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[1:]
    data = pd.read_csv(filename, sep = " ", header = None, comment="#", names = headers)
    return data

def read_plumed_hills_file(hills_file):
    """
    Function for reading HILLS output file from plumed driver
    """
    with open(hills_file) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    hills_data = pd.read_csv(hills_file, skiprows = 3, delim_whitespace=True, header = None, comment="#", names = headers)
    return hills_data

def read_plumed_fes_file(fes_filename):
    """
    Function for reading FES output from plumed sum_hills function
    """
    with open(fes_filename) as f:
        headers = f.readline().strip()
    headers = headers.split(" ")[2:]
    fes_data = pd.read_csv(fes_filename, skiprows = 5, delim_whitespace=True, header = None, comment="#", names = headers)
    return fes_data