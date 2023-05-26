import signac
import os
import sys
import glob
import flow
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import glob
import mdtraj as md
from natsort import natsorted
from flow import FlowProject
import terphenyl_simulations as hs

# Scripts for running the submit_all.slurm script which submits all simulations
# at once with dependencies linking them. If this fails individual submission
# flow operations can be found below.




@FlowProject.label
def sw_check_nvt_equilibration(job):
    return(os.path.exists(job.fn("berendsen_nvt.gro")))

@FlowProject.label
def sw_check_npt_equilibration(job):
    return(os.path.exists(job.fn("berendsen_npt.gro")))

@FlowProject.label
def sw_check_npt_production(job):
    return(os.path.exists(job.fn("npt_new.gro")))

@FlowProject.post(sw_check_nvt_equilibration)
@FlowProject.post(sw_check_npt_equilibration)
@FlowProject.post(sw_check_npt_production)
@FlowProject.operation
def run_all_simulations(job):
    os.chdir(job.fn(""))
    print("Working on", job.fn(""))
    process = subprocess.Popen("bash run_local.sh".split(" "))
    process.wait()

if __name__ == '__main__':
    FlowProject().main()
