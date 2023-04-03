import signac
import os
import glob
import sys
import flow
import subprocess
import matplotlib.pyplot as plt
import plumed
from flow import FlowProject

# Scripts for running the submit_all.slurm script which submits all simulations
# at once with dependencies linking them. If this fails individual submission
# flow operations can be found below.

@FlowProject.label
def check_submitted(job):

    if job.isfile("status.txt"):
        with open(job.fn("status.txt"), "r") as f:
            lines = f.readlines()
        
        if "SUBMITTED" in lines:
            return True
    return False


@FlowProject.post(check_submitted)
@FlowProject.operation
def submit_all_simulations(job):
    os.chdir(job.path)
    n_jobs_old = len(subprocess.check_output(["squeue", "-u", "tfobe"]).splitlines()) - 1
    subprocess.run(["bash", "submit_all.slurm"])
    n_jobs = len(subprocess.check_output(["squeue", "-u", "tfobe"]).splitlines()) - 1
    jobs_submitted = n_jobs - n_jobs_old
    if jobs_submitted == 3:
        with open(job.fn("status.txt"), "w") as f:
            f.write("SUBMITTED")
    else:
        with open(job.fn("status.txt"), "w") as f:
            f.write("FAILED")

# Labels for individual simulation components

@FlowProject.label
def check_berendsen_nvt_start(job):
    if job.isfile("berendsen_nvt.log"):
        return True
    return False

@FlowProject.label
def check_berendsen_nvt_finish(job):
    if job.isfile("berendsen_nvt.gro"):
        return True
    return False

@FlowProject.label
def check_berendsen_npt_start(job):
    if job.isfile("berendsen_npt.log"):
        return True
    return False

@FlowProject.label
def check_berendsen_npt_finish(job):
    if job.isfile("berendsen_npt.gro"):
        return True
    return False

@FlowProject.label
def check_production_npt_start(job):
    if job.isfile("npt_new.log"):
        return True
    return False

@FlowProject.label
def check_production_npt_finish(job):
    if job.isfile("npt_new.gro"):
        return True
    return False

# Flow operations to run specific simulation parts

@FlowProject.post(check_berendsen_nvt_finish)
@FlowProject.operation
def submit_berendsen_nvt_simulations(job):
    os.chdir(job.path)
    subprocess.run(["sbatch", "submit.berendsen_nvt.slurm"])

@FlowProject.pre(check_berendsen_nvt_finish)
@FlowProject.post(check_berendsen_npt_finish)
@FlowProject.operation
def submit_berendsen_npt_simulations(job):
    os.chdir(job.path)
    subprocess.run(["sbatch", "submit.berendsen_npt.slurm"])

@FlowProject.pre(check_berendsen_npt_finish)
@FlowProject.post(check_production_npt_finish)
@FlowProject.operation
def submit_production_npt_simulations(job):
    os.chdir(job.path)
    slurm_id = subprocess.check_output(["sbatch", "submit.production.slurm"])

@FlowProject.pre(check_production_npt_finish)
@FlowProject.operation
def continue_production_npt_simulation(job):
    os.chdir(job.path)
    slurm_id = subprocess.check_output(["sbatch", "submit.continue.slurm"])

# Utility operation to remove all backup files
@FlowProject.label
def has_backup_files(job):
    backup_files = glob.glob(job.fn("#*#"))
    backup_files += glob.glob(job.fn("bck.*"))
    if len(backup_files) == 0:
        return False
    return True

@FlowProject.label
def has_failed_step_files(job):
    step_files = glob.glob(job.fn("step*.pdb"))
    if len(step_files) == 0:
        return False
    return True

@FlowProject.pre(has_backup_files)
@FlowProject.post.not_(has_backup_files)
@FlowProject.operation
def remove_backup_files(job):
    os.chdir(job.path)
    backup_files = glob.glob("#*#")
    backup_files += glob.glob("bck.*")
    for file in backup_files:
        print("Removing", os.path.abspath(file))
        os.remove(file)

@FlowProject.pre(has_failed_step_files)
@FlowProject.post.not_(has_failed_step_files)
@FlowProject.operation
def remove_failed_step_files(job):
    os.chdir(job.path)
    step_files = glob.glob("step*.pdb")
    for file in step_files:
        print("Removing", os.path.abspath(file))
        os.remove(file)

@FlowProject.pre(check_production_npt_start)
@FlowProject.post.isfile("CV_bias_plot.png")
@FlowProject.operation
def plot_CV_bias(job):
    plt.figure(dpi=300)
    stride = 100
    data = plumed.read_as_pandas(job.fn("HBOND_SUMS"))
    data_2 = plumed.read_as_pandas(job.fn("HILLS"))
    fig, ax = plt.subplots(6,1, figsize = [22.5,15])
    ax[0].set_title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    ax[0].plot(1/1000 * data.time[::stride], data.values[:,1][::stride])
    ax[0].set_ylabel("Hydrogen Bond 1")
    ax[1].plot(1/1000 * data.time[::stride], data.values[:,2][::stride])
    ax[1].set_ylabel("Hydrogen Bond 2")
    ax[2].plot(1/1000 * data.time[::stride], data.values[:,3][::stride])
    ax[2].set_ylabel("Hydrogen Bond 3")
    ax[3].plot(1/1000 * data.time[::stride], data.values[:,4][::stride])
    ax[3].set_ylabel("Sum of H-Bonds")
    ax[4].plot(1/1000 * data.time[::stride], data.values[:,5][::stride],'r')
    ax[4].set_ylabel("Bias Energy")
    ax[5].plot(1/1000 * data_2.time, data_2.values[:,3],'b')
    ax[5].set_ylabel("Gaussian Height")
    ax[5].set_xlabel("Time (ns)")
    plt.savefig(job.fn("CV_bias_plot.png"), dpi = 300)
    plt.close()

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("CV_bias_plot_zoomed.png")
@FlowProject.operation
def plot_CV_bias_zoomed(job):
    plt.figure(dpi=300)
    stride = 100
    data = plumed.read_as_pandas(job.fn("HBOND_SUMS"))
    data_2 = plumed.read_as_pandas(job.fn("HILLS"))
    fig, ax = plt.subplots(6,1, figsize = [22.5,15])
    print(len(data.values[:,1]))
    print(len(data_2.values[:,3]))
    ax[0].set_title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    ax[0].plot(1/1000 * data.time[2000000:2666667], data.values[:,1][2000000:2666667], linewidth=0.5)
    ax[0].set_ylabel("Hydrogen Bond 1")
    ax[1].plot(1/1000 * data.time[2000000:2666667], data.values[:,2][2000000:2666667], linewidth=0.5)
    ax[1].set_ylabel("Hydrogen Bond 2")
    ax[2].plot(1/1000 * data.time[2000000:2666667], data.values[:,3][2000000:2666667], linewidth=0.5)
    ax[2].set_ylabel("Hydrogen Bond 3")
    ax[3].plot(1/1000 * data.time[2000000:2666667], data.values[:,4][2000000:2666667], linewidth=0.5)
    ax[3].set_ylabel("Sum of H-Bonds")
    ax[4].plot(1/1000 * data.time[2000000:2666667], data.values[:,5][2000000:2666667],'r', linewidth=0.5)
    ax[4].set_ylabel("Bias Energy")
    ax[5].plot(1/1000 * data_2.time[40000:53333], data_2.values[:,3][40000:53333],'b', linewidth=0.5)
    ax[5].set_ylabel("Gaussian Height")
    ax[5].set_xlabel("Time (ns)")
    plt.savefig(job.fn("CV_bias_plot_zoomed.png"), dpi = 300)
    plt.close()

@FlowProject.pre(check_production_npt_finish)
@FlowProject.post.isfile("sum_hills_FE.png")
@FlowProject.operation
def calculate_sum_hills_FE(job):
    current_dir = os.path.abspath("")
    kt = 300 * 8.314462618 * 10 ** -3
    os.chdir(job.fn(""))
    subprocess.run(["plumed", "sum_hills", "--hills", "HILLS", "--kt", str(kt)]) # Silences output from sum_hills
    fes_data = plumed.read_as_pandas("fes.dat")
    print(fes_data)
    filtered = fes_data[fes_data["n_hbonds"] >= -0.1]
    filtered = filtered[filtered["n_hbonds"] <= 3.2]
    plt.figure()
    plt.xlim([0,3])
    plt.plot(filtered.values[:,0], filtered.values[:,1])
    plt.title("SIGMA: " +  str(job.sp.sigma) +  " HEIGHT:" + str(job.sp.height) + " BF:" + str(job.sp.bf))
    plt.xlabel("Number of Hydrogen Bonds")
    plt.ylabel("Free Energy (kJ/mol)")
    plt.grid(visible=True, which="both", axis="both")
    plt.savefig("sum_hills_FE.png")
    plt.close()
    os.chdir(current_dir)
    

if __name__ == '__main__':
    FlowProject().main()