import signac
import os
import glob
import flow
import subprocess
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
    if jobs_submitted == 5:
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

# Utility operation to remove all backup files
@FlowProject.label
def has_backup_files(job):
    backup_files = glob.glob(job.fn("#*#"))
    backup_files += glob.glob(job.fn("bck.*"))
    if len(backup_files) == 0:
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


if __name__ == '__main__':
    FlowProject().main()