#!/bin/bash

# NVT Equilibration
SLURM_OUT=$(sbatch submit.berendsen_nvt.slurm)
JOB_ID=$(echo "$SLURM_OUT" | awk '{print $NF}')

# NPT Equilibration
SLURM_OUT=$(sbatch --dependency=afterok:$JOB_ID submit.berendsen_npt.slurm)
JOB_ID=$(echo "$SLURM_OUT" | awk '{print $NF}')

# Production Run
SLURM_OUT=$(sbatch --dependency=afterok:$JOB_ID submit.production.slurm)
JOB_ID=$(echo "$SLURM_OUT" | awk '{print $NF}')
SLURM_OUT=$(sbatch --dependency=afterok:$JOB_ID submit.continue.slurm)
JOB_ID=$(echo "$SLURM_OUT" | awk '{print $NF}')
SLURM_OUT=$(sbatch --dependency=afterok:$JOB_ID submit.continue.slurm)

