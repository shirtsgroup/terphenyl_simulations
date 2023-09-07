#!/bin/bas

# Small bash script to submit 300 ns of MD simulations

# Continue production run
SLURM_OUT=$(sbatch submit.continue.slurm)
JOB_ID=$(echo "$SLURM_OUT" | awk '{print $NF}')
SLURM_OUT=$(sbatch --dependency=afterok:$JOB_ID submit.continue.slurm)
JOB_ID=$(echo "$SLURM_OUT" | awk '{print $NF}')
SLURM_OUT=$(sbatch --dependency=afterok:$JOB_ID submit.continue.slurm)
