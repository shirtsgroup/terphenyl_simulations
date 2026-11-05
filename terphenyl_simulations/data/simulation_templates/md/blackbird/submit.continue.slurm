#!/bin/sh -l
#SBATCH -N 1 -n 128
#SBATCH --ntasks-per-node=128
#SBATCH -t 168:00:00
#SBATCH -p batch
#SBATCH -J 'SIMULATION_NAME'
#SBATCH --hint=nomultithread

module load openmpi/4.1.5/gcc-12
source /wrk/tlf3/software/gromacs-2023.5/bin/GMXRC

# Extend npt_new.tpr file using convert-tpr
mpirun -np 1 gmx_mpi convert-tpr -s production_npt.tpr -extend 1000000 -o production_npt.tpr

# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -deffnm production_npt -ntomp 1 -cpi production_npt.cpt
