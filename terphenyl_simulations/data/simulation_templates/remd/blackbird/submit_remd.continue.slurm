#!/bin/sh -l
#SBATCH -N 1 -n 316
#SBATCH --ntasks-per-node=316
#SBATCH -t 48:00:00
#SBATCH -p batch
#SBATCH -J 'mop-octamer-remd-continue'

module load gromacs/2023.4/openmpi-4.1.5-gcc-9

# Extend npt_new.tpr file using convert-tpr
for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi convert-tpr -s production_npt.tpr -extend 100000 -o production_npt.tpr; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm production_npt -ntomp 1 -nex 262144 -replex 100 -cpi production_npt.cpt
