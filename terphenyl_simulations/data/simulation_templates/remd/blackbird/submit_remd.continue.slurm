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
for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi convert-tpr -s production_npt.tpr -extend 1000000 -o production_npt.tpr; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm production_npt -ntomp 1 -nex 262144 -replex 100 -cpi production_npt.cpt
