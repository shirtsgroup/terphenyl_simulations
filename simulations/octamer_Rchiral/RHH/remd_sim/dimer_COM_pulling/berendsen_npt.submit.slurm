#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 32
#SBATCH -J berendsen_nvt 
#SBATCH -t 6:00:00 


module load gromacs/2020.2-cpu

# for dir in sim{0..19}; do cd $dir; gmx grompp -f berendsen.mdp -c em_solvated.gro -p topol.top -o berendsen; cd ..; done

mpirun -np 32 gmx_mpi mdrun -v -deffnm berendsen_npt -ntomp 1
