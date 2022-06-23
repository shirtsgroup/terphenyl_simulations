#!/bin/bash
#SBATCH -N 20 --ntasks-per-node 128
#SBATCH -t 6:00:00
#SBATCH -p RM

# Run grompp on all berendsen_nvt files
for dir in sim{0..79}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f berendsen_nvt.mdp -c berendsen_npt.gro -p topol.top -o berendsen_nvt; cd ..; done

# Run all berendsen_nvt simulations in parallel

mpirun -np 2560 gmx_mpi mdrun -v -multidir sim{0..79} -deffnm berendsen_nvt -ntomp 1