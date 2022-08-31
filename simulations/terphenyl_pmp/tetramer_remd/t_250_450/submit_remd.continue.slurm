#!/bin/bash
#SBATCH -N 2 --ntasks-per-node 128
#SBATCH -t 36:00:00
#SBATCH -p RM 
#SBATCH -J 'pmp-tetramer-production'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
source /jet/home/tfobe/software/gromacs-2022.2/bin/GMXRC


# Extend npt_new.tpr file using convert-tpr
for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi convert-tpr -s npt_new.tpr -extend 100000 -o npt_new.tpr; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 256 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm npt_new -ntomp 1 -nex 262144 -replex 100 -cpi npt_new.cpt
