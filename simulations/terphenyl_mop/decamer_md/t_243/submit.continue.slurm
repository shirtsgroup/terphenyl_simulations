#!/bin/bash
#SBATCH -N 1 --ntasks-per-node 8
#SBATCH -t 12:00:00
#SBATCH -p RM-shared
#SBATCH -J 'terphenyl-mop-equilibration'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
source /jet/home/tfobe/software/gromacs-2022.2/bin/GMXRC


# Run grompp on all berendsen_nvt files
mpirun -np 1 gmx_mpi convert-tpr -s npt_new.tpr -extend 100000 -o npt_new.tpr

# Run all berendsen_nvt simulations in parallel
mpirun -np 8  gmx_mpi mdrun -v -deffnm npt_new -ntomp 1 -cpi npt_new.cpt
