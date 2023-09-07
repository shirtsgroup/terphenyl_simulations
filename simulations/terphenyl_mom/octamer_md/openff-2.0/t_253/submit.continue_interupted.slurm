#!/bin/sh -l
#SBATCH -N 1 --ntasks-per-node 8
#SBATCH -t 4:00:00
#SBATCH -p amilan
#SBATCH -J 'mop-octamer-production'
#SBATCH --account=ucb368_asc1
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module purge
module load anaconda
conda activate plumed_mpi

module load gcc/11.2.0
module load openmpi/4.1.1
source /projects/thfo9888/software/gromacs-2022.5/bin/GMXRC

# Run all berendsen_nvt simulations in parallel
mpirun -np 8  gmx_mpi mdrun -v -deffnm npt_new -ntomp 1 -cpi npt_new.cpt
