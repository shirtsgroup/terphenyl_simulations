#!/bin/bash
#SBATCH -N 2 --ntasks-per-node 64
#SBATCH -t 18:00:00
#SBATCH -p amilan
#SBATCH -J 'mop-octamer-production'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


module load anaconda
conda activate plumed_mpi

module load gcc/11.2.0
module load openmpi/4.1.1
source /projects/thfo9888/software/gromacs-2022.5/bin/GMXRC


# Run grompp on all berendsen_nvt files
for dir in WALKER*; do cd $dir; mpirun -np 1 gmx_mpi convert-tpr -s npt_new.tpr -extend 100000 -o npt_new.tpr; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 128  gmx_mpi mdrun -v -deffnm npt_new -ntomp 1 -cpi npt_new.cpt -multidir WALKER_DIRS -plumed plumed_hbond_dist.dat
