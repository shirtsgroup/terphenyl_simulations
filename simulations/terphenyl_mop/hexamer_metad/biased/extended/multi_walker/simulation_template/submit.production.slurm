#!/bin/bash
#SBATCH -N 1 --ntasks-per-node 32
#SBATCH -t 14:00:00
#SBATCH -p RM-shared
#SBATCH -J 'terphenyl-mop-equilibration'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load anaconda3
conda activate plumed_mpi

module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
source /jet/home/tfobe/software/gromacs-2021.5_plumed/bin/GMXRC


# Run grompp on all berendsen_nvt files
for dir in WALKER*; do cd $dir; mpirun -np 1 gmx_mpi grompp -f npt_new.mdp -c berendsen_npt.gro -p system_hmr.top -o npt_new; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 32  gmx_mpi mdrun -v -deffnm npt_new -ntomp 1 -multidir WALKER_DIRS -plumed plumed.dat
