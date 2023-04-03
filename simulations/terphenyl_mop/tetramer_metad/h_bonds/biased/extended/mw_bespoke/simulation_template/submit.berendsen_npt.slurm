#!/bin/bash
#SBATCH -N 1 --ntasks-per-node 32
#SBATCH -t 3:00:00
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
for dir in WALKER*; do cd $dir; mpirun -np 1 gmx_mpi grompp -f berendsen_npt.mdp -c berendsen_nvt.gro -p system_bespoke_hmr.top -o berendsen_npt -maxwarn 2; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 32  gmx_mpi mdrun -v -deffnm berendsen_npt -ntomp 1 -multidir WALKER_DIRS
