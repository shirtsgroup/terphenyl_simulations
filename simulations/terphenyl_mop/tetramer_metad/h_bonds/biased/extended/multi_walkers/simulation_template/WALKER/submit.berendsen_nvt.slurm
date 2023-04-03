#!/bin/bash
#SBATCH -N 1 --ntasks-per-node 8
#SBATCH -t 1:00:00
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
mpirun -np 1 gmx_mpi grompp -f berendsen_nvt.mdp -c terphenyl_extended.gro -p system_hmr.top -o berendsen_nvt -maxwarn 1

# Run all berendsen_nvt simulations in parallel
mpirun -np 8  gmx_mpi mdrun -v -deffnm berendsen_nvt -ntomp 1
