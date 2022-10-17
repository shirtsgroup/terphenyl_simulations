#!/bin/bash
#SBATCH -N 2 --ntasks-per-node 128
#SBATCH -t 8:00:00
#SBATCH -p RM 
#SBATCH -J 'terphenyl-pom-equilibration'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

export HWLOC_HIDE_ERRORS=1

module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
source /jet/home/tfobe/software/gromacs-2022.2/bin/GMXRC

# Run grompp on all berendsen_nvt files
# for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f berendsen_npt.mdp -c berendsen_nvt.gro -p system.top -o berendsen_npt -maxwarn 1; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 256  gmx_mpi mdrun -v -multidir sim{0..63} -deffnm berendsen_npt -ntomp 1
