#!/bin/sh -l
#SBATCH -N 1 -n 128
#SBATCH --ntasks-per-node=128
#SBATCH -t 24:00:00
#SBATCH -p batch
#SBATCH -J 'SIMULATION_NAME'
#SBATCH --hint=nomultithread

module load openmpi/4.1.5/gcc-12
source /wrk/tlf3/software/gromacs-2023.5/bin/GMXRC

# Run grompp on all berendsen_nvt files
mpirun -np 1 gmx_mpi grompp -f production_npt.mdp -c berendsen_npt.gro -p TOPOLOGY_FILE -o production_npt -maxwarn 1

# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -deffnm production_npt -ntomp 1