#!/bin/sh -l
#SBATCH -N 1 -n 128
#SBATCH --ntasks-per-node=128
#SBATCH -t 4:00:00
#SBATCH -p batch
#SBATCH -J 'SIMULATION_NAME'
#SBATCH --hint=nomultithread

module load openmpi/4.1.5/gcc-12
source /wrk/tlf3/software/gromacs-2023.5/bin/GMXRC

# Run grompp on all berendsen_nvt files
mpirun -np 1 gmx_mpi grompp -f berendsen_nvt.mdp -c INITIAL_STRUCTURE -p TOPOLOGY_FILE -o berendsen_nvt -maxwarn 1

# Run all berendsen_nvt simulations in parallel
mpirun -np 128  gmx_mpi mdrun -v -deffnm berendsen_nvt -ntomp 1
