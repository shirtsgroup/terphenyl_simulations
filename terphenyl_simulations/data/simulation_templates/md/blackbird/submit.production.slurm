#!/bin/sh -l
#SBATCH -N 1 -n 128
#SBATCH --ntasks-per-node=128
#SBATCH -t 24:00:00
#SBATCH -p batch
#SBATCH -J 'SIMULATION_NAME'
#SBATCH --hint=nomultithread

module load openmpi/4.1.5/gcc-12
source /wrk/tlf3/software/gromacs-2023.5/bin/GMXRC

# Berendsen NVT Equilibration
mpirun -np 1 gmx_mpi grompp -f berendsen_nvt.mdp -c INITIAL_STRUCTURE -p TOPOLOGY_FILE -o berendsen_nvt -maxwarn 1

# Run simulation
mpirun -np 128  gmx_mpi mdrun -v -deffnm berendsen_nvt -ntomp 1

# Berendsen NPT Equilibration
mpirun -np 1 gmx_mpi grompp -f berendsen_npt.mdp -c berendsen_nvt.gro -p TOPOLOGY_FILE -o berendsen_npt -maxwarn 2

# Run simulation
mpirun -np 128  gmx_mpi mdrun -v -deffnm berendsen_npt -ntomp 1


# Production NPT Run
mpirun -np 1 gmx_mpi grompp -f production_npt.mdp -c berendsen_npt.gro -p TOPOLOGY_FILE -o production_npt -maxwarn 1

# Run simulation
mpirun -np 128 gmx_mpi mdrun -v -deffnm production_npt -ntomp 1