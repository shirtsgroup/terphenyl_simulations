#!/bin/sh -l
#SBATCH -N 1 -n 316
#SBATCH --ntasks-per-node=316
#SBATCH -t 24:00:00
#SBATCH -p batch
#SBATCH -J 'mop-octamer-remd'

module load gromacs/2023.4/openmpi-4.1.5-gcc-9

# Run grompp on all berendsen_nvt files
for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f production_npt.mdp -c berendsen_npt.gro -p TOPOLOGY_FILE -o production_npt -maxwarn 1; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm production_npt -ntomp 1 -nex 262144 -replex 100
# mpirun -np 256 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm npt_new -ntomp 1 -nex 262144 -replex 100 -cpi npt_new.cpt
