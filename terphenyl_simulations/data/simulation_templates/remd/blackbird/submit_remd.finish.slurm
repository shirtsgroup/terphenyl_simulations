#!/bin/sh -l
#SBATCH -N 1 -n 316
#SBATCH --ntasks-per-node=316
#SBATCH -t 48:00:00
#SBATCH -p batch
#SBATCH -J 'mop-octamer-remd-finish'

module load gromacs/2023.4/openmpi-4.1.5-gcc-9

# Run grompp on all berendsen_nvt files
# for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f npt_new.mdp -c berendsen_npt.gro -p system_hmr.top -o npt_new -maxwarn 1; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm production_npt -ntomp 1 -nex 262144 -replex 100 -cpi production_npt.cpt
