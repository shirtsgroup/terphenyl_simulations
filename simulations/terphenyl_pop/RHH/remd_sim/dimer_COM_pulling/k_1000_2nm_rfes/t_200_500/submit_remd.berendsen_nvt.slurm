#!/bin/bash
#SBATCH -N 7 --ntasks-per-node 128
#SBATCH -t 08:00:00
#SBATCH -p RM 
#SBATCH -J 'dimer-nvt-equilibration-G1'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
source /jet/home/susa/pkgs/gromacs/2020.5_cpu/bin/GMXRC


# Run grompp on all berendsen_nvt files
# for dir in sim{0..55}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f berendsen_nvt.mdp -c em.gro -p topol.top -o berendsen_nvt -n index.ndx; cd ..; done


# Run all berendsen_nvt simulations in parallel
mpirun -np 896  gmx_mpi mdrun -v -multidir sim{0..55} -deffnm berendsen_nvt -ntomp 1
