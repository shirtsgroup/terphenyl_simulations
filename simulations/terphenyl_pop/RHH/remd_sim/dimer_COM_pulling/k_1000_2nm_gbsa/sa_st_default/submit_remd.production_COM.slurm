#!/bin/bash
#SBATCH -N 1 --ntasks-per-node 128
#SBATCH -t 24:00:00
#SBATCH -p RM 
#SBATCH -J 'dimer-implicit-solvent-production'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=theodore.fobe@colorado.edu

export HWLOC_HIDE_ERRORS=1

module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
source /jet/home/tfobe/software/gromacs-2018.2_mpi/bin/GMXRC
# /jet/home/susa/pkgs/gromacs/2020.5_cpu/bin/GMXRC

# Run grompp on all berendsen_nvt files
# for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f nvt_COM_pulling_gbsa.mdp -o nvt_COM_pulling_gbsa -p topol.top -c octamer_dimer_box.gro -n index.ndx; cd ..; done


# Run all berendsen_nvt simulations in parallel
mpirun -np 128 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm nvt_COM_pulling_gbsa -nex 8000 -replex 100
