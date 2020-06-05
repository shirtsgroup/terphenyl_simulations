#!/bin/bash

#SBATCH -p sgpu
#SBATCH -N 1 
#SBATCH -t 24:00:00  # max wall time is 24 hrs
#SBATCH -J terphenyl-230K

module load intel
module load impi
module load mkl

source /projects/beco4952/pkgs/gromacs/2018_gpu/bin/GMXRC

# mpirun -np 1 gmx_mpi grompp -f npt_berendsen.mdp -p topol.top -c em_solvated.gro -o berendsen
# mpirun -np 4 gmx_mpi mdrun -v -deffnm berendsen -ntomp 6

mpirun -np 1 gmx_mpi grompp -f npt.mdp -p topol.top -c berendsen.gro -o npt 
mpirun -np 4 gmx_mpi mdrun -v -deffnm npt -ntomp 6
