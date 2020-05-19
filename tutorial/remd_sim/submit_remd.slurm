#!/bin/bash

#SBATCH -p sgpu
#SBATCH -N 5
#SBATCH -t 24:00:00  # max wall time is 24 hrs

module load intel
module load impi
module load mkl

source /projects/beco4952/pkgs/gromacs/2018_gpu/bin/GMXRC

sbatch submit_remd.slurm -d=afterok:$SLURM_JOB_ID

if [ -f sim0/npt.cpt ]; then
    mpirun -np 20 gmx_mpi mdrun -v -multidir sim{0..19} -deffnm npt -cpi npt.cpt -append -ntomp 6 -nex 8000 -replex 100
else
    for dir in sim{0..19};
        do cd $dir;
	echo $dir;
        mpirun -np 1 gmx_mpi grompp -f npt.mdp -p topol.top -c berendsen.gro -o npt
        cd ..;
    done 
    mpirun -np 20 gmx_mpi mdrun -v -multidir sim{0..19} -deffnm npt -ntomp 6 -nex 8000 -replex 100
fi
