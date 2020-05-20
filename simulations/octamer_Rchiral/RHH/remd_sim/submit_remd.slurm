#!/bin/bash

#SBATCH -p sgpu
#SBATCH -N 5
#SBATCH -t 24:00:00  # max wall time is 24 hrs
#SBATCH --ntasks-per-node 4

module load intel
module load impi
module load mkl

source /projects/beco4952/pkgs/gromacs/2018_gpu/bin/GMXRC

sbatch submit_remd.slurm -d=afterok:$SLURM_JOB_ID

# IF statement to see if there is a checkpoint file in the sim0 dir
if [ -f sim0/npt.cpt ]; then
    # if checkpoints are present, we continue the simulation
    for dir in sim{0..19};
        do cd $dir;
        echo $dir;
	# extend simulation by 15 ns (time is determined empirically)
        mpirun -np 1 gmx_mpi convert-tpr -s npt.tpr -extend 15000.000 -o npt.tpr 
        cd ..;
    done
    # Run simulation from top directory
    mpirun -np 20 gmx_mpi mdrun -v -multidir sim{0..19} -deffnm npt -cpi npt.cpt -noappend -ntomp 6 -nex 8000 -replex 100
else
    # For loop for initializing simulation
    for dir in sim{0..19};
        do cd $dir;
	echo $dir;
        mpirun -np 1 gmx_mpi grompp -f npt.mdp -p topol.top -c npt_restart.gro -o npt 
        cd ..;
    done
    # Initial run of simulation
    mpirun -np 20 gmx_mpi mdrun -v -multidir sim{0..19} -deffnm npt -ntomp 6 -nex 8000 -replex 100
fi
