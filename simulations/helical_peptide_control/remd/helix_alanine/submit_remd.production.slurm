#!/bin/sh -l
#SBATCH -N 1 --ntasks-per-node 64
#SBATCH -t 24:00:00
#SBATCH -p amilan
#SBATCH -J 'mop-octamer-equilibration'
#SBATCH --account=ucb368_asc1
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --constraint=ib 

module purge
module load anaconda
conda activate plumed_mpi

module load gcc/11.2.0
module load openmpi/4.1.1
export SLURM_EXPORT_ENV=ALL

source /projects/thfo9888/software/gromacs-2022.5/bin/GMXRC


# Run grompp on all berendsen_nvt files
# for dir in sim{0..31}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f npt_new.mdp -c berendsen_npt.gro -p system_hmr.top -o npt_new -maxwarn 1; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 64 gmx_mpi mdrun -v -multidir sim{0..31} -deffnm npt_new -ntomp 1 -nex 32768 -replex 100 
