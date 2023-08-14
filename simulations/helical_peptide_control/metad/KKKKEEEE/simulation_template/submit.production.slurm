#!/bin/sh -l
#SBATCH -N 1 --ntasks-per-node 16
#SBATCH -t 18:00:00
#SBATCH -p amilan
#SBATCH -J 'peptide-production'
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
for dir in WALKER*; do cd $dir; mpirun -np 1 gmx_mpi grompp -f npt_new.mdp -c berendsen_npt.gro -p system_hmr.top -o npt_new; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 16  gmx_mpi mdrun -v -deffnm npt_new -ntomp 1 -multidir WALKER_DIRS -plumed plumed_hbond_dist.dat
