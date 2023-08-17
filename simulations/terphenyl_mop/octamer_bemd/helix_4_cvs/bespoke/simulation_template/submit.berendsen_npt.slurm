#!/bin/sh -l
#SBATCH -N 1 --ntasks-per-node 64
#SBATCH -t 4:00:00
#SBATCH -p amilan
#SBATCH -J 'mop-octamer-equilibration'
#SBATCH --account=ucb368_asc1
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --constraint=ib 


module load anaconda
conda activate plumed_mpi

module load gcc/11.2.0
module load openmpi/4.1.1
SLURM_EXPORT_ENV=ALL

source /projects/thfo9888/software/gromacs-2022.5/bin/GMXRC

# Run grompp on all berendsen_nvt files
for dir in WALKER*; do cd $dir; mpirun -np 1 gmx_mpi grompp -f berendsen_npt.mdp -c berendsen_nvt.gro -p system_hmr_rest.top -o berendsen_npt -maxwarn 2 -r em_solvated_helix.gro; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 32  gmx_mpi mdrun -v -deffnm berendsen_npt -ntomp 1 -multidir WALKER_DIRS