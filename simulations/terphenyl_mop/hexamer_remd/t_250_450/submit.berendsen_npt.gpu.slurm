#!/bin/bash
#SBATCH -N 2
#SBATCH -t 2:00:00
#SBATCH -p GPU
#SBATCH --gres=gpu:v100-32:8
#SBATCH -J 'terphenyl-pom-equilibration'
#SBATCH --mail-user=theodore.fobe@colorado.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0
module load cuda/11.1.1
source /jet/home/tfobe/software/gromacs-2022.2_gpu/bin/GMXRC


# Run grompp on all berendsen_nvt files
for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f berendsen_npt.mdp -c berendsen_nvt.gro -p system_hmr.top -o berendsen_npt -maxwarn 2; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 64  gmx_mpi mdrun -v -multidir sim{0..63} -deffnm berendsen_npt -ntomp 1
