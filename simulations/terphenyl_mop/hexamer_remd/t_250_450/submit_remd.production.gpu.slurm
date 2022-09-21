#!/bin/bash
#SBATCH -N 2
#SBATCH -t 20:00:00
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
for dir in sim{0..63}; do cd $dir; mpirun -np 1 gmx_mpi grompp -f npt_new.mdp -c berendsen_npt.gro -p system_hmr.top -o npt_new -maxwarn 1; cd ..; done

# Run all berendsen_nvt simulations in parallel
mpirun -np 64 gmx_mpi mdrun -v -multidir sim{0..63} -deffnm npt_new -ntomp 1 -nex 262144 -replex 100
