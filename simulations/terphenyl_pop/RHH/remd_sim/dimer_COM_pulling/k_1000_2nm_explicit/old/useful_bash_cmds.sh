# apply grompp to berendsen.mdp files for all replicas
for dir in sim{0..63}; do cd $dir; gmx grompp -f berendsen.mdp -c em_solvated.gro -p topol.top -o berendsen; cd ..; done

# apply grompp to all npt.mdp files for all replicas
for dir in sim{0..63}; do cd $dir; gmx grompp -f npt.mdp -c berendsen.gro -p topol.top -o PR_remd; cd ..; done

# change step number in npt.mdp files
for dir in sim{0..63}; do cd $dir; sed -i 's/nsteps\ =\ -1/nsteps\ =\ 55000000/' npt.mdp; cd ..; done

# Remove cpt files for all replicas
for dir in sim{0..63}; do cd $dir; echo $dir; rm npt.cpt; cd ..; done
for dir in sim{0..63}; do cd $dir; rm npt.cpt; cd ..; done

# Convert all .tpr files to xtc files
for dir in sim{0..63}; do cd $dir; echo $dir; for filename in npt*.trr; do echo $filename; gmx trjconv -f $filename -o "${filename%.*}.xtc"; done; cd ..; done

for dir in sim{0..63}; do cd $dir; echo $dir; for filename in npt*.trr; do echo $filename; gmx trjconv -f $filename -o "${filename%.*}.whole.xtc" -s npt.tpr -pbc whole; done; cd ..; done

# Updated for mpi version of gromacs
for dir in sim{0..63}; do cd $dir; echo $dir; for filename in npt*.trr; do echo $filename; mpirun -np 1 gmx_mpi trjconv -f $filename -o "${filename%.*}.whole.xtc" -s npt_new.tpr -pbc whole <<< 0; done; cd ..; done

# extend replica exhchange simulation
for dir in sim{0..63}; do cd $dir; pwd; mpirun -np 1 gmx_mpi convert-tpr -s npt_new.tpr -extend 50000 -o npt_new.tpr; cd ..; done

# Delete last 34 lines of a file and replace with tail of a file

for dir in sim{0..63}; do cd $dir; sed -i '34,$ d' npt_COM_pulling.mdp; tail -n 20 ../npt_COM_pulling.mdp >> npt_COM_pulling.mdp; cd ..; done

