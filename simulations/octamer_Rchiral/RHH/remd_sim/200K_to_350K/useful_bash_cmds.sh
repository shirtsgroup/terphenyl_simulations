# apply grompp to berendsen.mdp files for all replicas
for dir in sim{0..39}; do cd $dir; gmx grompp -f berendsen.mdp -c em_solvated.gro -p topol.top -o berendsen; cd ..; done

# apply grompp to all npt.mdp files for all replicas
for dir in sim{0..39}; do cd $dir; gmx grompp -f npt.mdp -c berendsen.gro -p topol.top -o PR_remd; cd ..; done

# change step number in npt.mdp files
for dir in sim{0..39}; do cd $dir; sed -i 's/nsteps\ =\ -1/nsteps\ =\ 55000000/' npt.mdp; cd ..; done

# Remove cpt files for all replicas
for dir in sim{0..39}; do cd $dir; echo $dir; rm npt.cpt; cd ..; done
for dir in sim{0..19}; do cd $dir; rm npt.cpt; cd ..; done

# Convert all .tpr files to xtc files
for dir in sim{0..39}; do cd $dir; echo $dir; for filename in npt*.trr; do echo $filename; gmx trjconv -f $filename -o "${filename%.*}.xtc"; done; cd ..; done

for dir in sim{0..39}; do cd $dir; echo $dir; for filename in npt*.trr; do echo $filename; gmx trjconv -f $filename -o "${filename%.*}.whole.xtc" -s npt.tpr -pbc whole; done; cd ..; done
