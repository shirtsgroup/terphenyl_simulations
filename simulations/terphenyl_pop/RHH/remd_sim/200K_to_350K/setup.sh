for dir in sim{0..39}; do cd $dir; gmx grompp -f berendsen.mdp -c em_solvated.gro -p topol.top -o berendsen; cd ..; done
for dir in sim{0..39}; do cd $dir; gmx grompp -f npt.mdp -c berendsen.gro -p topol.top -o PR_remd; cd ..; done
