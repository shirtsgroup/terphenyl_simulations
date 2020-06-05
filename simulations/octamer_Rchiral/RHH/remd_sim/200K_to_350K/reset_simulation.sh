#!/bin/bash
min_time=`for dir in sim{0..19}; do cd $dir; gmx check -f npt.cpt; cd ..; done 2>&1 | grep Last\ frame | awk '{print $5}' | awk 'BEGIN{a=100000}{if ($1<0+a) a =$1} END{print a}'`

echo $min_time

# for dir in sim{0..19}; do cd $dir; gmx trjconv -f npt.part0002.trr -s npt.tpr -o npt_restart.gro -dump $min_time <<< 0; cd ..; done

for dir in sim{0..19}; do cd $dir; gmx grompp -f npt.mdp -c npt_restart.gro -p topol.top -o npt; cd ..; done


