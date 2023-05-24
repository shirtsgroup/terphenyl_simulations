#!/bin/bash

# Run berendsenv NVT equilibration simulation
gmx grompp -f berendsen_nvt.mdp -c em_solvated_helix.gro -p system_hmr_rest.top -o berendsen_nvt -maxwarn 2 -r em_solvated_helix.gro
gmx mdrun -v -deffnm berendsen_nvt 

# Run berendsenv NPT equilibration simulation
gmx grompp -f berendsen_npt.mdp -c berendsen_nvt.gro -p system_hmr_rest.top -o berendsen_npt -maxwarn 2 -r em_solvated_helix.gro
gmx mdrun -v -deffnm berendsen_npt 

# Run production Metadynamics simulation
gmx grompp -f npt_new.mdp -c berendsen_npt.gro -p system_hmr.top -o npt_new
gmx mdrun -v -deffnm npt_new -plumed plumed_hbond_dist.dat