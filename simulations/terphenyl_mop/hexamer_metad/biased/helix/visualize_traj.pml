# Read in trajectory of individaul jobs

load em_solvated_helix.gro, traj
load_traj npt_new.whole.xtc, traj

hide all
show sticks, resn HEX+CAP
intra_fit resn HEX+CAP

center