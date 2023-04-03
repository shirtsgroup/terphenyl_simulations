# load visual of mop-terphenyl tetramer

load terphenyl_extended.gro
load_traj npt_new.whole.xtc

hide all
show sticks, resn TET+CAP
intra_fit resn TET+CAP

center