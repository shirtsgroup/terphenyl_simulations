# script for reading in and aligning multiple walker simulations in PyMOL

load WALKER0/terphenyl_extended.gro, walker_0
load WALKER1/terphenyl_extended.gro, walker_1
load WALKER2/terphenyl_extended.gro, walker_2
load WALKER3/terphenyl_extended.gro, walker_3

load_traj   WALKER0/npt_new.whole.xtc, walker_0
load_traj   WALKER1/npt_new.whole.xtc, walker_1
load_traj   WALKER2/npt_new.whole.xtc, walker_2
load_traj   WALKER3/npt_new.whole.xtc, walker_3

intra_fit resn TET+CAP and walker_0
intra_fit resn TET+CAP and walker_1
intra_fit resn TET+CAP and walker_2
intra_fit resn TET+CAP and walker_3

hide sticks, resn TCM

alignto walker_0, align

translate [-40, 0, 0], walker_0, 0
translate [-20, 0, 0], walker_1, 0
translate [0,   0, 0], walker_2, 0
translate [20,  0, 0], walker_3, 0

center