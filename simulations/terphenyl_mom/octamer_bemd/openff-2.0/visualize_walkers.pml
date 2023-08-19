# script for reading in and aligning multiple walker simulations in PyMOL

load WALKER0/berendsen_npt.gro, walker_0
load WALKER1/berendsen_npt.gro, walker_1
load WALKER2/berendsen_npt.gro, walker_2
load WALKER3/berendsen_npt.gro, walker_3
load WALKER4/berendsen_npt.gro, walker_4

load_traj   WALKER0/npt_new.whole.xtc, walker_0,interval=1
load_traj   WALKER1/npt_new.whole.xtc, walker_1,interval=1
load_traj   WALKER2/npt_new.whole.xtc, walker_2,interval=1
load_traj   WALKER3/npt_new.whole.xtc, walker_3,interval=1
load_traj   WALKER4/npt_new.whole.xtc, walker_4,interval=1

intra_fit resn OCT+CAP and walker_0
intra_fit resn OCT+CAP and walker_1
intra_fit resn OCT+CAP and walker_2
intra_fit resn OCT+CAP and walker_3
intra_fit resn OCT+CAP and walker_4

hide sticks, resn TCM

alignto walker_0, align

translate [30, 0, 0], walker_1, 0
translate [60, 0, 0], walker_2, 0
translate [90, 0, 0], walker_3, 0
translate [120, 0, 0], walker_4, 0
