# script for reading in and aligning multiple walker simulations in PyMOL

load WALKER0/em_solvated_helix.gro, walker_0
load WALKER1/em_solvated_helix.gro, walker_1
load WALKER2/em_solvated_helix.gro, walker_2
load WALKER3/em_solvated_helix.gro, walker_3
load WALKER4/em_solvated_helix.gro, walker_4
load WALKER5/em_solvated_helix.gro, walker_5
load WALKER6/em_solvated_helix.gro, walker_6
load WALKER7/em_solvated_helix.gro, walker_7
load WALKER8/em_solvated_helix.gro, walker_8
load WALKER9/em_solvated_helix.gro, walker_9
load WALKER10/em_solvated_helix.gro, walker_10
load WALKER11/em_solvated_helix.gro, walker_11
load WALKER12/em_solvated_helix.gro, walker_12
load WALKER13/em_solvated_helix.gro, walker_13
load WALKER14/em_solvated_helix.gro, walker_14
load WALKER15/em_solvated_helix.gro, walker_15

load_traj   WALKER0/npt_new.whole.xtc, walker_0,interval=10
load_traj   WALKER1/npt_new.whole.xtc, walker_1,interval=10
load_traj   WALKER2/npt_new.whole.xtc, walker_2,interval=10
load_traj   WALKER3/npt_new.whole.xtc, walker_3,interval=10
load_traj   WALKER4/npt_new.whole.xtc, walker_4,interval=10
load_traj   WALKER5/npt_new.whole.xtc, walker_5,interval=10
load_traj   WALKER6/npt_new.whole.xtc, walker_6,interval=10
load_traj   WALKER7/npt_new.whole.xtc, walker_7,interval=10
load_traj   WALKER8/npt_new.whole.xtc, walker_8,interval=10
load_traj   WALKER9/npt_new.whole.xtc, walker_9,interval=10
load_traj   WALKER10/npt_new.whole.xtc, walker_10,interval=10
load_traj   WALKER11/npt_new.whole.xtc, walker_11,interval=10
load_traj   WALKER12/npt_new.whole.xtc, walker_12,interval=10
load_traj   WALKER13/npt_new.whole.xtc, walker_13,interval=10
load_traj   WALKER14/npt_new.whole.xtc, walker_14,interval=10
load_traj   WALKER15/npt_new.whole.xtc, walker_15,interval=10

intra_fit resn OCT+CAP and walker_0
intra_fit resn OCT+CAP and walker_1
intra_fit resn OCT+CAP and walker_2
intra_fit resn OCT+CAP and walker_3
intra_fit resn OCT+CAP and walker_4
intra_fit resn OCT+CAP and walker_5
intra_fit resn OCT+CAP and walker_6
intra_fit resn OCT+CAP and walker_7
intra_fit resn OCT+CAP and walker_8
intra_fit resn OCT+CAP and walker_9
intra_fit resn OCT+CAP and walker_10
intra_fit resn OCT+CAP and walker_11
intra_fit resn OCT+CAP and walker_12
intra_fit resn OCT+CAP and walker_13
intra_fit resn OCT+CAP and walker_14
intra_fit resn OCT+CAP and walker_15

hide sticks, resn TCM

alignto walker_0, align