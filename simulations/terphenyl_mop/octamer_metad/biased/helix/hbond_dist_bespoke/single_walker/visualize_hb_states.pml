load ../npt_new.gro, hb_0
load ../npt_new.gro, hb_1
load ../npt_new.gro, hb_2
load ../npt_new.gro, hb_3
load ../npt_new.gro, hb_4
load ../npt_new.gro, hb_5
load ../npt_new.gro, hb_6
load ../npt_new.gro, hb_7

hide all
show sticks, resn OCT+CAP

load_traj hb_state_0.xtc, hb_0
load_traj hb_state_1.xtc, hb_1
load_traj hb_state_2.xtc, hb_2
load_traj hb_state_3.xtc, hb_3
load_traj hb_state_4.xtc, hb_4
load_traj hb_state_5.xtc, hb_5
load_traj hb_state_6.xtc, hb_6
load_traj hb_state_7.xtc, hb_7

intra_fit hb_0 and resn OCT+CAP 
intra_fit hb_1 and resn OCT+CAP 
intra_fit hb_2 and resn OCT+CAP 
intra_fit hb_3 and resn OCT+CAP 
intra_fit hb_4 and resn OCT+CAP 
intra_fit hb_5 and resn OCT+CAP 
intra_fit hb_6 and resn OCT+CAP 
intra_fit hb_7 and resn OCT+CAP 

alignto hb_0, align

center
