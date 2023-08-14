import pymol
import glob


def main():
    n_simulations = len(glob.glob("sim*"))

    pymol.finish_launching(["pymol"])
    
    # Load individual structure files
    for i in range(n_simulations):
        pymol.cmd.load("sim" + str(i) + "/berendsen_npt.gro", "sim" + str(i))
        pymol.cmd.load_traj("sim"+ str(i) + "/npt_new.whole.xtc", "sim" + str(i), interval = 20)
        pymol.cmd.intra_fit("not resn SOL+K+F and sim" + str(i))

    pymol.cmd.hide("all")
    pymol.cmd.show("sticks", "not resn SOL+K+F")
    pymol.cmd.show("cartoon", "not resn SOL+K+F")
    pymol.cmd.alignto("sim0", "align")
   


if __name__ == "__main__":
    main() 
