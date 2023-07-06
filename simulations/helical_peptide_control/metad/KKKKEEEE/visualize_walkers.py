import pymol
import glob

def main():
    n_simulations = len(glob.glob("WALKER*"))

    pymol.finish_launching(["pymol"])
    
    # Load individual structure files
    for i in range(n_simulations):
        pymol.cmd.load("WALKER" + str(i) + "/berendsen_npt.gro", "walker_" + str(i))
        pymol.cmd.load_traj("WALKER" + str(i) + "/npt_new.whole.xtc", "walker_" + str(i), interval = 10)
        pymol.cmd.intra_fit("not resn SOL+K+F and walker_" + str(i))

    pymol.cmd.hide("all")
    pymol.cmd.show("sticks", "not resn SOL+K+F")
    pymol.cmd.show("cartoon", "not resn SOL+K+F")
    pymol.cmd.alignto("walker_0", "align")
    
if __name__ == "__main__":
    main()    
