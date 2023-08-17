import terphenyl_simulations as ts
import mdtraj as md
import numpy as np


def main():
    gromacs_log = ts.utils.GromacsLogFile("sim0/npt_new.log")

    print(len(gromacs_log.states))


if __name__ == "__main__":
    main()