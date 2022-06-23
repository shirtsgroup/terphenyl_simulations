import physical_validation
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt

def main():
    t1 = time.time()

    # physical_validation.util.gromacs_interface.GromacsInterface(exe = "")

    parser = physical_validation.data.GromacsParser()
    res = parser.get_simulation_data(
        mdp = "mdout.mdp",
        top = "topol.top",
        trr = "nvt_COM_pulling_implicit_solv.part0001.trr",
        edr = "nvt_COM_pulling_implicit_solv.part0001.edr"
    )

    res.system.ndof_reduction_tra = 0

    results = physical_validation.kinetic_energy.equipartition(
        data=res, 
        strict=False, 
        filename='ke_equipartition',
    )

    print(results)

    t2 = time.time()
    print(f'Time elapsed: {t2-t1} seconds.')


if __name__ == "__main__":
    main()