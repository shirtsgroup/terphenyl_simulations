import physical_validation
import mdtraj as md
import matplotlib.pyplot as plt

def main():
    parser = physical_validation.data.GromacsParser()
    res = parser.get_simulation_data(
        mdp = "sim0/nvt_COM_pulling_implicit_solv.mdp",
        top = "sim0/topol.top",
        gro = "sim0/nvt_COM_pulling_implicit_solv.part0001.xtc",
        edr = "sim0/nvt_COM_pulling_implicit_solv.part0001.edr"
    )


if __name__ == "__main__":
    main()