import mdtraj
import numpy as np


class ReadItpFile:
    """
    object used to extract various atom selections from an itp file
    """

    def __init__(self, filename, verbose = False):
        """
        ReadItpFile constructor, given an .itp file name
        """
        self.lines = []
        self.verbose = verbose
        with open(filename, 'r') as fh:
            for line in fh:
                self.lines.append(line.rstrip())

    def extract_annotated_atoms(self, annotation):
        atom_names = []
        for line in self.lines:
            if annotation in line.split():
                atom_names.append(line.split()[4])
        return(atom_names)

    def construct_dihe_selection(self, anote_list):
        assert len(anote_list) == 4
        atom_names = []
        for anote in anote_list:
            atoms = self.extract_annotated_atoms(anote)
            atom_names.append(atoms)
        dihe_names = []
        for a1, a2, a3, a4 in zip(*atom_names):
            if self.verbose:
                print(a1, a2, a3, a4)
            dihe_names.append(" ".join([a1, a2, a3, a4]))
        return " ".join(dihe_names)


def get_dihedrals(selector, traj, verbose = False):
    if verbose:
        print(selector)
    top = traj.topology
    i_dihes = []
    for atom in selector.split(' '):
        i_dihe = top.select("name "+atom)
        i_dihes.append(i_dihe)
    i_dihes = np.array(i_dihes)
    i_dihes = i_dihes.reshape(int(len(i_dihes)/4), 4)
    dihes = mdtraj.compute_dihedrals(traj, i_dihes, periodic=True).flatten()
    return(dihes)