import numpy as np
import os
import matplotliblib.pyplot as plt
import utils


utils.ReadItpFile("../OCT_dihes.itp")

def get_remd_dihedrals(torsion_ids, torsion_annotations, remd_traj, itp_file):
    torsion_matrix = []
    itp_reader = utils.ReadItpFile(itp_file)
    for i_torsions in range(len(torsion_annotations)):
        selector = itp_reader.construct_dihe_selection(torsion_annotations[i_torsions])
        dihes = utils.get_dihedrals(selector, remd_trajs.trajs[0])
        torsions_temp = []
        for i_rep in range(remd_traj.n_replicas):
            dihes_t = utils.get_dihedrals(selector, remd_traj.trajs[i_rep])
            torsions_temp.append(dihes_t)
        torsion_matrix.append(torsions_temp)
    return(torsion_matrix)
