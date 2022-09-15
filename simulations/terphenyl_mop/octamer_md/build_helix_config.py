import heteropolymer_simulations as hs
import MDAnalysis as mda
import mdtraj as md
import numpy as np
import scipy as sp
import sys


def main():

    # Torsion ids for torsions of interest

    hexamer_torsions = {
        "a1" : [
            ["C131", "C66", "C61", "C65"],
            ["C70", "C70", "C33", "C32"],
            ["C5", "C7", "C99", "C25"],
            ["C123", "C50", "C44", "C11"],
            ["C24", "C106", "C100", "C105"],
            ["C91", "C86", "C83", "C84"],
        ],
        "a2" : [
            ["C61", "C65", "C74", "C73"],
            ["C33", "C32", "C72", "C9"],
            ["C99", "C57", "C53", "C1"],
            ["C44", "C11", "C47", "C48"],
            ["C100", "C105", "C108", "C109"],
            ["C83", "C84", "C92", "C93"],
        ],
        "p1" : [
            ["C73", "C130", "C136", "N6"],
            ["C9", "C114", "C116", "N1"],
            ["C1", "C56", "C31", "N3"],
            ["C48", "C35", "C124", "N4"],
            ["C21", "C15", "C36", "N2"],
            ["C95", "C94", "C98", "N5"],
        ],
        "p2" : [
            ["C130", "C136", "N6", "C28"],
            ["C114", "C116", "N1", "C43"],
            ["C56", "C31", "N3", "C2"],
            ["C35", "C124", "N4", "C37"],
            ["C15", "C36", "N2", "C97"],
            ["C94", "C98", "N5", "C16"],
        ],
        "p3" : [
            ["N6", "C28", "C68", "C20"],
            ["N1", "C43", "C119", "C117"],
            ["N3", "C2", "C8", "C34"],
            ["N4", "C37", "C40", "C39"],
            ["N2", "C97", "C89", "C88"],
            ["N3", "C2", "C8", "C34"],
        ],
    }

    octamer_torsions = {
        "a1" : [
            ["C101", "C174", "C168", "C173"],
            ["C160", "C155", "C152", "C153"],
            ["C121", "C114", "C108", "C113"],
            ["C140", "C135", "C132", "C133"],
            ["C10", "C86", "C80", "C85"],
            ["C72", "C67", "C64", "C65"],
            ["C33", "C26", "C20", "C25"],
            ["C52", "C47", "C44", "C45"], 
        ],
        "a2" : [
            ["C168", "C173", "C176", "C177"],
            ["C152", "C153", "C161", "C162"],
            ["C108", "C113", "C116", "C117"],
            ["C132", "C133", "C141", "C142"],
            ["C80", "C85", "C88", "C89"],
            ["C64", "C65", "C73", "C74"],
            ["C20", "C25", "C28", "C29"],
            ["C44", "C45", "C53", "C54"], 
        ],
        "p1" : [
            ["C177", "C98", "C104", "N5"],
            ["C162", "C163", "C167", "N8"],
            ["C117", "C118", "C126", "N6"],
            ["C142", "C143", "C147", "N7"],
            ["C89", "C4", "C15", "N1"],
            ["C74", "C75", "C79", "N4"],
            ["C29", "C30", "C38", "N2"],
            ["C54", "C55", "C59", "N3"], 
        ],
        "p2" : [
            ["C98", "C104", "N5", "C166"],
            ["C163", "C167", "N8", "C127"],
            ["C118", "C126", "N6", "C146"],
            ["C143", "C147", "N7", "C16"],
            ["C4", "C15", "N1", "C78"],
            ["C75", "C79", "N4", "C39"],
            ["C30", "C38", "N2", "C58"],
            ["C55", "C59", "N3", "C5"], 
        ],
        "p3" : [
            ["N5", "C166", "C158", "C159"],
            ["N8", "C127", "C123", "C122"],
            ["N6", "C146", "C138", "C139"],
            ["N7", "C16", "C19", "C13"],
            ["N1", "C78", "C70", "C71"],
            ["N4", "C39", "C35", "C34"],
            ["N2", "C58", "C50", "C51"],
            ["O19", "C105", "C107", "C103"], 
        ],
    }

    # Build ICE object
    ice = hs.edit_conf.InternalCoordinateEditor("mop_octamer.gro", "terphenyl_mop_octamer.itp")
    
    # Load hexamer helical cluster
    print("Loading Cluster trajectory files...")
    cluster_traj = md.load("clustering_output/cluster_9.gro")


    for torsion_type in list(hexamer_torsions.keys()):

        # if torsion_type == "p1":
        #   continue

        print("Working on torsion", torsion_type + "...")
        torsion_atom_names = hexamer_torsions[torsion_type]

        mirror_sym = False
        if torsion_type == "a1":
            mirror_sym = True

        hs.plotting.plot_torsions_distributions(cluster_traj,
                                torsion_atom_names,
                                torsion_type,
                                torsion_type,
                                torsion_type,
                                legend = None,
                                mirror_sym = mirror_sym
                                )
        
        # Get torsions from cluster of interest
        torsions = hs.observables.get_torsions(cluster_traj, torsion_atom_names, mirror_sym = mirror_sym)

        # Histogram torsions
        bin_edges = np.linspace(-np.pi, np.pi, 50 + 1)
        bin_centers = np.array([(bin_edges[i] + bin_edges[i+1]) * 0.5 for i in range(len(bin_edges) - 1) ])
        hist, bin_edges_out = np.histogram(np.array(torsions), bins = bin_edges, density = True)
        # For p2 take the 2nd largest peak
        torsion_max_density = bin_centers[np.argmax(hist)]

        if torsion_type == "a1":
            torsion_max_density -= np.pi
        
        # if torsion_type == "p1":
        #    torsion_max_density += 50 * np.pi / 180


        print("Setting", torsion_type, "to", torsion_max_density*180/np.pi)

        # Build InternalCoordinateEditor Object
        center_torsion_atoms = [torsion[1:3] for torsion in octamer_torsions[torsion_type]]
        skip = 0
        for i in range(len(center_torsion_atoms)):
            torsion_ids, torsions = ice.find_torsions(center_torsion_atoms[i], positions = [1, 2])
            torsion_id = ice.identify_chain_prop_torsion(torsion_ids)
            set_non_prop_torsion(ice, torsion_id, octamer_torsions[torsion_type][i], torsion_max_density)
            ice.update_internal_coordinates()
            if torsion_type == "p1":
                if skip == 0:
                    break
                else:
                    skip += 1
    ice.write_structure("octamer_adjusted.gro")


def set_non_prop_torsion(ice_obj, prop_torsion_id, non_prop_torsion, value):
    print("Target torsion:", value * 180 / np.pi)
    def loss(x):
        ice_obj.set_torsion(prop_torsion_id, x[0])
        ice_obj.update_internal_coordinates()
        atom_coords = [ice_obj.universe.select_atoms("name " + atom)[0].position for atom in non_prop_torsion]
        torsion = mda.lib.distances.calc_dihedrals(*atom_coords)
        return (value - torsion) ** 2

    sp.optimize.minimize(loss, np.array([value]), options={'eps':0.00001})
    atom_coords = [ice_obj.universe.select_atoms("name " + atom)[0].position for atom in non_prop_torsion]
    torsion = mda.lib.distances.calc_dihedrals(*atom_coords)
    print("Final torsion:",  torsion * 180 / np.pi)
    # print(res)



if __name__ == "__main__":
    main()