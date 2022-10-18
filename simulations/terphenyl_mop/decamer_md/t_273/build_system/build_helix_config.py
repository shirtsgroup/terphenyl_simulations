import heteropolymer_simulations as hs
import MDAnalysis as mda
import mdtraj as md
import numpy as np
import scipy as sp
import sys

def get_torsion_ids(universe, resname, torsion_id, template_residue_i = 0):
    """
    Using an MDAnalysis universe file with proper residue definitions, this function
    will extract the torsion ids of all torsions propagated along the chain. Specifically
    torsions should try to be fully defined within a single residue.
    """

    # Get index in residue
    atoms_in_residue = [a.name for a in universe.residues[template_residue_i].atoms]
    residue_atom_index  = [atoms_in_residue.index(a) for a in torsion_id if a in atoms_in_residue ]    
    dihedral_ids = []
    for residue in universe.residues:
        if residue.resname == resname:
            torsion_atoms = [residue.atoms[i] for i in residue_atom_index]
            dihedral_ids.append([ta.name for ta in torsion_atoms])
    
    return(dihedral_ids)


def main():

    # Torsion ids for torsions of interest

    hexamer_torsions_old = {
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


    decamer_u = mda.Universe("terphenyl_mop_decamer.itp", "em_decamer.gro")

    decamer_r1_torsions = {
        "a1" : ["C7", "C6", "C5", "C13"],
        "a2" : ["C5", "C13", "C14", "C15"],
        "p1" : ["C15", "C16", "C17", "N1"],
        "p2" : ["C16", "C17", "N1", "H15"],
        "p3" : ["O2", "C10", "C9", "C11"],
    }

    decamer_torsions = {}
    for torsion_type in decamer_r1_torsions.keys():
        t_ids = get_torsion_ids(decamer_u, "HEX", decamer_r1_torsions[torsion_type])
        decamer_torsions[torsion_type] = t_ids

    # Build ICE object
    ice = hs.edit_conf.InternalCoordinateEditor("em_decamer.gro", "terphenyl_mop_decamer.itp")
    
    # Load hexamer helical cluster
    print("Loading Cluster trajectory files...")
    cluster_traj = md.load("clustering_output/cluster_9.gro")


    for torsion_type in list(decamer_torsions.keys()):

        # if torsion_type == "p1":
        #   continue

        print("Working on torsion", torsion_type + "...")
        torsion_atom_names = hexamer_torsions_old[torsion_type]

        print("New atoms Torsion Atoms:")
        for t_id in decamer_torsions[torsion_type]:
            print(t_id)

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
        if torsion_type == "a2":
            torsion_max_density += - 25 * np.pi/180
        if torsion_type == "p2":
           torsion_max_density += -np.pi + 45 * np.pi/180
#        if torsion_type == "p3":
#           torsion_max_density += np.pi


        print("Setting", torsion_type, "to", torsion_max_density*180/np.pi)

        # Build InternalCoordinateEditor Object
        center_torsion_atoms = [torsion[1:3] for torsion in decamer_torsions[torsion_type]]
        skip = 0
        for i in range(len(center_torsion_atoms)):
            torsion_ids, torsions = ice.find_torsions(center_torsion_atoms[i], positions = [1, 2])
            prop_torsion_id = ice.identify_chain_prop_torsion(torsion_ids)
            set_non_prop_torsion(ice, prop_torsion_id, decamer_torsions[torsion_type][i], torsion_max_density)
            ice.update_internal_coordinates()
            # if torsion_type == "p1":
            #     if skip == 0:
            #        break
            #     else:
            #         skip += 1
    ice.write_structure("decamer_adjusted.gro")


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