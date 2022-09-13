import heteropolymer_simulations as hs
import mdtraj as md
import numpy as np


def main():

    # Torsion ids for torsions of interest

    hexamer_torsions = {
        "a1" : [
            ["C133", "C66", "C61", "C65"],
            ["C69", "C70", "C33", "C32"],
            ["C111", "C7", "C99", "C25"],
            ["C123", "C50", "C44", "C45"],
            ["C24", "C106", "C100", "C105"],
            ["C87", "C86", "C83", "C82"],
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

    octamer_torsion_bonds = {
        "a1" : [
            ["C174", "C168"],
            ["C155", "C152"],
            ["C108", "C114"],
            ["C132", "C135"],
            ["C86", "C80"],
            ["C64", "C67"],
            ["C20", "C26"],
            ["C47", "C44"], 
        ],
        "a2" : [
            ["C176", "C173"],
            ["C161", "C153"],
            ["C113", "C116"],
            ["C133", "C141"],
            ["C85", "C88"],
            ["C73", "C65"],
            ["C25", "C28"],
            ["C53", "C45"], 
        ],
        "p1" : [
            ["C104", "C98"],
            ["C163", "C167"],
            ["C126", "C118"],
            ["C147", "C143"],
            ["C15", "C4"],
            ["C79", "C75"],
            ["C38", "C30"],
            ["C59", "C55"], 
        ],
        "p2" : [
            ["C104", "N5"],
            ["C167", "N8"],
            ["C126", "N6"],
            ["C147", "N7"],
            ["C15", "N1"],
            ["C79", "N4"],
            ["C38", "N2"],
            ["C59", "N3"], 
        ],
        "p3" : [
            ["C158", "C166"],
            ["C123", "C127"],
            ["C138", "C146"],
            ["C16", "C19"],
            ["C70", "C78"],
            ["C39", "C35"],
            ["C50", "C58"],
            ["C107", "C105"], 
        ],
    }

    # Build ICE object
    ice = hs.edit_conf.InternalCoordinateEditor("mop_octamer.gro", "terphenyl_mop_octamer.itp")
    
    # Load hexamer helical cluster
    print("Loading Cluster trajectory files...")
    cluster_traj = md.load("clustering_output/cluster_22.gro")


    for torsion_type in list(hexamer_torsions.keys())[2:4]:

        print("Working on torsion", torsion_type + "...")
        print(hexamer_torsions[torsion_type])
        torsion_atom_names = hexamer_torsions[torsion_type]
        
        # Get torsions from cluster of interest
        torsions = hs.observables.get_torsions(cluster_traj, torsion_atom_names)

        # Histogram torsions
        bin_edges = np.linspace(-np.pi, np.pi, 50 + 1)
        bin_centers = np.array([(bin_edges[i] + bin_edges[i+1]) * 0.5 for i in range(len(bin_edges) - 1) ])
        hist, bin_edges_out = np.histogram(np.array(torsions), bins = bin_edges, density = True)
        torsion_max_density = bin_centers[np.argmax(hist)]

        print("Setting torsion", torsion_type, "to", torsion_max_density * 180 / np.pi)

        # Build InternalCoordinateEditor Object
        octamer_a1_pairs = octamer_torsion_bonds[torsion_type]

        for atom_pair in octamer_a1_pairs:
            torsion_ids, torsions = ice.find_torsions(atom_pair, positions = [1, 2])
            torsion_id = ice.identify_chain_prop_torsion(torsion_ids)
            ice.set_torsion(torsion_id, torsion_max_density)
            ice.update_internal_coordinates()
    
    ice.write_structure("octamer_adjusted.gro")



if __name__ == "__main__":
    main()