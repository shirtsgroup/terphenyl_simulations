import heteropolymer_simulations as hs
import mdtraj as md
import numpy as np



def main():
    # Get high probability torsion values for torsions of interest
    print("Loading Cluster trajectory files...")
    cluster_traj = md.load("clustering_output/cluster_22.gro")


    # A1 Torsion
    print("Working on Aromatic Torsions 1...")
    torsion_atom_names = [
        ["C133", "C66", "C61", "C65"],
        ["C69", "C70", "C33", "C32"],
        ["C111", "C7", "C99", "C25"],
        ["C123", "C50", "C44", "C45"],
        ["C24", "C106", "C100", "C105"],
        ["C87", "C86", "C83", "C82"],
    ]
    
    # Get torsions from cluster of interest
    torsions = hs.observables.get_torsions(cluster_traj, torsion_atom_names)

    # Histogram torsions
    bin_edges = np.linspace(-np.pi, np.pi, 50 + 1)
    bin_centers = np.array([(bin_edges[i] + bin_edges[i+1]) * 0.5 for i in range(len(bin_edges) - 1) ])
    hist, bin_edges_out = np.histogram(np.array(torsions), bins = bin_edges, density = True)
    max_density = bin_centers[np.argmax(hist)]

    # Build InternalCoordinateEditor Object
    ice = hs.edit_conf.InternalCoordinateEditor("mop_octamer.gro", "terphenyl_mop_octamer.itp")

    octamer_a1_pairs = [
        ["C174", "C168"],
        ["C155", "C152"],
        ["C108", "C114"],
        ["C132", "C135"],
        ["C86", "C80"],
        ["C64", "C67"],
        ["C20", "C26"],
        ["C47", "C44"],    
    ]

    all_torsion_ids = []
    all_torsions = []

    for atom_pair in octamer_a1_pairs:
        torsion_ids, torsions = ice.find_torsions(atom_pair, set_op = "and")
        
        # List comprehension black magic
        # filter_torsion_ids = [
        #     t_id for i_id in torsion_ids if atom_pair[0] not in t_id.split[0] \
        #                                     and atom_pair[0] not in t_id.split[3] \
        #                                     and atom_pair[1] not in t_id.split[0] \
        #                                     and atom_pair[1] not in t_id.split[3] \                                                            
        # ]

        # filter_torsions = [
        #     t_id for i_id in torsion_ids if atom_pair[0] not in t_id.split[0] \
        #                                     and atom_pair[0] not in t_id.split[3] \
        #                                     and atom_pair[1] not in t_id.split[0] \
        #                                     and atom_pair[1] not in t_id.split[3] \                                                            
        # ]

        for t_id, t in zip(torsion_ids, torsions):
            if atom_pair[0] not in t_id.split(" ")[0] and atom_pair[0] not in t_id.split(" ")[3]:
                if atom_pair[1] not in t_id.split(" ")[0] and atom_pair[1] not in t_id.split(" ")[3]:
                    all_torsion_ids.append(t_id)
                    all_torsions.append(t)
        
    print(all_torsion_ids)
    print(all_torsions)



if __name__ == "__main__":
    main()