import heteropolymer_simulations as hs
import MDAnalysis as mda
import mdtraj as md
import numpy as np
import scipy as sp#
import sys
import yaml

def main():
    # Torsion ids for torsions of interest

    with open("old_hexamer_torsion_ids.yml", 'r') as yaml_file:
        old_hexamer_torsions = yaml.safe_load(yaml_file)

    octamer_u = mda.Universe("terphenyl_mop_octamer.itp", "em_octamer.gro")

    octamer_r1_torsions = {
        "a1" : ["C4", "C5", "C6", "C7"],
        "a2" : ["C6", "C7", "C13", "C14"],
        "p1" : ["C18", "C17", "C19", "N1"],
        "p2" : ["C17", "C19", "N1", "H14"],
        "p3" : ["O1", "C1", "C2", "C3"],
    }

    # Get all torsions defined in residue 1
    

    octamer_torsions = {}
    for torsion_type in octamer_r1_torsions.keys():
        t_ids = hs.utils.get_torsion_ids(octamer_u, "OCT", octamer_r1_torsions[torsion_type])
        octamer_torsions[torsion_type] = t_ids

    # Build ICE object
    ice = hs.edit_conf.InternalCoordinateEditor("em_octamer.gro", "terphenyl_mop_octamer.itp")
    
    # Load hexamer helical cluster
    print("Loading Cluster trajectory files...")
    cluster_traj = md.load("cluster_9.gro")


    for torsion_type in list(octamer_torsions.keys()):

        # if torsion_type == "p1":
        #   continue

        print("Working on torsion", torsion_type + "...")
        torsion_atom_names = old_hexamer_torsions[torsion_type]

        print("Old Torsion Atoms:")
        for t_id in torsion_atom_names:
            print(t_id)

        print("New atoms Torsion Atoms:")
        for t_id in octamer_torsions[torsion_type]:
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

        #if torsion_type == "a1":
        #    torsion_max_density -= np.pi
        #if torsion_type == "a2":
        #    torsion_max_density += np.pi - 25 * np.pi/180
        #if torsion_type == "p2":
        #   torsion_max_density += -np.pi + 45 * np.pi/180
        #if torsion_type == "p3":
        #   torsion_max_density += np.pi


        print("Setting", torsion_type, "to", torsion_max_density*180/np.pi)

        # Build InternalCoordinateEditor Object
        center_torsion_atoms = [torsion[1:3] for torsion in octamer_torsions[torsion_type]]
        skip = 0
        for i in range(len(center_torsion_atoms)):
            torsion_ids, torsions = ice.find_torsions(center_torsion_atoms[i], positions = [1, 2])
            prop_torsion_id = ice.identify_chain_prop_torsion(torsion_ids)
            set_non_prop_torsion(ice, prop_torsion_id, octamer_torsions[torsion_type][i], torsion_max_density)
            # if torsion_type == "p1":
            #     if skip == 0:
            #        break
            #     else:
            #         skip += 1
    
    shift = ice.torsions[ice.bat._primary_torsion_indices]
    shift[ice.bat._unique_primary_torsion_indices] = 0
    ice.torsions -= shift
    ice.update_internal_coordinates()
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