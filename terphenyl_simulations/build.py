import yaml
import os
import sys
import numpy as np
from abc import ABC, abstractclassmethod
from subprocess import Popen, PIPE
import warnings
from openbabel import openbabel
from mbuild import load
from mbuild.lib.recipes.polymer import Polymer
from mbuild.utils.geometry import calc_dihedral
from openff.toolkit import ForceField
from openff.toolkit.topology import Molecule, Topology
from openff.interchange import Interchange
from openff.interchange.components._packmol import pack_box, UNIT_CUBE
from openff.units import unit
from .assign_parameters import FoldamerOFFDefault, FoldamerOFFBespoke, SystemOFFDefault
from .gromacs_wrapper import GromacsWrapper
from .topology_manager import TopologyManager
from .utils import ROOT_DIR, replace_all_pattern, make_path, renumber_pdb_atoms, get_solvent_structure_file


class FoldamerBuilder:
    """
    The FoldamerBuilder object is used to build foldamers with MBuild based on
    instructions provided in a foldamer build YML file. These YML input files
    provide the monomer unit SMILE string, connection atoms, and capping residues.
    It may take some hand-tuning of the Polymer object's parameters to get the
    foldamer be chemically correct. These parameters can be modified in the input
    build.yml file.

    Parameters
    ----------
    build_file_yml : str
        String specifying the build yaml file.
    path : str
        Path to directory where to write the output files
    """

    def __init__(self, build_file_yml, path="", topology_manager=TopologyManager()):
        self.build_file = build_file_yml
        self.label = "molecule"
        self.topology_manager = topology_manager
        with open(build_file_yml, "r") as f:
            self.build_params = yaml.safe_load(f)
        self.path = path
        if not os.path.isdir(self.path):
            make_path(path)

    def get_foldamer(self):
        # Check DB of entries first
        if self.topology_manager.check_file_type(self.build_file, self.label, "pdb") \
             and self.topology_manager.check_file_type(self.build_file, self.label, "mol"):
            print("Using database structure_file...")
            self.topology_manager.get_structure(
                self.build_file, self.label, self.path, filetype="pdb"
            )
            self.topology_manager.get_structure(
                self.build_file, self.label, self.path, filetype="mol"
            )
        else:
            print("Building Foldamer from", self.build_file + "...")
            self.build_foldamer()
            self.fix_peptide_bonds()
            self.write_pdb()
            self.write_mol()

    def build_foldamer(self):
        smile_str = self.build_params["foldamer_smile"]
        subunit = load(
            smile_str, smiles=True, name=self.build_params["residue_name"]
        )
        head_cap = load(self.build_params["cap_smiles"]["upper"], smiles=True)
        tail_cap = load(self.build_params["cap_smiles"]["lower"], smiles=True)

        # We use the MBuild Polymer object to build our foldamer model
        self.chain = Polymer()
        self.chain.add_monomer(
            compound=subunit,
            indices=[
                self.build_params["upper_connect"],
                self.build_params["lower_connect"],
            ],
            separation=0.15,
            replace=True,
            # orientation = [[0,-1,0],[1,0,0]]
        )
        self.chain.add_end_groups(
            compound=head_cap, index=-1, separation=0.15, label="head", duplicate=False
        )

        self.chain.add_end_groups(
            compound=tail_cap, index=-1, separation=0.15, label="tail", duplicate=False
        )
        self.chain.build(n=self.build_params["foldamer_length"], sequence="A")
        self.chain.name = self.build_params["residue_name"]

        # Adjust peptide bonds to be trans
        for bond in self.chain.bonds(return_bond_order=True):
            atom_1 = bond[0]
            atom_2 = bond[1]

            if (atom_1.name == "C" and atom_2.name == "N") or (atom_2.name == "C" and atom_1.name == "N"):
                if atom_1.n_direct_bonds == 3 and atom_2.n_direct_bonds == 3:
                    bonded_1 = list(atom_1.direct_bonds())
                    n_bonded_1 = [atom.n_direct_bonds for atom in bonded_1]
                    bonded_2 = list(atom_2.direct_bonds())
                    n_bonded_2 = [atom.n_direct_bonds for atom in bonded_2]

                    # Dihedral between singly bonded H and O
                    hydro = bonded_1[n_bonded_1.index(1)]
                    carboxyl = bonded_2[n_bonded_2.index(1)]

                    # Get dihedral of peptide bond
                    dihe = calc_dihedral(hydro.pos, atom_1.pos, atom_2.pos, carboxyl.pos)
                    adjust = np.pi - dihe

                    # This corrects dihedrals to be trans
                    self.chain.rotate_dihedral(bond[0:2], adjust)

                    # Minimize after each adjustment
                    # self.chain.energy_minimize(forcefield="MMFF94")


        self.chain.save("test.pdb", overwrite=True)
        self.chain.energy_minimize(forcefield="MMFF94")

        # Change residue names in chain object
        for label in self.chain.labels["monomer"]:
            label.name = self.build_params["residue_name"]
        for label in self.chain.labels["Compound"]:
            label.name = "CAP"

    def fix_peptide_bonds(self):
        pass

    def write_pdb(self):
        filename = os.path.join(self.path, self.build_params["structure_file"] + ".pdb")
        self.chain.save(
            filename,
            overwrite=True,
            residues=[self.build_params["residue_name"], "CAP"],
        )
        self.topology_manager.add_structure(filename, self.build_file, self.label)

    def write_mol(self):
        pdb_fn = os.path.join(self.path, self.build_params["structure_file"] + ".pdb")
        mol_fn = os.path.join(self.path, self.build_params["structure_file"] + ".mol")
        if not os.path.exists(pdb_fn):
            self.write_pdb()
        ob_convert = openbabel.OBConversion()
        ob_convert.SetInAndOutFormats("pdb", "mol")
        mol = openbabel.OBMol()
        ob_convert.ReadFile(mol, pdb_fn)
        ob_convert.WriteFile(mol, mol_fn)
        self.topology_manager.add_structure(mol_fn, self.build_file, self.label)


class SystemBuilderOpenFF:
    """
    """

    def __init__(
        self,
        solute_pdb,
        solvent_id,
        build_file_yml,
        path="",
        topology_manager=TopologyManager(),
        label = "system"
    ):
        self.build_file = build_file_yml
        self.topology_manager = topology_manager
        self.label = label
        self.path = path
        with open(build_file_yml, "r") as f:
            self.build_params = yaml.safe_load(f)

        self.filename = self.build_params["structure_file"] + "_" + self.label

        self.solute = Molecule.from_file(solute_pdb)
        self.solvent = Molecule.from_file(get_solvent_structure_file(solvent_id))
        self.md_engine = GromacsWrapper()

    def get_system_structure(self):
        if self.topology_manager.check_file_type(
            self.build_file, self.label, "gro"
        ) and self.topology_manager.check_file_type(
            self.build_file, self.label, "pdb"
        ):
            print("Using database structure_file...")   
            self.gro_file = self.topology_manager.get_structure(
                self.build_file, self.label, self.path, filetype="gro"
            )
            self.pdb_file = self.topology_manager.get_structure(
                self.build_file, self.label, self.path, filetype="pdb"
            )
        else:
            self.build_system()
            self.write_pdb()
            self.write_gro()

    def build_system(self):
        self.system = pack_box(
            molecules=[self.solute, self.solvent],
            number_of_copies=[1, self.build_params["system"]["n_solvent"]],
            box_vectors=self.build_params["system"]["box_size"]
            * UNIT_CUBE
            * unit.angstrom,
        )

    def write_pdb(self):
        self.pdb_file = self.filename + ".pdb"
        self.system.to_file(self.pdb_file)
        self.topology_manager.add_structure(
            self.pdb_file, self.build_file, label = self.label
        )

    def write_gro(self):
        self.gro_file = self.filename + ".gro"
        gmx = GromacsWrapper()
        gmx.edit_conf(f = self.filename + ".pdb", o = self.filename + ".gro")
        self.topology_manager.add_structure(
            self.gro_file, self.build_file, label = self.label
        )



    # def minimize_system(self):
    #     centered_gro = self.gro_file.split(".gro")[0] + "_box.gro"
    #     self.md_engine.center_configuration(
    #         self.gro_file,
    #         centered_gro,
    #     )
    #     self.md_engine.minimize(centered_gro, self.top_file, prefix="em_solvated")
    #     self.gro_file = "em_solvated.gro"


# I probably best to get rid of this class
class MoleculeTopologyGenerator:
    def __init__(
        self,
        molecule_file,
        pdb_file,
        build_file,
        ff_method,
        path="",
        ff_name="openff-2.0.0",
        topology_manager=TopologyManager(),
        topology_label="molecule",
    ):
        self.path = path
        self.build_file = build_file
        self.topology_manager = topology_manager
        self.topology_label = topology_label
        self.ff_method = ff_method
        if not os.path.isdir(self.path):
            make_path(path)

        self._ff_generation_methods = {
            "openff": FoldamerOFFDefault,
            "bespoke": FoldamerOFFBespoke,
        }

        if ff_method in self._ff_generation_methods.keys():
            self.ff_generator = self._ff_generation_methods[ff_method](
                molecule_file, pdb_file, path=self.path, ff_str=ff_name, build_file = self.build_file, topology_manager = self.topology_manager
            )
        else:
            warnings.warn(
                "WARNING: "
                + ff_method
                + " is not one of the available "
                + "force field parameter generation methods. Please pick from:\n"
                + " ".join(self._ff_generation_methods.keys())
            )
            sys.exit(1)

        # Define other attributes populated by other functions
        self.md_engine = None
        self.top_file = None
        self.gro_file = None

    def get_ff_parameters(self):
        if self.topology_manager.check_file_type(
            self.build_file, self.topology_label, "top"
        ):
            self.top_file = self.topology_manager.get_topology(
                self.build_file, self.topology_label, self.path
            )
            self.gro_file = self.topology_manager.get_structure(
                self.build_file, self.topology_label, self.path, filetype="gro"
            )
            self.sdf_file = self.topology_manager.get_structure(
                self.build_file, self.topology_label, self.path, filetype="sdf"
            )
            if self.ff_method == "bespoke":
                print("Bespoke parameters detected, re-generating force-field offxml...")
                if self.topology_manager.check_file_type(self.build_file, self.topology_label, "offxml"):
                    self.topology_manager.get_force_field(self.build_file, self.topology_label, self.path)
                else:
                    self.assign_parameters()
        else:
            self.assign_parameters()

    def set_simulation_engine(self, md_engine_object):
        self.md_engine = md_engine_object

    def assign_parameters(self):
        # Need to pass on existing topology manager to keep track of added strucutres/ff files
        top_file, gro_file, sdf_file = self.ff_generator.assign_parameters()
        self.top_file = top_file
        self.gro_file = gro_file
        self.sdf_file = sdf_file

        self.topology_manager.add_topology(
            self.top_file, self.build_file, self.topology_label
        )
        self.topology_manager.add_structure(
            self.gro_file, self.build_file, self.topology_label
        )
        self.topology_manager.add_structure(
            self.sdf_file, self.build_file, self.topology_label
        )

    def minimize(self):
        self.md_engine.center_configuration(
            self.gro_file, self.gro_file.split(".gro")[0] + "_box.gro"
        )
        self.gro_file = self.gro_file.split(".gro") + "_box.gro"
        self.md_engine.minimize(self.gro_file, self.top_file)


class SystemTopologyGenerator:
    def __init__(
        self,
        system_pdb,
        system_molecule_files,
        system_charge_files,
        output_file,
        ff_method,
        build_file,
        path="",
        ff_names=["openff-2.0.0"],
        topology_manager=TopologyManager(),
    ):
        self.build_file = build_file
        self.label = "system"
        self.path = path
        self.name = output_file
        self.topology_manager = topology_manager
        if not os.path.isdir(self.path):
            make_path(path)

        self._ff_generation_methods = {
            "openff": SystemOFFDefault,
            "bespoke" : SystemOFFDefault,
        }

        if ff_method in self._ff_generation_methods.keys():
            self.ff_generator = self._ff_generation_methods[ff_method](
                system_molecule_files,
                system_charge_files,
                system_pdb, self.name,
                path=self.path,
                force_field_strings = ff_names,
                ff_id=ff_method,
                topology_manager = self.topology_manager
            )
        else:
            print(
                "WARNING: "
                + ff_method
                + " is not one of the available "
                + "force field parameter generation methods. Please pick from:\n"
                + " ".join(self._ff_generation_methods.keys())
            )
            sys.exit()

        # Define other attributes populated by other functions
        self.md_engine = None
        self.top_file = None
        self.gro_file = None

    def set_simulation_engine(self, md_engine_object):
        self.md_engine = md_engine_object

    def get_parameters(self):
        if self.topology_manager.check_file_type(self.build_file, self.label, "top"):
            print(
                "Getting Topology Files from",
                self.topology_manager.topology_object + "...",
            )
            self.top_file = self.topology_manager.get_topology(
                self.build_file, self.label, self.path
            )
            self.gro_file = self.topology_manager.get_structure(
                self.build_file, self.label, self.path, filetype="gro"
            )
        else:
            self.assign_parameters() 
            self.minimize()
            self.topology_manager.add_structure(self.gro_file, self.build_file, label = self.label)
            self.topology_manager.add_topology(self.top_file, self.build_file, label = self.label)

    def assign_parameters(self):
        top_file, gro_file = self.ff_generator.assign_parameters()
        self.top_file = top_file
        self.gro_file = gro_file

    def minimize(self):
        self.md_engine.center_configuration(
            self.gro_file, self.gro_file.split(".gro")[0] + "_box.gro"
        )
        self.gro_file = self.gro_file.split(".gro")[0] + "_box.gro"
        self.md_engine.minimize(self.gro_file, self.top_file, prefix = "em_" + self.label)
        self.gro_file = "em_" + self.label + ".gro"
