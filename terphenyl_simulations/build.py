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
from openff.toolkit import ForceField
from openff.toolkit.topology import Molecule, Topology
from openff.interchange import Interchange
from openff.interchange.components._packmol import pack_box, UNIT_CUBE
from openff.units import unit
from .assign_parameters import FoldamerOFFDefault, FoldamerOFFBespoke, SystemOFFDefault
from .gromacs_wrapper import GromacsWrapper
from .topology_manager import TopologyManager
from .utils import ROOT_DIR, replace_all_pattern, make_path, renumber_pdb_atoms


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
        if self.topology_manager.check_file_type(self.build_file, self.label, "pdb"):
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
                    # This corrects dihedrals to be trans
                    self.chain.rotate_dihedral(bond[0:2], np.pi + (50 * np.pi / 180))


        self.chain.energy_minimize(forcefield="MMFF94")

        # Change residue names in chain object
        for label in self.chain.labels["monomer"]:
            label.name = self.build_params["residue_name"]
        for label in self.chain.labels["Compound"]:
            label.name = "CAP"

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
    """ """

    def __init__(
        self,
        solute_sdf,
        solvent_smiles,
        build_file_yml,
        path="",
        force_field="openff-2.0.0.offxml",
        topology_manager=TopologyManager(),
    ):
        self.build_file = build_file_yml
        self.force_field = ForceField(force_field)
        self.topology_manager = topology_manager
        self.topology_label = "system"
        self.path = path
        with open(build_file_yml, "r") as f:
            self.build_params = yaml.safe_load(f)

        self.solute = Molecule.from_file(solute_sdf)
        self.solvent = Molecule.from_smiles(solvent_smiles)
        self.md_engine = GromacsWrapper()

    def get_system_topology(self):
        if self.topology_manager.check_file_type(
            self.build_file, self.topology_label, "top"
        ):
            self.top_file = self.topology_manager.get_topology(
                self.build_file, self.topology_label, self.path
            )
            self.gro_file = self.topology_manager.get_structure(
                self.build_file, self.topology_label, self.path, filetype="gro"
            )
        else:
            self.build_system_topology()

    def build_system_topology(self):
        self.topology = pack_box(
            molecules=[self.solute, self.solvent],
            number_of_copies=[1, self.build_params["system"]["n_solvent"]],
            box_vectors=self.build_params["system"]["box_size"]
            * UNIT_CUBE
            * unit.angstrom,
        )

        interchange = Interchange.from_smirnoff(
            force_field=self.force_field,
            topology=self.topology,
            charge_from_molecules=[self.solute],
        )
        self.gro_file = self.build_params["structure_file"] + "_system.gro"
        self.top_file = self.build_params["structure_file"] + "_system.top"
        interchange.to_gro(self.gro_file)
        interchange.to_top(self.top_file)
        self.topology_manager.add_structure(
            self.gro_file, self.build_file, self.topology_label
        )
        self.topology_manager.add_topology(
            self.top_file, self.build_file, self.topology_label
        )

    def minimize_system(self):
        centered_gro = self.gro_file.split(".gro")[0] + "_box.gro"
        self.md_engine.center_configuration(
            self.gro_file,
            centered_gro,
        )
        self.md_engine.minimize(centered_gro, self.top_file, prefix="em_solvated")
        self.gro_file = "em_solvated.gro"


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
        if not os.path.isdir(self.path):
            make_path(path)

        self._ff_generation_methods = {
            "openff-foldamer": FoldamerOFFDefault,
            "bespoke-foldamer": FoldamerOFFBespoke,
        }

        if ff_method in self._ff_generation_methods.keys():
            self.ff_generator = self._ff_generation_methods[ff_method](
                molecule_file, pdb_file, path=self.path, ff_str=ff_name, build_file = self.build_file
            )
        else:
            warnings.warn(
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
        else:
            self.assign_parameters()

    def set_simulation_engine(self, md_engine_object):
        self.md_engine = md_engine_object

    def assign_parameters(self):
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
        molecule_files,
        charge_files,
        pdb_file,
        output_file,
        ff_method,
        path="",
        ff_name="openff-2.0.0.offxml",
        topology_manager=TopologyManager(),
    ):
        self.label = "system"
        self.path = path
        self.name = output_file
        if not os.path.isdir(self.path):
            make_path(path)

        self._ff_generation_methods = {
            "openff-system": SystemOFFDefault,
        }

        if ff_method in self._ff_generation_methods.keys():
            self.ff_generator = self._ff_generation_methods[ff_method](
                molecule_files, charge_files, pdb_file, path=self.path, ff_str=ff_name
            )
        else:
            warnings.warn(
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
                self.build_file, self.topology_label, self.path
            )
            self.gro_file = self.topology_manager.get_structure(
                self.build_file, self.topology_label, self.path
            )
        else:
            self.assign_parameters()
            self.minimize()
            self.topology_manager.add_topology(self.top_file)
            self.topology_manager.add_structure(self.gro_file)
            self.topology_manager.add_structure(self.em_gro_file)

    def assign_parameters(self):
        top_file, gro_file = self.ff_generator.assign_parameters()
        self.top_file = top_file
        self.gro_file = gro_file

    def minimize(self):
        self.md_engine.center_configuration(
            self.gro_file, self.gro_file.split(".gro") + "_box.gro"
        )
        self.gro_file = self.gro_file.split(".gro") + "_box.gro"
        self.output_file
        self.md_engine.minimize(self.gro_file, self.top_file)
        self.em_gro_file = "em_" + self.gro_file
