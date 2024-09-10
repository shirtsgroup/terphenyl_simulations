import yaml
import os
import sys
import shutil
import pickle
import json
import uuid
import glob
from subprocess import Popen, PIPE
import mbuild as mb
import warnings
from openbabel import openbabel
from abc import ABC, abstractclassmethod
from mbuild.lib.recipes.polymer import Polymer
from openff.toolkit import ForceField
from openff.toolkit.topology import Molecule, Topology
from .utils import ROOT_DIR, replace_all_pattern, make_path, renumber_pdb_atoms
from .force_fields import FoldamerOFFDefault, FoldamerOFFBespoke, SystemOFFDefault
from .gromacs_wrapper import GromacsWrapper
from openff.interchange import Interchange
from openff.interchange.components._packmol import pack_box, UNIT_CUBE
from openff.toolkit import unit


class TopologyManager:
    """
    The TopologyManager class is responsible for storing and keeping track of
    existing molecule topologies. Since parameter assigment can be a costly
    computation, this object stores topology files in a local file system,
    and provides these files to simulations when needed.
    """

    def __init__(
        self,
        topology_dir=None,
        topology_object="top_manager.pkl"
    ):
        if topology_dir is None:
            topology_dir = os.path.join(ROOT_DIR, "data/topology_manager")
        self.topology_dir = topology_dir
        self.topology_object = os.path.join(topology_dir, topology_object)

        # Generate path for saving topologies
        if not os.path.isdir(topology_dir) and len(topology_dir) > 0:
            os.makedirs(topology_dir)
        if not os.path.exists(self.topology_object):
            print("Creating a new TopologyManager:", self.topology_object)
            # Initialize object once
            # This dictionary will hold the paths to topology files
            self.topology_dictionary = {}
            self.save()
        else:
            # print("Loading TopologyManager:", self.topology_object)
            # If a pickled object exists for this object
            # Load it into this object
            self.load()
            # Need to overwrite old topology dirs on new machines
            # So path to ROOT_DIR is correct
            self.topology_dir = topology_dir
            self.topology_object = os.path.join(topology_dir, topology_object)

    def save(self):
        with open(self.topology_object, "wb") as fw:
            pickle.dump(self.__dict__, fw)

    def load(self):
        with open(self.topology_object, "rb") as fr:
            tmp_dict = pickle.load(fr)
        self.__dict__.update(tmp_dict)

    def get_build_json(self, build_file):
        with open(build_file, "r") as stream:
            topology_dict = yaml.safe_load(stream)

        topology_json = json.dumps(topology_dict)
        return topology_json

    def get_entry_dir_id(self, build_file):
        build_json = self.get_build_json(build_file)
        unique_dir_str = str(uuid.uuid5(uuid.NAMESPACE_X500, str(build_json)))
        unique_dir = "".join([a for a in unique_dir_str if a != "-"])
        return unique_dir

    def check_file_type(self, build_file, label, filetype):
        build_file_id = self.get_entry_dir_id(build_file)
        files = glob.glob(
            os.path.join(self.topology_dir, build_file_id, label, "*." + filetype)
        )
        return len(files) > 0

    def add_buildfile_entry(self, build_file):
        print("Adding", build_file, "to TopologyManager...")
        build_file_id = self.get_entry_dir_id(build_file)

        # Generate entry for internal dictionary
        self.topology_dictionary[build_file_id] = {}

        # Make filesystem reflecting new entry
        # Generate unique key for JSON file

        output_dir = os.path.join(self.topology_dir, build_file_id)
        os.makedirs(output_dir)
        self.save()

    def add_entry_label(self, build_file, label):
        # Add label to entry
        build_file_key = self.get_entry_dir_id(build_file)

        self.topology_dictionary[build_file_key][label] = {
            "structure_files": [],
            "topology_files": [],
        }

        # Add directory to local storage
        label_dir = os.path.join(self.topology_dir, build_file_key, label)
        os.makedirs(label_dir)
        self.save()

    def add_structure(self, structure_file, build_file, label):
        build_file_id = self.get_entry_dir_id(build_file)
        filename = structure_file.split("/")[-1]
        if build_file_id not in self.topology_dictionary.keys():
            self.add_buildfile_entry(build_file)
        if label not in self.topology_dictionary[build_file_id].keys():
            self.add_entry_label(build_file, label)

        # Save files internally
        label_directory = os.path.join(self.topology_dir, build_file_id, label)
        if os.path.exists(os.path.join(label_directory, filename)):
            return
        else:
            shutil.copy(structure_file, os.path.join(label_directory, filename))

        # Add to internal dictionary
        # But check if extension already exists
        self.topology_dictionary[build_file_id][label]["structure_files"].append(
            structure_file
        )
        self.save()

    def get_structure(self, build_file, label, path, filetype=None):
        build_file_id = self.get_entry_dir_id(build_file)
        # Get most recently stored file
        if filetype is None:
            database_file = self.topology_dictionary[build_file_id][label][
                "structure_files"
            ][-1]
        else:
            stored_file_types = [
                structure.split(".")[-1]
                for structure in self.topology_dictionary[build_file_id][label][
                    "structure_files"
                ]
            ]
            index = stored_file_types.index(filetype)
            database_file = self.topology_dictionary[build_file_id][label][
                "structure_files"
            ][index]
        database_file = os.path.join(self.topology_dir, build_file_id, label, database_file)
        output_file = os.path.join(path, database_file.split("/")[-1])
        shutil.copy(database_file, output_file)
        return output_file

    def add_topology(self, topology_file, build_file, label):
        build_file_id = self.get_entry_dir_id(build_file)
        # Save files internally
        label_directory = os.path.join(self.topology_dir, build_file_id, label)
        stored_filename = os.path.join(label_directory, topology_file.split("/")[-1])
        shutil.copy(topology_file, stored_filename)

        # Add to internal dictionary
        if not stored_filename.split("/")[-1] in self.topology_dictionary[build_file_id][label]["topology_files"]:
            self.topology_dictionary[build_file_id][label]["topology_files"].append(
                stored_filename.split("/")[-1]
            )
        self.save()

    def get_topology(self, build_file, label, path, filetype=None):
        build_file_id = self.get_entry_dir_id(build_file)
        if filetype is None:
            database_file = self.topology_dictionary[build_file_id][label][
                "topology_files"
            ][0]
        else:
            stored_file_types = [
                structure.split(".")[-1]
                for structure in self.topology_dictionary[build_file_id][label][
                    "structure_files"
                ]
            ]
            index = stored_file_types.index(filetype)
            database_file = self.topology_dictionary[build_file_id][label][
                "topology_files"
            ][index]

        database_file = os.path.join(self.topology_dir, build_file_id, label, database_file)
        output_file = os.path.join(path, database_file.split("/")[-1])
        shutil.copy(database_file, output_file)
        return output_file


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
        subunit = mb.load(
            smile_str, smiles=True, name=self.build_params["residue_name"]
        )
        head_cap = mb.load(self.build_params["cap_smiles"]["upper"], smiles=True)
        tail_cap = mb.load(self.build_params["cap_smiles"]["lower"], smiles=True)

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
        self.chain.name = "MOP"
        # self.chain.energy_minimize()

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
                molecule_file, pdb_file, path=self.path, ff_str=ff_name
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
