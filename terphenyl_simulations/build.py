import yaml
import os
import sys
import shutil
import pickle
import json
import uuid
from glob import glob
from subprocess import Popen, PIPE
import mbuild as mb
import warnings
from openbabel import openbabel
from abc import ABC, abstractclassmethod
from mbuild.lib.recipes.polymer import Polymer
from .utils import ROOT_DIR, replace_all_pattern, make_path, renumber_pdb_atoms
from .force_fields import FoldamerOFFDefault, FoldamerOFFBespoke, SystemOFFDefault

PACKMOL = shutil.which("packmol")

class TopologyManager:
    """
    The TopologyManager class is responsible for storing and keeping track of
    existing molecule topologies. Since parameter assigment can be a costly
    computation, this object stores topology files in a local file system, 
    and provides these files to simulations when needed.
    """

    def __init__(self, topology_dir = os.path.join(ROOT_DIR, "data/topology_files"), topology_object = os.path.join(ROOT_DIR, "data/topology_files/top_manager.pkl")):
        # Save values to object
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
            print("Loading TopologyManager:", self.topology_object)
            # If a pickled object exists for this object
            # Load it into this object
            self.load()

    def save(self):
        with open(self.topology_object, 'wb') as fw:
            pickle.dump(self.__dict__, fw)

    def load(self):
        with open(self.topology_object, 'rb') as fr:
            tmp_dict = pickle.load(fr)
        self.__dict__.update(tmp_dict)

    def get_build_json(self, build_file):
        with open(build_file, "r") as stream:
            topology_dict = yaml.safe_load(stream)

        topology_json = json.dumps(topology_dict)
        return topology_json

    def get_entry_dir_id(self, top_json):
        unique_dir_str = str(uuid.uuid5(uuid.NAMESPACE_X500, str(top_json)))
        unique_dir = "".join([a for a in unique_dir_str if a != "-"])
        return unique_dir

    def check_file_type(self, build_file, label, filetype):
        topology_json = self.get_build_json(build_file)
        unique_dir = self.get_entry_dir_id(topology_json)
        files = glob(os.path.join(self.topology_dir, unique_dir, label, "*." + filetype))
        return len(files) > 0

    def add_buildfile_entry(self, build_file):
        print("Adding", build_file, "to TopologyManager...")
        topology_json = self.get_build_json(build_file)
        
        # Generate entry for internal dictionary
        self.topology_dictionary[topology_json] = {}

        # Make filesystem reflecting new entry
        # Generate unique key for JSON file

        unique_dir = self.get_entry_dir_id(topology_json)
        output_dir = os.path.join(self.topology_dir, unique_dir)
        os.makedirs(output_dir)
        self.save()

    def add_entry_label(self, build_file, label):
        # Add label to entry
        topology_json = self.get_build_json(build_file)
        self.topology_dictionary[topology_json][label] = {
            "structure_files" :[],
            "topology_files" : [],
        }

        # Add directory to local storage
        unique_dir = self.get_entry_dir_id(topology_json)
        label_dir = output_dir = os.path.join(self.topology_dir, unique_dir, label)
        os.makedirs(label_dir)
        self.save()
        
    def add_structure(self, structure_file, build_file, label):
        topology_json = self.get_build_json(build_file)
        if topology_json not in self.topology_dictionary.keys():
            self.add_buildfile_entry(build_file)
        if label not in self.topology_dictionary[topology_json].keys():
            self.add_entry_label(build_file, label)

        # Save files internally
        entry_id = self.get_entry_dir_id(topology_json)
        label_directory = os.path.join(self.topology_dir, entry_id, label)
        if os.path.exists(os.path.join(label_directory, structure_file)):
            return
        else:
            shutil.copy(structure_file, label_directory)

        # Add to internal dictionary
        # But check if extension already exists
        self.topology_dictionary[topology_json][label]["structure_files"].append(os.path.join(label_directory, structure_file))

        self.save()

    def get_structure(self, build_file, label, filetype = None):
        topology_json = self.get_build_json(build_file)
        if filetype is None:
            return self.topology_dictionary[topology_json][label]["structure_files"][0]
        else:
            stored_file_types = [structure.split(".")[-1] for structure in self.topology_dictionary[topology_json][label]["structure_files"]]
            index = stored_file_types.index(filetype)
            return self.topology_dictionary[topology_json][label]["structure_files"][index]

    def add_topology(self, topology_file, build_file, label):
        topology_json = self.get_build_json(build_file)

        # Add to internal dictionary
        self.topology_dictionary[topology_json][label]["topology_file"] = topology_file

        # Save files internally
        entry_id = self.get_entry_dir_id(topology_json)
        label_directory = os.path.join(self.topology_dir, entry_id, label)
        shutil.copy(topology_file, label_directory)
        self.save()

    def get_topology(self):
        pass

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

    def __init__(self, build_file_yml, path="", topology_manager = TopologyManager(), topology_label = "molecule"):
        self.build_file = build_file_yml
        self.label = topology_label
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
            self.chain = mb.load(self.topology_manager.get_structure(self.build_file, self.label))
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


# Requires packmol
class SystemBuilder:
    """
    """
    def __init__(self, solute_pdb, build_file_yml, path=""):
        with open(build_file_yml, "r") as f:
            self.build_params = yaml.safe_load(f)
        self.packmol = PACKMOL
        self.path = path
        if ".pdb" not in solute_pdb:
            solute_pdb += ".pdb"
        self.solute_pdb = solute_pdb
        if not os.path.isdir(self.path):
            make_path(path)

    def build_packmol_inp(self):
        self.inp_file = os.path.join(self.path, "solvate.inp")

        # Get template solvation script from terphenyl_simulations/data/solvents
        shutil.copy(
            os.path.join(ROOT_DIR, "data/solvents/packmol_solvate_template.inp"),
            self.inp_file,
        )

        # Get solvent PDB location

        # Check top signac directory
        solvent_pdb = os.path.join(
            self.path, self.build_params["system"]["solvent"] + ".pdb"
        )

        # Check stored solvent pdbs

        if not os.path.exists(solvent_pdb):
            if solvent_pdb.split("/")[-1] in os.listdir(
            os.path.join(ROOT_DIR, "data/solvents/")):
                print(
                    "Using "
                    + solvent_pdb
                    + " from the internal library of solvents. If this "
                    + "is not the behavior you intend, please include the solvent pdb you wish "
                    + "to solvate the system with in your working directory."
                )
                stored_solvent = os.path.join(
                    ROOT_DIR,
                    "data/solvents/",
                    self.build_params["system"]["solvent"] + ".pdb",
                )

                # Copy to path directory
                shutil.copy(
                    stored_solvent, os.path.join(self.path, solvent_pdb.split("/")[-1])
                )
            else:
                warnings.warn(
                    "Warning! Unable to find "
                    + solvent_pdb
                    + " in the working directory "
                    + "or the internal library. "
                )
                sys.exit()

        # Make replacements to the template inp file
        replace_all_pattern(
            "OUTPUT_FILENAME",
            os.path.join(
                self.path, "solvated_" + self.build_params["structure_file"] + ".pdb"
            ),
            self.inp_file,
        )

        solute_pdb_str = os.path.join(self.path, self.solute_pdb)

        replace_all_pattern(
            "SOLUTE_PDB",
            solute_pdb_str,
            self.inp_file,
        )

        # Create string for solute position
        solute_position = [
            self.build_params["system"]["box_size"] / 2,
            self.build_params["system"]["box_size"] / 2,
            self.build_params["system"]["box_size"] / 2,
            45,
            45,
            45,
        ]
        solute_position_str = [str(round(n, 1)) for n in solute_position]
        replace_all_pattern(
            "SOLUTE_POSITION", " ".join(solute_position_str), self.inp_file
        )
        replace_all_pattern("SOLVENT_PDB", solvent_pdb, self.inp_file)
        replace_all_pattern(
            "N_SOLVENT", str(self.build_params["system"]["n_solvent"]), self.inp_file
        )

        # Create string for simulation box specification
        box_parameters = [
            0,
            0,
            0,
            self.build_params["system"]["box_size"],
            self.build_params["system"]["box_size"],
            self.build_params["system"]["box_size"],
        ]
        box_parameters_str = [str(round(n, 1)) for n in box_parameters]

        replace_all_pattern("SOLVENT_BOX", " ".join(box_parameters_str), self.inp_file)

    def solvate_system(self):
        # Since we are chdir into the output path, we need to remove
        cmdline_entry = self.packmol + " < " + self.inp_file
        with open("temp", "w") as f:
            f.write(cmdline_entry)

        print(cmdline_entry)
        process = Popen(["bash temp"], shell=True, stdin=PIPE)
        process.wait()
        os.remove("temp")

# I probably best to get rid of this class
class MoleculeTopologyGenerator:
    def __init__(self, molecule_file, pdb_file, build_file, ff_method, path="", ff_name="openff-2.0.0", topology_manager = TopologyManager(), topology_label = "molecule"):
        self.path = path
        self.build_file = build_file
        self.topology_manager = topology_manager
        self.topology_label = topology_label
        if not os.path.isdir(self.path):
            make_path(path)

        self._ff_generation_methods = {
            "openff-foldamer" : FoldamerOFFDefault,
            "bespoke-foldamer" : FoldamerOFFBespoke,
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
        if self.topology_manager.check_file_type(self.build_file, self.topology_label, "top"):
            self.top_file = self.topology_manager.get_topology(self.build_file, self.topology_label)
            self.gro_file = self.topology_manager.get_structure(self.build_file, self.topology_label)
        else:
            self.assign_parameters()

    def set_simulation_engine(self, md_engine_object):
        self.md_engine = md_engine_object

    def assign_parameters(self):
        top_file, gro_file = self.ff_generator.assign_parameters()
        self.top_file = top_file
        self.gro_file = gro_file
        self.topology_manager.add_topology(self.top_file, self.build_file, self.topology_label)

    def minimize(self):
        self.md_engine.center_configuration(
            self.gro_file, self.gro_file.split(".gro") + "_box.gro"
        )
        self.gro_file = self.gro_file.split(".gro") + "_box.gro"
        self.md_engine.minimize(self.gro_file, self.top_file)

class SystemTopologyGenerator:
    def __init__(self, molecule_files, charge_files, pdb_file, output_file, ff_method, path="", ff_name="openff-2.0.0"):
        self.path = path
        self.name = output_file
        if not os.path.isdir(self.path):
            make_path(path)

        self._ff_generation_methods = {
            "openff-system" : SystemOFFDefault,
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

    def assign_parameters(self):
        top_file, gro_file = self.ff_generator.assign_parameters()
        self.top_file = top_file
        self.gro_file = gro_file

    def minimize(self):
        self.md_engine.center_configuration(
            self.gro_file, self.gro_file.split(".gro") + "_box.gro"
        )
        self.gro_file = self.gro_file.split(".gro") + "_box.gro"
        self.md_engine.minimize(self.gro_file, self.top_file)