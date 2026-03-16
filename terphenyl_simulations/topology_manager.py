import os
import yaml
import shutil
import pickle
import json
import uuid
import glob
from terphenyl_simulations.utils import ROOT_DIR

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
        if self.topology_dir != os.path.join(ROOT_DIR, "data/topology_manager"):
            self.topology_dir = os.path.join(ROOT_DIR, "data/topology_manager")

        topology_entries = os.listdir(self.topology_dir)
        topology_keys = list(self.topology_dictionary.keys())
        # Remove topologies from internal dict
        # That are not present in filesystem
        update = False
        for key in topology_keys:
            if key not in topology_entries:
                print("Removing", key, "from TopologyManager...")
                update = True
                del self.topology_dictionary[key]
                continue
            # Check if labels are present too
            label_entries = os.listdir(os.path.join(self.topology_dir, key))
            for label in list(self.topology_dictionary[key].keys()):
                if label not in label_entries:
                    print("Removing label", label, "from", key, "TopologyManager entry...")
                    update = True
                    del self.topology_dictionary[key][label]
            # If no labels are under build_id, remove it
            if len(self.topology_dictionary[key].keys()) == 0:
                print("Removing", key, "from TopologyManager...")
                update = True
                del self.topology_dictionary[key]
        # if any physical files are missing, update internal object and save
        if update:
            self.save()

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
        # See if build_file_id exists
        if build_file_id in self.topology_dictionary.keys() and \
            os.path.isdir(os.path.join(self.topology_dir, build_file_id)):
                # See if label exists
                if label in self.topology_dictionary[build_file_id].keys() and \
                    os.path.isdir(os.path.join(self.topology_dir, build_file_id, label)):

                    # Aggregate all files from dictionary entry
                    files = []
                    for ftype in ["structure_files", "topology_files", "force_field_files"]:
                        files += self.topology_dictionary[build_file_id][label][ftype]
                    
                    # Check if files are present in filesystem
                    files = [file for file in files if os.path.isfile(os.path.join(self.topology_dir, build_file_id, label, file))]
                    return filetype in [f.split(".")[-1] for f in files] 
        # If file id or label or filetype is not present return false
        return False

    def add_buildfile_entry(self, build_file):
        print("Adding", build_file, "to TopologyManager...")
        build_file_id = self.get_entry_dir_id(build_file)

        # Generate entry for internal dictionary
        self.topology_dictionary[build_file_id] = {}

        # Make filesystem reflecting new entry
        # Generate unique key for JSON file

        output_dir = os.path.join(self.topology_dir, build_file_id)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        self.save()

    def add_entry_label(self, build_file, label):
        # Add label to entry
        print("Adding label", label, "to", build_file)
        build_file_key = self.get_entry_dir_id(build_file)

        if label not in self.topology_dictionary[build_file_key].keys():
            self.topology_dictionary[build_file_key][label] = {
                "structure_files": [],
                "topology_files": [],
                "force_field_files" : [],
            }

            # Add directory to local storage
            label_dir = os.path.join(self.topology_dir, build_file_key, label)
            if not os.path.isdir(label_dir):
                os.makedirs(label_dir)
            self.save()
        else:
            print("Label", label, "is already defined for", build_file)

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
                    "topology_files"
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


    def add_force_field(self, ff_file, build_file, label):
        print("Addiing force-field", ff_file, "to TopologyManager...")
        build_file_id = self.get_entry_dir_id(build_file)
        # Save files internally
        label_directory = os.path.join(self.topology_dir, build_file_id, label)
        stored_filename = os.path.join(label_directory, ff_file.split("/")[-1])
        shutil.copy(ff_file, stored_filename)

        # Add to internal dictionary
        if not "force_field_files" in self.topology_dictionary[build_file_id][label].keys():
            self.topology_dictionary[build_file_id][label]["force_field_files"] = []
        if not stored_filename.split("/")[-1] in self.topology_dictionary[build_file_id][label]["force_field_files"]:
            self.topology_dictionary[build_file_id][label]["force_field_files"].append(
                stored_filename.split("/")[-1]
            )
        self.save()


    def get_force_field(self, build_file, label, path, filetype = None):
        build_file_id = self.get_entry_dir_id(build_file)
        if filetype is None:
            database_file = self.topology_dictionary[build_file_id][label][
                "force_field_files"
            ][0]
        else:
            stored_file_types = [
                structure.split(".")[-1]
                for structure in self.topology_dictionary[build_file_id][label][
                    "force_field_files"
                ]
            ]
            index = stored_file_types.index(filetype)
            database_file = self.topology_dictionary[build_file_id][label][
                "force_field_files"
            ][index]

        database_file = os.path.join(self.topology_dir, build_file_id, label, database_file)
        output_file = os.path.join(path, database_file.split("/")[-1])
        shutil.copy(database_file, output_file)
        return output_file