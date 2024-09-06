import unittest
from unittest.mock import patch, mock_open, MagicMock
import os
import pickle
import shutil
import uuid
import yaml
import json
import pytest
from glob import glob
from terphenyl_simulations.build import (
    TopologyManager, ROOT_DIR
)  # Assuming this class is in a module named 'topology_manager'
import terphenyl_simulations


class TestTopologyManager(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def navigate_to_test_dir(self):
        top_dir = os.path.abspath("")
        os.chdir(os.path.join(ROOT_DIR, "tests/test_topology_manager"))
        yield
        shutil.rmtree("test_dir")
        os.chdir(top_dir)

    @patch("os.path.isdir", return_value=False)
    @patch("os.path.exists", return_value=False)
    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.dump")
    def test_init_new_topology_manager(
        self, mock_pickle_dump, mock_open, mock_exists, mock_isdir
    ):
        # Test initialization when the topology manager does not exist
        with patch("builtins.print") as mock_print:
            tm = TopologyManager(
                topology_dir="test_dir", topology_object="test.pkl"
            )
            mock_isdir.assert_called_once_with("test_dir")
            mock_open.assert_called_once_with("test_dir/test.pkl", "wb")
            mock_pickle_dump.assert_called_once()
            mock_print.assert_any_call(
                "Creating a new TopologyManager:", "test_dir/test.pkl"
            )

    @patch("os.path.exists", return_value=True)
    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.load")
    def test_init_load_existing_topology_manager(
        self, mock_pickle_load, mock_open, mock_exists
    ):
        # Test initialization when the topology manager already exists
        with patch("builtins.print") as mock_print:
            mock_pickle_load.return_value = {"topology_dictionary": {}}
            tm = TopologyManager(
                topology_dir="test_dir", topology_object="test.pkl"
            )
            mock_open.assert_called_once_with("test_dir/test.pkl", "rb")
            mock_pickle_load.assert_called_once()
            mock_print.assert_any_call(
                "Loading TopologyManager:", "test_dir/test.pkl"
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load", return_value={"mock_key": "mock_value"})
    @patch("pickle.dump")
    @patch("json.dumps", return_value='{"mock_key": "mock_value"}')
    def test_get_build_json(
        self, mock_json_dumps, mock_pickle_dump, mock_yaml_load, mock_open
    ):
        # Test converting YAML to JSON
        tm = TopologyManager(
                topology_dir="test_dir", topology_object="test.pkl"
            )
        build_file = "mock_build.yaml"
        result = tm.get_build_json(build_file)
        mock_open.assert_called_with(build_file, "r")
        mock_yaml_load.assert_called_once()
        mock_json_dumps.assert_called_once()
        mock_pickle_dump.assert_called()
        self.assertEqual(result, '{"mock_key": "mock_value"}')

    @patch("uuid.uuid5", return_value=uuid.UUID("12345678-1234-5678-1234-567812345678"))
    @patch("json.dumps", return_value='{"mock_key": "mock_value"}')
    @patch("yaml.safe_load", return_value={"mock_key": "mock_value"})
    @patch("pickle.dump")
    @patch("builtins.open", new_callable=mock_open)
    def test_get_entry_dir_id(
        self, mock_open, mock_pickle_dump, mock_yaml_load, mock_json_dumps, mock_uuid5
    ):
        # Test getting the unique directory ID
        tm = TopologyManager(topology_dir="test_dir", topology_object="test.pkl")
        build_file = "mock_build.yaml"
        result = tm.get_entry_dir_id(build_file)
        expected = "12345678123456781234567812345678"
        self.assertEqual(result, expected)
        mock_pickle_dump.assert_called()
        mock_yaml_load.assert_called()

    @patch("glob.glob", return_value=["test_file.pdb"])
    @patch(
        "terphenyl_simulations.build.TopologyManager.get_entry_dir_id",
        return_value="mock_dir_id",
    )
    @patch("shutil.copy")
    def test_check_file_type(
        self, mock_copy, mock_get_entry_dir_id, mock_glob
    ):
        # Test checking file type in directory
        tm = TopologyManager(topology_dir="test_dir", topology_object="test.pkl")
        tm.add_structure("test_file.pdb", "mock_build.yml", "label")
        result = tm.check_file_type("mock_build.yaml", "label", "pdb")
        mock_get_entry_dir_id.assert_called_with("mock_build.yaml")
        mock_copy.assert_called_once()
        self.assertTrue(result)

    @patch(
        "terphenyl_simulations.build.TopologyManager.get_entry_dir_id",
        return_value="mock_dir_id",
    )
    @patch("terphenyl_simulations.build.TopologyManager.save")
    def test_add_buildfile_entry(self, mock_save, mock_get_entry_dir_id):
        # Test adding a buildfile entry
        tm = TopologyManager(topology_dir="test_dir", topology_object="test.pkl")
        tm.add_buildfile_entry("mock_build.yaml")
        mock_get_entry_dir_id.assert_called_with("mock_build.yaml")
        self.assertIn("mock_dir_id", tm.topology_dictionary)
        mock_save.assert_called()

    @patch("shutil.copy")
    @patch("os.path.exists", return_value=False)
    @patch("terphenyl_simulations.build.TopologyManager.save")
    def test_add_structure(self, mock_save, mock_exists, mock_copy):
        # Test adding structure files
        tm = TopologyManager(topology_dir="test_dir", topology_object="test.pkl")
        tm.topology_dictionary = {"mock_dir_id": {"label": {"structure_files": []}}}
        with patch(
            "terphenyl_simulations.build.TopologyManager.get_entry_dir_id",
            return_value="mock_dir_id",
        ):
            tm.add_structure("mock_structure.pdb", "mock_build.yaml", "label")
            mock_copy.assert_called_once_with(
                "mock_structure.pdb", "test_dir/mock_dir_id/label/mock_structure.pdb"
            )
            self.assertIn(
                "test_dir/mock_dir_id/label/mock_structure.pdb",
                tm.topology_dictionary["mock_dir_id"]["label"]["structure_files"],
            )
            mock_save.assert_called()

    @patch("shutil.copy")
    @patch("os.path.exists", return_value=False)
    @patch(
        "terphenyl_simulations.build.TopologyManager.get_entry_dir_id",
        return_value="mock_dir_id",
    )
    def test_get_structure(self, mock_get_entry_dir_id, mock_exists, mock_copy):
        # Test retrieving the structure file
        tm = TopologyManager(topology_dir="test_dir")
        tm.topology_dictionary = {
            "mock_dir_id": {
                "label": {
                    "structure_files": ["test_dir/mock_dir_id/label/mock_structure.pdb"]
                }
            }
        }
        result = tm.get_structure("mock_build.yaml", "label", "output_dir")
        mock_copy.assert_called_once_with(
            "test_dir/mock_dir_id/label/mock_structure.pdb",
            "output_dir/mock_structure.pdb",
        )
        self.assertEqual(result, "output_dir/mock_structure.pdb")


if __name__ == "__main__":
    unittest.main()
