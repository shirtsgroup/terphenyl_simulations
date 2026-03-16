"""
Testcases for build.py objects

Many of these objects should be able to generate simulation input files
using just the input foldamer.build file
"""

from terphenyl_simulations.build import FoldamerBuilder
import terphenyl_simulations
from unittest.mock import patch, mock_open, MagicMock
from terphenyl_simulations.utils import ROOT_DIR
import pytest
import shutil
import os

# These fixtures setup and breakdown test cases
# Setup involves navigating to the correct test directory
# and removing any exisiting output


@pytest.fixture
def setup_foldamer_builder_tests():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_build_workflow"))
    if os.path.isdir("output"):
        shutil.rmtree("output")
    yield FoldamerBuilder("mop_tetramer.build", path="output")
    os.chdir(top_dir)

@patch("terphenyl_simulations.build.TopologyManager.add_structure")
def test_foldamer_builder_chain(mock_add_structure, setup_foldamer_builder_tests):
    builder = setup_foldamer_builder_tests
    builder.build_foldamer()
    assert builder.chain.n_particles == 197
    for label in builder.chain.labels["monomer"]:
        assert label.name == "MOP"
    for label in builder.chain.labels["Compound"]:
        assert label.name == "CAP"

@patch("terphenyl_simulations.build.TopologyManager.add_structure")
def test_foldamer_builder_file_writing(mock_add_stricture, setup_foldamer_builder_tests):
    builder = setup_foldamer_builder_tests
    builder.build_foldamer()
    builder.write_pdb()
    builder.write_mol()
    assert os.path.exists("output/mop_tetramer.pdb")
    assert os.path.exists("output/mop_tetramer.mol")
