"""
Testcases for build.py objects

Many of these objects should be able to generate simulation input files
using just the input foldamer.build file
"""

from terphenyl_simulations.build import FoldamerBuilder, SystemBuilder
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


@pytest.fixture
def setup_system_builder_tests():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_build_workflow"))
    if os.path.isdir("output"):
        shutil.rmtree("output")
    fb = FoldamerBuilder("mop_tetramer.build", path="output")
    fb.build_foldamer()
    fb.write_pdb()
    yield SystemBuilder("mop_tetramer.build", path="output")
    os.chdir(top_dir)


def test_foldamer_builder_chain(setup_foldamer_builder_tests):
    builder = setup_foldamer_builder_tests
    builder.build_foldamer()
    assert builder.chain.n_particles == 197
    for label in builder.chain.labels["monomer"]:
        assert label.name == "MOP"
    for label in builder.chain.labels["Compound"]:
        assert label.name == "CAP"


def test_foldamer_builder_file_writing(setup_foldamer_builder_tests):
    builder = setup_foldamer_builder_tests
    builder.build_foldamer()
    builder.write_pdb()
    builder.write_mol()
    assert os.path.exists("output/mop_tetramer.pdb")
    assert os.path.exists("output/mop_tetramer.mol")


def test_system_builder_inp(setup_system_builder_tests):
    builder = setup_system_builder_tests
    builder.build_packmol_inp()
    assert os.path.exists("output/solvate.inp")
    with open("output/solvate.inp", "r") as f:
        for line in f.readlines():
            assert "OUTPUT_FILENAME" not in line
            assert "SOLUTE_PDB" not in line
            assert "SOLUTE_POSITION" not in line
            assert "SOLVENT_PDB" not in line
            assert "N_SOLVENT" not in line
            assert "SOLVENT_BOX" not in line


def test_system_builder_packmol(setup_system_builder_tests):
    system_builder = setup_system_builder_tests
    system_builder.build_packmol_inp()
    system_builder.solvate_system()
    assert os.path.exists("output/solvated_mop_tetramer.pdb")
    assert os.path.exists("output/TCM.pdb")
