"""
Testcases for force_fileds.py objects

These objects are used to generate force field topology files using the
output from the FoldamerBuilder/SystemBuilder objects
"""

from terphenyl_simulations.force_fields import FoldamerOFFDefault, FoldamerOFFBespoke
from terphenyl_simulations.build import (
    MoleculeTopologyGenerator,
    SystemTopologyGenerator,
    TopologyManager,
)
from terphenyl_simulations.utils import ROOT_DIR, make_path
import pytest
import os
import shutil


@pytest.fixture
def navigate_to_test_dir():
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_topology_generator"))
    yield
    os.chdir(top_dir)


@pytest.fixture
def setup_default_ff_tests():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_topology_generator"))
    if os.path.isdir("output"):
        shutil.rmtree("output")
    yield FoldamerOFFDefault("mop_dimer.mol", "mop_dimer.pdb", path="output")
    shutil.rmtree("output")
    os.chdir(top_dir)


@pytest.fixture
def setup_default_ff_tests_post_charges():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_topology_generator"))
    if os.path.isdir("output"):
        make_path("output")
    top_generator = FoldamerOFFDefault("mop_dimer.mol", "mop_dimer.pdb", path="output")
    shutil.copy("mop_dimer_charges.sdf", "output/mop_dimer_charges.sdf")
    yield top_generator
    shutil.rmtree("output")
    os.chdir(top_dir)


def test_default_topology_generator(setup_default_ff_tests):
    default_ff_object = setup_default_ff_tests
    assert os.path.exists("output/mop_dimer_renum.pdb")
    assert default_ff_object.omm_topology.getNumAtoms() == 109
    assert default_ff_object.omm_topology.getNumBonds() == 114
    assert default_ff_object.omm_topology.getNumChains() == 1


@pytest.mark.skip(reason="this test is very slow unless OpenEye is installed")
def test_default_tg_charges(setup_default_ff_tests):
    default_ff_object = setup_default_ff_tests
    default_ff_object._get_partial_charges()
    assert os.path.exists("output/mop_dimer_charges.sdf")


def test_default_tg_charges_from_file(setup_default_ff_tests_post_charges):
    default_ff_object = setup_default_ff_tests_post_charges
    default_ff_object._get_partial_charges()
    assert default_ff_object.molecule


def test_default_tg_output(setup_default_ff_tests_post_charges):
    default_ff_object = setup_default_ff_tests_post_charges
    default_ff_object._get_partial_charges()
    default_ff_object._generate_ff_topologies()
    assert os.path.exists("output/mop_dimer_openff-2.0.0.top")
    assert os.path.exists("output/mop_dimer_openff-2.0.0.gro")


def test_topology_manager_save_load(navigate_to_test_dir):
    assert navigate_to_test_dir is None
    tp_manager = TopologyManager(topology_dir="", topology_object="top_manager.pkl")
    tp_manager.topology_dictionary["test_key"] = "test_value"
    tp_manager.save()

    tp_manager_load = TopologyManager(
        topology_dir="", topology_object="top_manager.pkl"
    )

    assert os.path.exists("top_manager.pkl")
    assert "test_key" in tp_manager_load.topology_dictionary.keys()
    assert tp_manager_load.topology_dictionary["test_key"] == "test_value"
    os.remove("top_manager.pkl")


def test_topology_manager_save_top_file(setup_default_ff_tests_post_charges):
    default_ff_object = setup_default_ff_tests_post_charges
    default_ff_object._get_partial_charges()
    default_ff_object._generate_ff_topologies()

    tp_manager = TopologyManager(
        topology_dir="output_2", topology_object="top_manager.pkl"
    )
    tp_manager.add_topology(
        "mop_tetramer.build",
        [
            "output/mop_dimer_charges.sdf",
            "output/mop_dimer_openff-2.0.0.gro",
            "output/mop_dimer_openff-2.0.0.top",
            "output/mop_dimer_renum.pdb",
        ],
    )

    assert os.path.exists("output_2/top_manager.pkl")
    assert os.path.exists("output_2/mop_dimer_charges.sdf")
    assert os.path.exists("output_2/mop_dimer_openff-2.0.0.gro")
    assert os.path.exists("output_2/mop_dimer_openff-2.0.0.top")
    assert os.path.exists("output_2/mop_dimer_renum.pdb")

    shutil.rmtree("output_2")


def test_topology_manager_save_top_file(setup_default_ff_tests_post_charges):
    default_ff_object = setup_default_ff_tests_post_charges
    default_ff_object._get_partial_charges()
    default_ff_object._generate_ff_topologies()

    tp_manager = TopologyManager(
        topology_dir="output_2", topology_object="top_manager.pkl"
    )
    tp_manager.add_topology(
        "mop_tetramer.build",
        [
            "output/mop_dimer_charges.sdf",
            "output/mop_dimer_openff-2.0.0.gro",
            "output/mop_dimer_openff-2.0.0.top",
            "output/mop_dimer_renum.pdb",
        ],
    )

    loaded_tp = TopologyManager(
        topology_dir="output_2", topology_object="top_manager.pkl"
    )
    print(loaded_tp)

    assert loaded_tp.topology_dictionary.keys() == tp_manager.topology_dictionary.keys()
    shutil.rmtree("output_2")
