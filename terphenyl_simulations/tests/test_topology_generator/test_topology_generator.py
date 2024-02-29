"""
Testcases for force_fileds.py objects

These objects are used to generate force field topology files using the
output from the FoldamerBuilder/SystemBuilder objects
"""

from terphenyl_simulations.force_fields import FoldamerOFFDefault, FoldamerOFFBespoke
from terphenyl_simulations.build import SimulationTopologyGenerator
from terphenyl_simulations.utils import ROOT_DIR, make_path
import pytest
import os
import shutil

@pytest.fixture
def setup_default_ff_tests():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_topology_generator"))
    if os.path.isdir("output"):
        make_path("output")
    yield FoldamerOFFDefault("mop_dimer.mol", "mop_dimer.pdb", path="output")
    os.chdir(top_dir)

@pytest.fixture
def setup_default_ff_tests_post_charges():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_topology_generator"))
    if os.path.isdir("output"):
        make_path("output")
    top_generator = FoldamerOFFDefault("mop_dimer.mol", "mop_dimer.pdb", path="output")
    shutil.copy('mop_dimer_charges.sdf', 'output/mop_dimer_charges.sdf')
    yield top_generator
    os.chdir(top_dir)


def test_default_topology_generator(setup_default_ff_tests):
    default_ff_object = setup_default_ff_tests
    assert os.path.exists('output/mop_dimer_renum.pdb')
    assert default_ff_object.omm_topology.getNumAtoms() == 109
    assert default_ff_object.omm_topology.getNumBonds() == 114
    assert default_ff_object.omm_topology.getNumChains() == 1

def test_default_tg_charges(setup_default_ff_tests):
    default_ff_object = setup_default_ff_tests
    default_ff_object._get_partial_charges()
    assert os.path.exists('output/mop_dimer_charges.sdf')

def test_default_tg_charges_from_file(setup_default_ff_tests_post_charges):
    default_ff_object = setup_default_ff_tests_post_charges
    default_ff_object._get_partial_charges()
    assert default_ff_object.molecule

def test_default_tg_output(setup_default_ff_tests_post_charges):
    default_ff_object = setup_default_ff_tests_post_charges
    default_ff_object._get_partial_charges()
    default_ff_object._generate_ff_topologies()
    assert os.path.exists('output/mop_dimer_openff-2.0.0.top')
    assert os.path.exists('output/mop_dimer_openff-2.0.0.gro')


