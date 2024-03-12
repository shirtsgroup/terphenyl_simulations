"""
Test cases for GromacsWrapper object
"""

import os
import pytest
import shutil
from terphenyl_simulations.utils import ROOT_DIR, make_path
from terphenyl_simulations.gromacs_wrapper import GromacsWrapper


@pytest.fixture
def setup_default_gmx_wrapper():
    # Navigate to specific test directory
    top_dir = os.path.abspath("")
    os.chdir(os.path.join(ROOT_DIR, "tests/test_gromacs_wrapper"))
    if os.path.isdir("output"):
        shutil.rmtree("output")
    yield GromacsWrapper(path="output")
    os.chdir(top_dir)


def test_gromacs_wrapper_object(setup_default_gmx_wrapper):
    gmx_wrapper = setup_default_gmx_wrapper
    assert os.path.isdir("output")


def test_gromacs_wrapper_minimize(setup_default_gmx_wrapper):
    gmx_wrapper = setup_default_gmx_wrapper
    if gmx_wrapper.gmx is None:
        pytest.skip("Unable to find GMX executable")
    gmx_wrapper.minimize("mop_dimer.gro", "mop_dimer.top")
    assert os.path.exists("output/em.tpr")
    assert os.path.exists("output/em.trr")
    assert os.path.exists("output/em.gro")
    assert os.path.exists("output/em.edr")
