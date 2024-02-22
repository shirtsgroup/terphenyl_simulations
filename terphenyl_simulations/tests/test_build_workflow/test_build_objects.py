"""
Testcases for build.py objects

Many of these objects should be able to generate simulation input files
using just the input foldamer.build file
"""

from terphenyl_simulations.build import FoldamerBuilder, SystemBuilder
import pytest
import shutil
import os

@pytest.fixture(autouse=True)
def setup_test_directory():
    if os.path.isdir('output'):
        shutil.rmtree('output')


def test_foldamer_builder_chain():
    builder = FoldamerBuilder('mop_tetramer.build', path = 'output')
    builder.build_foldamer()
    assert builder.chain.n_particles == 197
    for label in builder.chain.labels['monomer']:
        assert label.name == 'MOP'
    for label in builder.chain.labels['Compound']:
        assert label.name == 'CAP'
    

def test_foldamer_builder_file_writing():
    builder = FoldamerBuilder('mop_tetramer.build', path = 'output')
    builder.build_foldamer()
    builder.write_pdb()
    builder.write_mol()
    assert os.path.exists('output/mop_tetramer.pdb')
    assert os.path.exists('output/mop_tetramer.mol')


def test_system_builder_inp():
    builder = SystemBuilder('mop_tetramer.build', path = 'output')
    builder.build_packmol_inp()
    assert os.path.exists('output/solvate.inp')
    with open('output/solvate.inp', 'r') as f:
        for line in f.readlines():
            assert 'OUTPUT_FILENAME' not in line
            assert 'SOLUTE_PDB' not in line
            assert 'SOLUTE_POSITION' not in line
            assert 'SOLVENT_PDB' not in line
            assert 'N_SOLVENT' not in line
            assert 'SOLVENT_BOX' not in line

def test_system_builder_packmol():
    builder = FoldamerBuilder('mop_tetramer.build', path = 'output')
    builder.build_foldamer()
    builder.write_pdb()
    system_builder = SystemBuilder('mop_tetramer.build', path = 'output')
    system_builder.build_packmol_inp()
    system_builder.solvate_system()
    assert os.path.exists('output/solvated_mop_tetramer.pdb')
    assert os.path.exists('output/TCM.pdb')
