from .utils import renumber_pdb_atoms
from openmm import app
from openff.toolkit.topology import FrozenMolecule, Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange.components.interchange import Interchange
import pdb
import os

class OpenForceFieldDefaultFF:
    def __init__(self, mol_file, pdb_file, prefix = ''):
        self.molecule = Molecule.from_file(mol_file)
        self.name = mol_file.split('/')[-1].split('.mol')[0]
        self.prefix = prefix
        pdb_prefix = pdb_file.split('.pdb')[0]
        renumber_pdb_atoms(pdb_file, os.join(prefix, pdb_prefix + '_renum.pdb'))
        self.pdb_file = app.PDBFile(os.join(prefix, pdb_prefix + '_renum.pdb'))
        self.omm_topology = self.pdb_file.topology
        self.off_topology = Topology.from_openmm(self.omm_topology, unique_molecules = [])
        self.force_field = ForceField('openff-2.0.0.offxml')

    def get_partial_charges(self, method = 'am1bcc'):
        sdf_file = os.path.join(self.prefix, self.molecule.name.split('.mol')[0] + '_charges.sdf')
        if not os.path.exists(sdf_file):
            self.molecule.assign_partial_charges(partial_charge_method = method)
            self.molecule.to_file(sdf_file, file_format = 'sdf')
        else:
            self.molecule.from_file(sdf_file)
    
    def generate_ff_topologies(self):
        interchange = Interchange.from_smirnoff(
            forcefield = self.force_field,
            topology = self.off_topology,
            charge_from_molecule = [self.molecule]
        )
        interchange.positions = self.pdb_file.getPositions()

        interchange.to_top(os.path.join(self.prefix, self.name + '_openff-2.0.0.top'))
        interchange.to_gro(os.path.join(self.prefix, self.name + '_openff-2.0.0.gro'))


class BespokeFitFF:
    def __init__(self, mol_file, ff_str = 'openff-2.0'):
        pass
