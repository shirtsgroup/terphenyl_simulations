from .utils import renumber_pdb_atoms, make_path
from openmm import app
from openff.toolkit.topology import FrozenMolecule, Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange.components.interchange import Interchange
import pdb
import os


class FoldamerOFFDefault:
    def __init__(self, mol_file, pdb_file, path="", ff_str="openff-2.0.0"):
        if not os.path.isdir(path):
            make_path(path)
        self.molecule = Molecule.from_file(mol_file)
        self.name = mol_file.split("/")[-1].split(".mol")[0]
        self.path = path
        pdb_path = pdb_file.split(".pdb")[0]
        renumber_pdb_atoms(pdb_file, os.path.join(path, pdb_path + "_renum.pdb"))
        self.pdb_file = app.PDBFile(os.path.join(path, pdb_path + "_renum.pdb"))
        self.omm_topology = self.pdb_file.topology
        self.off_topology = Topology.from_openmm(self.omm_topology, unique_molecules=[self.molecule])
        self.force_field = ForceField(ff_str + ".offxml")

    def assign_parameters(self, charge_method = 'am1bcc'):
        self._get_partial_charges(charge_method)
        top_file, gro_file = self._generate_ff_topologies()
        return top_file, gro_file

    def _get_partial_charges(self, method="am1bcc"):
        sdf_file = os.path.join(
            self.path, self.name + "_charges.sdf"
        )
        print(sdf_file)
        if not os.path.exists(sdf_file):
            # Expensive step
            self.molecule.assign_partial_charges(partial_charge_method=method)
            self.molecule.to_file(sdf_file, file_format="sdf")
        else:
            self.molecule = Molecule.from_file(sdf_file)

    def _generate_ff_topologies(self):
        interchange = Interchange.from_smirnoff(
            force_field=self.force_field,
            topology=self.off_topology,
            charge_from_molecules=[self.molecule],
        )
        interchange.positions = self.pdb_file.getPositions()

        top_file = os.path.join(self.path, self.name + "_openff-2.0.0.top")
        gro_file = os.path.join(self.path, self.name + "_openff-2.0.0.gro")
        interchange.to_top(top_file)
        interchange.to_gro(gro_file)

        return top_file, gro_file



class FoldamerOFFBespoke:
    def __init__(self, mol_file, ff_str="openff-2.0"):
        pass
