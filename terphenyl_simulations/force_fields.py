from .utils import renumber_pdb_atoms, make_path
from abc import ABC, abstractclassmethod
from openmm import app
from openff.toolkit.topology import FrozenMolecule, Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange.components.interchange import Interchange
import pdb
import os
import sys


class OFFMethod(ABC):
    @abstractclassmethod
    def assign_parameters(self):
        pass


class FoldamerOFFDefault(OFFMethod):
    def __init__(self, mol_file, pdb_file, output_file = None, path="", ff_str="openff-2.0.0"):
        if type(mol_file) is not str:
            print('FoldamerOFFDefault method does cannot take process multiple molecules, try' + \
                  'SystemOFFDefault method instead.')
            sys.exit()
        if not os.path.isdir(path):
            make_path(path)
        self.molecule = Molecule.from_file(mol_file)
        self.name = output_file
        if output_file is None:
            self.name = mol_file.split("/")[-1].split(".mol")[0]
        self.path = path
        pdb_path = pdb_file.split(".pdb")[0]
        renumber_pdb_atoms(pdb_file, os.path.join(path, pdb_path + "_renum.pdb"))
        self.pdb_file = app.PDBFile(os.path.join(path, pdb_path + "_renum.pdb"))
        self.omm_topology = self.pdb_file.topology
        self.off_topology = Topology.from_openmm(
            self.omm_topology, unique_molecules=[self.molecule]
        )
        self.force_field = ForceField(ff_str + ".offxml")

    def assign_parameters(self, charge_method="am1bcc"):
        self._get_partial_charges(charge_method)
        top_file, gro_file = self._generate_ff_topologies()
        return top_file, gro_file

    def _get_partial_charges(self, method="am1bcc"):
        sdf_file = os.path.join(self.path, self.name + "_charges.sdf")
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


class FoldamerOFFBespoke(OFFMethod):
    def __init__(self, mol_file, ff_str="openff-2.0"):
        pass


class SystemOFFDefault(OFFMethod):
    def __init__(
        self, system_molecules_list, system_pdb, path="", ff_str="openff-2.0.0"
    ):
        if not os.path.isdir(path):
            make_path(path)
        self.path = path
        self.molecules = []
        self.charges = []
        for molecule in system_molecules_list:
            off_molecule = Molecule.from_file(molecule)
            self.molecules.append(off_molecule)
            self.charges.append(off_molecule.partial_charges != None)
        self.pdb_file = system_pdb
        self.off_topology = Topology.from_pdb(
            self.pdb_file, unique_molecules=self.molecules
        )
        self.force_field = ForceField(ff_str + ".offxml")

    def _generate_ff_topology(self):
        interchange = Interchange.from_smirnoff(
            force_fild = self.force_field,
            topology = self.off_topology,
            charge_from_molecules = [mol for mol, charge in zip(self.molecules, self.charges) if charge]
        )
        interchange.positions = self.pdb_file.getPositions()

        top_file = os.path.join(self.path, self.name + "_openff-2.0.0.top")
        gro_file = os.path.join(self.path, self.name + "_openff-2.0.0.gro")
        interchange.to_top(top_file)
        interchange.to_gro(gro_file)

        return top_file, gro_file
    
    def assign_parameters(self):
        top_file, gro_file = self._generate_ff_topologies()
        return top_file, gro_file