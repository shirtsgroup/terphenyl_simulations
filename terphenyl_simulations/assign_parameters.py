from .utils import renumber_pdb_atoms, make_path
from abc import ABC, abstractmethod
from openmm import app
from openff.toolkit.topology import FrozenMolecule, Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange.components.interchange import Interchange
from openff.bespokefit.executor import BespokeExecutor, BespokeWorkerConfig, wait_until_complete
from openff.bespokefit.workflows import BespokeWorkflowFactory
from openff.fragmenter.fragment import WBOFragmenter
from openff.bespokefit.schema.targets import TorsionProfileTargetSchema
from openff.qcsubmit.common_structures import QCSpec
from openff.bespokefit.utilities.smirks import SMIRKSettings
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.schema.smirnoff import  ProperTorsionHyperparameters
from .topology_manager import TopologyManager
import pdb
import os
import sys
import yaml
import subprocess


class OFFMethod(ABC):
    @classmethod
    @abstractmethod
    def assign_parameters(self):
        pass


class FoldamerOFFBespoke(OFFMethod):
    def __init__(
        self, mol_file, pdb_file, build_file = None, output_file=None, path="", ff_str="openff-2.0.0"
    ):
        if type(mol_file) is not str:
            print(
                "FoldamerOFFDefault method does cannot process multiple molecules, try"
                + "SystemOFFDefault method instead."
            )
            sys.exit()
        if not os.path.isdir(path):
            make_path(path)
        self.molecule = Molecule.from_file(mol_file)
        self.name = output_file
        self.initial_ff = ff_str
        self.build_file_yml = build_file
        self.label = "molecule"
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
        self.force_field = None
        self.topology_manager = TopologyManager()

    def assign_parameters(self):
        if not self.topology_manager.check_file_type(self.build_file_yml, "molecule", "offxml"):
            self.generate_trimer_molecule(self.build_file_yml)
            self.assign_trimer_partial_charges()
            self.run_bespoke_fit_workflow()
            self._get_partial_charges()
        else:
            self.force_field = self.topology_manager.get_force_field(self.build_file_yml, "molecule", self.path)
            self.sdf_file = self.topology_manager.get_structure(self.build_file_yml, "molecule", self.path, filetype="sdf")
        top_file, gro_file = self.generate_ff_topologies()
        return top_file, gro_file, self.sdf_file

    def _get_partial_charges(self, method="am1bcc"):
        self.sdf_file = os.path.join(self.path, self.name + "_charges.sdf")
        if not os.path.exists(self.sdf_file):
            self.molecule.assign_partial_charges(partial_charge_method = method)
            self.molecule.to_file(self.sdf_file, file_format="sdf")
        else:
            self.molecule = Molecule.from_file(self.sdf_file)

    def generate_trimer_molecule(self, build_file_yml):
        from .build import FoldamerBuilder

        with open(build_file_yml, "r") as f:
            self.build_params = yaml.safe_load(f)
        self.build_params["foldamer_length"] = 3
        self.trimer_buildfile = self.build_params["structure_file"].split("_")[0] + "_trimer.build"
        self.build_params["structure_file"] = self.build_params["structure_file"].split("_")[0] + "_trimer"
        with open(self.trimer_buildfile, "w") as wf:
            yaml.dump(self.build_params, wf)
        foldamer_builder = FoldamerBuilder(self.trimer_buildfile)
        foldamer_builder.get_foldamer()

    def assign_trimer_partial_charges(self):
        with open(self.trimer_buildfile, "r") as f:
            trimer_build_params = yaml.safe_load(f)
        self.trimer_molecule = Molecule.from_file(trimer_build_params["structure_file"] + ".mol")
        pdb_file = trimer_build_params["structure_file"] + ".pdb"
        renumber_pdb_atoms(pdb_file, os.path.join(self.path, pdb_file + "_renum.pdb"))
        self.trimer_pdb  = app.PDBFile(os.path.join(self.path, pdb_file + "_renum.pdb"))

        self.topology_manager.load()
        self.trimer_sdf_file = trimer_build_params["structure_file"] + ".sdf"
        if not self.topology_manager.check_file_type(self.trimer_buildfile, "molecule", "sdf"):
            self.trimer_molecule.assign_partial_charges(partial_charge_method = "am1bcc")
            self.trimer_molecule.to_file(self.trimer_sdf_file, file_format="sdf")
            self.topology_manager.add_structure(self.trimer_sdf_file, self.trimer_buildfile, "molecule")
        else:
            self.topology_manager.get_structure(self.trimer_buildfile, "molecule", self.path, filetype="sdf")

    def run_bespoke_fit_workflow(self, 
                                   n_fragmenter_workers = 4,
                                   n_qc_compute_workers = 4,
                                   n_optimizer_workers = 4,
                            ):

        # Keep Bespoke Executor Files
        subprocess.run(["BEFLOW_KEEP_TMP_FILES=True"], shell=True)

        # Setup Executor object
        self.bespoke_fit_executor = BespokeExecutor(
            n_fragmenter_workers = n_fragmenter_workers,
            n_qc_compute_workers = n_qc_compute_workers,
            n_optimizer_workers = n_optimizer_workers,
            launch_redis_if_unavailable = True
        )

        # Setup Workflow
        self.factory = BespokeWorkflowFactory()
        self.factory.target_torsion_smirks = ['[!#1]~[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]~[!#1]']
        self.factory.fragmentation_engine = WBOFragmenter()
        self.factory.target_templates = [TorsionProfileTargetSchema()]
        self.factory.default_qc_specs = [
            QCSpec(
                method="gfn2xtb",
                basis=None,
                program="xtb",
                spec_name="xtb",
                spec_description="gfn2xtb",
            )
        ]
        self.factory.smirk_settings = SMIRKSettings(
            expand_torsion_terms=True,
            generate_bespoke_terms=True,
        )
        self.force_field = None
        self.factory.optimizer = ForceBalanceSchema()
        self.factory.parameter_hyperparameters = [ProperTorsionHyperparameters()]
        self.factory.to_file('bespoke_flow.json')

#     def run_bespoke_fit_executor(self):
        trimer_molecule = Molecule.from_file(self.trimer_sdf_file)
        bespoke_workflow_schema = self.factory.optimization_schema_from_molecule(trimer_molecule)

        with self.bespoke_fit_executor:
            task_id = self.bespoke_fit_executor.submit(bespoke_workflow_schema)
            output = wait_until_complete(task_id)
            self.force_field = output.bespoke_force_field

        ff_file_name = self.build_params["structure_file"] + "_bespoke_" + self.initial_ff + ".offxml"
        self.force_field.to_file(ff_file_name)
        self.topology_manager.add_force_field(ff_file_name, self.build_file_yml, self.label)



    def generate_ff_topologies(self):
        interchange = self.force_field.create_interchange(
            self.off_topology,
            charge_from_molecules = [self.molecule]
        )
        interchange.positions = self.pdb_file.getPositions()

        top_file = os.path.join(self.path, self.name + "_bespoke_" + self.initial_ff + ".top")
        gro_file = os.path.join(self.path, self.name + "_bespoke_" + self.initial_ff + ".gro")

        interchange.to_top(top_file)
        interchange.to_gro(gro_file)

        return top_file, gro_file

class FoldamerOFFDefault(OFFMethod):
    def __init__(
        self, mol_file, pdb_file, output_file=None, path="", ff_str="openff-2.0.0", build_file = None
    ):
        if type(mol_file) is not str:
            print(
                "FoldamerOFFDefault method does cannot process multiple molecules, try"
                + "SystemOFFDefault method instead."
            )
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
        return top_file, gro_file, self.sdf_file

    def _get_partial_charges(self, method="am1bcc"):
        self.sdf_file = os.path.join(self.path, self.name + "_charges.sdf")
        if not os.path.exists(self.sdf_file):
            # Expensive step
            self.molecule.assign_partial_charges(partial_charge_method=method)
            self.molecule.to_file(self.sdf_file, file_format="sdf")
        else:
            self.molecule = Molecule.from_file(self.sdf_file)

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


class SystemOFFDefault(OFFMethod):
    def __init__(
        self,
        system_molecules_list,
        charge_files,
        system_pdb,
        output_file,
        path="",
        force_field_strings=["openff-2.0.0"],
    ):
        if not os.path.isdir(path):
            make_path(path)
        self.path = path
        self.molecules = []
        self.charges = []
        self.name = output_file
        for molecule, charge_file in zip(system_molecules_list, charge_files):
            off_molecule = Molecule.from_file(molecule)
            off_molecule.name = molecule.split(".")[0]
            # If we have partial charges in a file use it
            if charge_file != None and os.path.exists(charge_file):
                print("Getting charges for", molecule, "from", charge_file)
                molecule_charges = Molecule.from_file(charge_file).partial_charges
                off_molecule.partial_charges = molecule_charges
                off_molecule.perceive_residues()
            self.molecules.append(off_molecule)
            self.charges.append(off_molecule.partial_charges != None)
        self.pdb_file = app.PDBFile(system_pdb)
        self.off_topology = Topology.from_openmm(
            self.pdb_file.topology, unique_molecules=self.molecules
        )
        
        self.force_fields = [ForceField(ff_str + ".offxml") for ff_str in force_field_strings]

    def _generate_ff_topology(self):
        
        if len(self.force_fields) == 1:
            self.force_fields *= len(self.molecules)
        
        interchange = None
        os.environ["INTERCHANGE_EXPERIMENTAL"] = "1"
        for mol in self.off_topology.molecules:
            ff_index = self.molecules.index(mol)
            ff =  self.force_fields[ff_index]
            # print("Parameterizing", mol.name, "with", ff.author,  ff.date, "release.")
            if interchange is None:
                interchange = ff.create_interchange(
                    mol.to_topology(),
                    charge_from_molecules = [mol])
            else:
                add_interchange = ff.create_interchange(mol.to_topology())
                interchange = interchange.combine(add_interchange)

        interchange.positions = self.pdb_file.getPositions()
        top_file = os.path.join(self.path, self.name + "_openff.top")
        gro_file = os.path.join(self.path, self.name + "_openff.gro")
        
        interchange.to_top(top_file)
        interchange.to_gro(gro_file)

        return top_file, gro_file

    def assign_parameters(self):
        top_file, gro_file = self._generate_ff_topology()
        return top_file, gro_file
