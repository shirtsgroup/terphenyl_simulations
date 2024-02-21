import yaml
import os
import sys
import shutil
from subprocess import Popen, PIPE
import mbuild as mb
import warnings
from mbuild.lib.recipes.polymer import Polymer
from .utils import ROOT_DIR, replace_all_pattern

PACKMOL = shutil.which('packmol')

class FoldamerBuilder:
    """
    The FoldamerBuilder object is used to build foldamers with MBuild based on
    instructions provided in a foldamer build YML file. These YML input files
    provide the monomer unit SMILE string, connection atoms, and capping residues.
    It may take some hand-tuning of the Polymer object's parameters to get the
    foldamer be chemically correct. These parameters can be modified in the input 
    build.yml file.
    """
    def __init__(self, build_file_yml):
        with open(build_file_yml, 'r') as f:
            self.build_params = yaml.safe_load(f)
    
    def build_foldamer(self, write = True, path = ""):
        smile_str = self.build_params['foldamer_smile']
        subunit = mb.load(smile_str, smiles = True, name = self.build_params['residue_name'])
        head_cap = mb.load(self.build_params['cap_smiles']['upper'], smiles = True)
        tail_cap = mb.load(self.build_params['cap_smiles']['lower'], smiles = True)
        
        # We use the MBuild Polymer object to build our foldamer model
        self.chain = Polymer()
        self.chain.add_monomer(compound=subunit,
                  indices=[self.build_params['upper_connect'], self.build_params['lower_connect']],
                  separation=.15,
                  replace=True,
                  # orientation = [[0,-1,0],[1,0,0]]
                 )
        self.chain.add_end_groups(compound = head_cap,
                            index = -1,
                            separation=0.15,
                            label="head",
                            duplicate = False
                            )

        self.chain.add_end_groups(compound = tail_cap,
                            index = -1,
                            separation=0.15,
                            label="tail",
                            duplicate = False
                            )
        self.chain.build(n=self.build_params['foldamer_length'], sequence='A')
        
        # Change residue names in chain object
        for label in self.chain.labels['monomer']:
            label.name = self.build_params['residue_name']
        for label in self.chain.labels["Compound"]:
            label.name = "CAP"

        if write is True:
            self._write_to_file(path = path)
    
    def _write_to_file(self, path = ''):
        filename = os.path.join(path, self.build_params['structure_file'] + '.pdb')
        print(filename)
        self.chain.save(filename, overwrite = True, residues = [self.build_params['residue_name'], 'CAP'])

# Requires packmol
class SystemBuilder:
    """
    """
    def __init__(self, build_file_yml, path = ''):
        with open(build_file_yml, 'r') as f:
            self.build_params = yaml.safe_load(f)
        self.packmol = PACKMOL
        self.path = path

    def build_packmol_inp(self):
        self.output_file =  os.path.join(self.path, "solvate.inp")

        # Get template solvation script from terphenyl_simulations/data/solvents
        shutil.copy(
            os.path.join(ROOT_DIR, 'data/solvents/packmol_solvate_template.inp'),
            self.output_file
        )

        # Make replacements to the template inp file
        replace_all_pattern('OUTPUT_FILENAME', os.path.join(self.path, 'solvated_' + self.build_params['structure_file'] + '.pdb'), self.output_file)
        replace_all_pattern('SOLUTE_PDB', os.path.join(self.path, self.build_params['structure_file'] + '.pdb'), self.output_file)
        
        # Create string for solute position
        solute_position = [
            self.build_params['system']['box_size'] / 2,
            self.build_params['system']['box_size'] / 2, 
            self.build_params['system']['box_size'] / 2, 
            45,
            45,
            45
        ]
        solute_position_str = [str(round(n, 1)) for n in solute_position]
        replace_all_pattern('SOLUTE_POSITION', ' '.join(solute_position_str), self.output_file)
        
        # Get solvent PDB location
        solvent_pdb = self.build_params['system']['solvent'] + '.pdb'
        print(os.listdir(os.path.join(ROOT_DIR, 'data/solvents/')))
        if not os.path.exists(solvent_pdb) and solvent_pdb in os.listdir(os.path.join(ROOT_DIR, 'data/solvents/')):
            print('Using ' + solvent_pdb + ' from the internal library of solvents. If this ' + \
                  'is not the behavior you intend, please include the solvent pdb you wish ' + \
                  'to solvate the system with in your working directory.')
            solvent_pdb = os.path.join(ROOT_DIR, 'data/solvents/', self.build_params['system']['solvent'] + '.pdb')
        else:
            warnings.warn('Warning! Unable to find ' + solvent_pdb + ' in the working directory \
                          or the internal library.')
        replace_all_pattern('SOLVENT_PDB', solvent_pdb, self.output_file)
        replace_all_pattern('N_SOLVENT', str(self.build_params['system']['n_solvent']), self.output_file)
        
        # Create string for simulation box specification
        box_parameters = [
            0,
            0,
            0,
            self.build_params['system']['box_size'],
            self.build_params['system']['box_size'], 
            self.build_params['system']['box_size'],
        ]
        box_parameters_str = [str(round(n, 1)) for n in box_parameters]
        
        replace_all_pattern('SOLVENT_BOX', ' '.join(box_parameters_str), self.output_file)

    def build_system(self):
        cmdline_entry = self.packmol + ' < ' + self.output_file
        with open('temp', 'w') as f:
            f.write(cmdline_entry)

        process = Popen(['bash temp'], shell = True, stdin=PIPE)
        process.wait()
        os.remove('temp')

class TopologyFactory:
    pass