import yaml
import os
import mbuild as mb
from mbuild.lib.recipes.polymer import Polymer


class FoldamerBuilder:
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
    
    def _write_to_file(self, path = ""):
        filename = os.path.join(path, self.build_params['structure_file'] + ".pdb")
        print(filename)
        self.chain.save(filename, overwrite = True, residues = [self.build_params['residue_name'], 'CAP'])

class SystemBuilder:
    pass

class TopologyFactory:
    pass