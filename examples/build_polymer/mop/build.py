import mbuild as mb
from mbuild.lib.recipes.polymer import Polymer

comp = mb.load('C(=O)c1ccc(c2c(cc(OC)cc2)c2cccc(c2)[C@@H](N)C)cc1', smiles = True)
chain = Polymer()

chain.add_monomer(compound=comp,
                  indices=[1, -2],
                  separation=.15,
                  replace=True)

chain.build(n=6, sequence='A')

chain.visualize(show_ports=True)