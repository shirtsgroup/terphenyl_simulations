from abc import ABC, abstractclassmethod
import os
import shutil
import subprocess
from .utils import ROOT_DIR, make_path

class MDEngineWrapper:
    @abstractclassmethod
    def minimize(self):
        pass


class GromacsWrapper(MDEngineWrapper):
    """
    Gromacs wrapper for small operations that can be run without submitting
    jobs to a cluster.
    """

    def __init__(self, gromacs_exe = None, path = '.'):
        self.path = path
        if not os.path.exists(path):
            make_path(path)
        if gromacs_exe is not None and os.path.exists(gromacs_exe):
            self.gmx = gromacs_exe
        else:
            self.gmx = shutil.which('gmx')

        self.mpi = False
        if 'mpi' in self.gmx:
            self.mpi = True
            self.gmx = 'mpirun -np 1 ' + self.gmx

    def center_configuration(self, gro_file, out_file):
        editconf_call = self.gmx + ' editconf -f ' + gro_file + \
                                            ' -c yes' + \
                                            ' -o ' + out_file
        process = subprocess.Popen(editconf_call.split(' '))
        process.wait()

    def add_box(self, gro_file, box_size, out_file):
        editconf_call = self.gmx + ' editconf -f ' + gro_file + \
                                            ' -box ' + str(box_size) +  \
                                            ' -o ' + out_file
        process = subprocess.Popen(editconf_call.split(' '))
        process.wait()


    def minimize(self, gro_file, top_file, prefix = 'em', mdp = 'em.mdp'):
        # Check for mdp file
        if not mdp in os.listdir(self.path):
            mdp = os.path.join(ROOT_DIR, 'data/simulation_templates/remd/em.mdp')

        # Run gmx grompp
        grompp_call = self.gmx + ' grompp -f ' + mdp + \
                                        ' -c ' + gro_file + \
                                        ' -p ' + top_file + \
                                        ' -o ' + os.path.join(self.path, prefix)
        
        process = subprocess.Popen(grompp_call.split(' '))
        process.wait()

        # Run gmx mdrun
        mdrun_call = self.gmx + ' mdrun -deffnm ' + os.path.join(self.path, prefix)

        process = subprocess.Popen(mdrun_call.split(' '))
        process.wait()



        
            




        
