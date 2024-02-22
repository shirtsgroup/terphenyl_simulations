from abc import ABC, abstractclassmethod


class MDEngineWrapper:
    @abstractclassmethod
    def minimize(self):
        pass


class GromacsWrapper(MDEngineWrapper):
    """
    Gromacs wrapper for small operations that can be run without submitting
    jobs to a cluster.
    """

    def __init__(self):
        pass
