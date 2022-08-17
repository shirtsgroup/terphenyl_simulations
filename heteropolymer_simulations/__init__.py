"""
Heteropolymer Simulations
A package for running terphenyl oligomery systems
"""

import heteropolymer_simulations.utils
import heteropolymer_simulations.observables
import heteropolymer_simulations.scripts
import heteropolymer_simulations.clustering
import heteropolymer_simulations.plotting

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions