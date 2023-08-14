"""
Heteropolymer Simulations
A package for running terphenyl oligomery systems
"""

import terphenyl_simulations.utils
import terphenyl_simulations.remd_utils
import terphenyl_simulations.observables
import terphenyl_simulations.scripts
import terphenyl_simulations.clustering
import terphenyl_simulations.plotting
import terphenyl_simulations.edit_conf

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
