#!/usr/bin/env python
"""
Heteropolymer Simulations
A package for running terphenyl oligomery systems
"""

import sys
import os
from importlib_metadata import entry_points
from setuptools import setup, find_packages
import versioneer
short_description = __doc__.split("\n")


# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:]),




setup(
    # Self-descriptive entries which should always be present
    name='heteropolymer_simulations',
    author='Theodore Fobe',
    author_email='theodore.fobe@colorado.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    # install_requires=['numpy', 'os', 'sys'],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

    # Entry points
    entry_points={
        'console_scripts': [
            'renumber_pdb_atoms = heteropolymer_simulations.scripts:renumber_pdb_atoms',
            'top_to_itp = heteropolymer_simulations.scripts:top_to_itp',
            'plot_edr_observable = heteropolymer_simulations.scripts:plot_edr_observables',
            'hmr_topology = heteropolymer_simulations.scripts:hmr_topology',
            'parameterize_foldamer = heteropolymer_simulations.scripts:parameterize_foldamer',
            'average_rtt = heteropolymer_simulations.scripts:calculate_average_rtt',
            'REMD_setup = heteropolymer_simulations.scripts:REMD_setup',
        ]
    }

)
