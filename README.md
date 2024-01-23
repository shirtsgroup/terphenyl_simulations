# Terphenyl Simulations

[![Python package](https://github.com/shirtsgroup/terphenyl_simulations/actions/workflows/python-package.yml/badge.svg)](https://github.com/shirtsgroup/terphenyl_simulations/actions/workflows/python-package.yml)

This repository stores simulation setup and analysis for various terphenyl foldamers. Using simulations we can identify potential foldamer chemistries to help guide experimental synthesis of new foldamers. In this repository there are input files for standard MD, temperature replica exchange MD, metadynamics MD, multiple walker metadynamics MD and bias-exchange metadynamics MD simulations. We currently have simulation parameters for POP, MOP, POM, MOM and PMP terphenyl foldamers.

*This repository is still in development and is subject to rapid changes*

# Installation

`terphenyl_simulations` can be installed over the included conda environment. If you don't have anaconda on your machine, you can install it [here](https://docs.anaconda.com/free/anaconda/install/index.html). First install the `environment.yml` file with:

```
$ conda env create -f environment_new.yml
```
Then activate that conda environment with:
```
$ conda activate ts_analysis
```

Now we can install `terphenyl_simulations` with:
```
$ pip install -e .
```

This installs a develop mode version of `terphenyl_simulations`, where changes to the package are reflected in the installed version. To install this package normally, please omit the `-e` flag.

All terphenyl simulations are run using GROMACS patched PLUMED for certain enhanced sampling simulations. So ensure your machine has access to a copy of [GROMACS]{https://manual.gromacs.org/current/install-guide/index.html} and [PLUMED]{https://www.plumed.org/doc-v2.8/user-doc/html/_installation.html}. Currently this repository is setup to run simulations with GROMACS 2022.5 patched with POLUMED 2.8.2.

# Usage

## Simulations

All simulation files are in the `simulations` directory, which is organized by what chemistry is being simulated. Simulation type, foldamer length and other simulation details are specified in the directory tree (e.g. `simulations/terphenyl_mom/octamer_remd/extended/t_250_450`). 

### Standard MD Simulations

Standard MD simulations are specified by directories with the `md` tag, for example `simulations/terphenyl_mop/tetramer_md`. In these directories there are all the files needed to run a standard MD simulation of the system specified by the directory path. The GROMACS simulations can be started with:

```
$ gmx grompp -f berendsen_npt.mdp -p system.top -o berendsen_npt -s em.gro
$ gmx mdrun -deffnm berendsen_npt -v
$ gmx grompp -f berendsen_nvt.mdp -p system.top -o berendsen_nvt -s berendsen_nvt.gro
$ gmx mdrun -deffnm berendsen_nvt -v
$ gmx grompp -f npt_production.mdp -p system.top -o npt_production -s berendsen_npt.gro
$ gmx mdrun -deffnm npt_production -v
```

This block of code sets up and runs a 5 ns NVT and a 5 ns NPT simulations to equilibrate the minimized structure (`em.gro`) in solution. The production simulation is set by default to run for 100 ns, however if you wish to extend this call:

```
$ gmx convert-tpr -s npt_production.tpr -extend 100000 -o npt_new.tpr
$ gmx mdrun -deffnm npt_production -v
```

The call to `gmx convert-tpr` adds 100 ns the GROMACS topology run file, so when we call `gmx mdrun` again it will continue the simulation from the last frame of the production simulation.

Depending on what system you're working on you may need to submit the simulation to a cluster queue system. Files are provided to do this, however for individual cluster submission scripts, you will likely need to change some header information and package locations so that your computer can find instances of GROMACS and PLUMED. This repository contains `submit_all.slurm` files, which submits the simulation code above to CU Boulder's Summit Cluster.

This bash file can be run using:
```
$ bash submit_all.slurm
```

### Temperature Replica Exchange Simulations

Temperature Replica Exchange Simulations are denoted with a `remd` in the directory path. These simulations run multiple copies of a simulation at different temperatures and periodically exchange configurations to improve sampling of the configuration space of low-temperature replicas. To setup T-REMD simulations we need to provide a few parameters to `REMD_setup`, the script responsible for setting up the replica directories. For instance here we setup an T-REMD simulation with 64 replicas with simulations ranging from 250 K to 450 K:

```
$ REMD_setup -N 64 --t_range 250 450 --mdps berendsen_nvt.mdp berendsen_npt.mdp npt_production.npt --extra_files system_hmr.top en_solvated.gro
```

We also specify which MDP files and other extra files are included in each simulation directory. Here we include MDP files for the two equilibration simulations and the production simulation, as well as the system topology file and an initial structure.

Once the directories are setup for the T-REMD simulation running the following will submit equilibration simulations and the production simulation:
```
$ sbatch submit_all.slurm
```

### Metadynamics Simulations

Metadynamics simulations are used to improve sampling along specified collective variables (CVs) in longer foldamer systems. In these simulations gaussian biasing potentials are added to a specified CV to push the simulation to sample different values of the CV over the course of the simulation. Before running a metadynamics simulation, there a few important considerations that we must make before running our simulation. 

#### Unbiased CV Simulations

First, in metadynamics simulations, CV selection is an difficult and nuanced step required to have a successful metadynamics simulation. In this repository, we bias the formation of native contacts determined from unbiased T-REMD simulations. In the figure below the native contacts of this helical configruation are illustrated in pink:

<p align="center">
  <img src="https://github.com/shirtsgroup/terphenyl_simulations/blob/update_readme/figures/mop_octamer_nc.png" alt="MOP-Octamer Native Contacts"/>
</p>

The CV we use here is defined as the sum of switching functions for the distance between native H-bonds in the folded state of the terphenyl foldamer.The switching functions implemented in PLUMED range from 0 to 1, based on the distance between the oxygen atom and the hydrogen atom in the hyrdogen bond that makes up the native contact. In the figure above, since we're calculating the native contacts for an octamer, there the CV can vary from 0 to 7 possible native contacts that can form.

Before running any metadynamics simulations we must first get an understanding of how a CV behaves in the abscence of metadynamics biases to determine some initial paramters for the metadynamics simulation. An example of these unbiased simulations can be found here: `terphenyl_simulations/simulations/terphenyl_mop/tetramer_metad/h_bonds/unbiased`. Simulation files in this directory run a standard MD simulation of a MOP-terphenyl tetramer at 300 K. Similar to other simulation directories these simulations can be run using:

```
$ bash submit_all.slurm
```

**After adjusting submission scripts to work correctly on the cluster you're working on.**

After the unbiased simulation has finished running we can analyze the CVs of interest using `plumed_hbonds.dat`. This file defines the atom's used in the native contact CV. We can use `plumed driver` to calculate the timeseries of the CV after the simulation has run by running:

```
$ plumed driver --plumed plumed_hbonds.dat --pdb build_system/solvated_tetramer.pdb --mf_xtc npt_new.xtc
```

This command will generate a plain text file, `HBOND_SUMS`, which has information 


#### Multiple Walkers

#### Bias Exchange Metadynamics

## Examples

