# terphenyl_simulations

Repository for storing simulation input files and keeping track of changes to simulation inputs/parameters over time. This repo will primarily be used for running and analyzing Gromacs simulations on CU Boulder's Summit supercomputer. If additional submission scripts are needed, be sure to specify which super computer the submission script is written for (i.e. submit\_remd\_sim.bridges.slurm)

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
