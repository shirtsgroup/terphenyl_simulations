#!/usr/bin/env python

try:
    from openmm import app
except ImportError:
    from simtk.openmm import app

from openff.toolkit.topology import FrozenMolecule, Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange.components.interchange import Interchange
import pdb
import os
import argparse
import panedr
import pandas as pd
import heteropolymer_simulations as hs
import parmed as pmd
import matplotlib.pyplot as plt


def renumber_pdb_atoms():
    def parse_args():
        parser = argparse.ArgumentParser(
            description="A quick script to rename atoms in a pdb \
                        file. This can be used to give unique atom \
                        to all atoms in the pdb file."
        )

        parser.add_argument(
            "-f", "--file", type=str, help="file name of original file to change"
        )

        parser.add_argument(
            "-o", "--output", type=str, help="output file name to write new pdb file"
        )

        return parser.parse_args()

    args = parse_args()
    hs.utils.renumber_pdb_atoms(args.file, args.output)


def top_to_itp():
    def parse_args():
        parser = argparse.ArgumentParser(
            description="A quick script to convert system .top \
                        files output from the openff-toolkit to  \
                        molecule specific .itp files."
        )

        parser.add_argument(
            "-t", "--top", type=str, help="file name of original .top file to convert"
        )

        parser.add_argument(
            "-o", "--output", type=str, help="output file name to write new .itp file"
        )

        return parser.parse_args()

    args = parse_args()
    top = hs.utils.TopFileObject(args.top)
    hs.utils.write_itp_file(top, args.output)


def plot_edr_observables():
    def parse_args():
        parser = argparse.ArgumentParser(
            description="A script to plot observables from \
                        gromacs simulations. Supply this function \
                        with which .edr obsevables you want plotted."
        )

        parser.add_argument(
            "--obs",
            type=str,
            nargs="+",
            required=True,
            help="List of observables to extract from .edr files provide \
                    possible entries include: 'Time', 'Bond', 'Angle', \
                    'Proper Dih.', 'LJ-14', 'Coulomb-14', 'LJ (SR)', \
                    'Coulomb (SR)', 'COM Pull En.', 'Potential', \
                    'Kinetic En.', 'Total Energy', 'Temperature', \
                    'Pressure', 'Constr. rmsd', 'Vir-XX', 'Vir-XY', \
                    'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', \
                    'Vir-ZY', ', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', \
                    'Vir-ZY', 'Pres-YZ', 'Pres-ZX', 'Pres-ZY', 'Pres-ZZ', \
                    '#Surf*SurfTen' and 'T-System'",
        )

        parser.add_argument(
            "--units",
            type=str,
            nargs="+",
            required=True,
            help="Units to be plotted with observables. If there are spaces \
                    in your units input, use quotation marks if spaces are in \
                    your input.",
        )

        parser.add_argument(
            "--time_series",
            action="store_true",
            help="flag to plot time series of provided data",
        )

        parser.add_argument(
            "--hist",
            action="store_true",
            help="flag to plot histogram of provided data",
        )

        parser.add_argument(
            "--base_file_name",
            help="basename of .edr, .trr, .log, files to search for in each directory",
            nargs="+",
        )

        parser.add_argument(
            "--directories",
            help="directories to search for specified file names. Example sim{0..63}",
            nargs="+",
        )

        parser.add_argument(
            "-o",
            "--output_base",
            type=str,
            help="extension to add to output figure names",
            default="",
        )

        return parser.parse_args()

    args = parse_args()
    edf_list = []
    for sim_dir in args.directories:
        all_files = os.listdir(sim_dir)
        for file_id in args.base_file_name:
            edr_files = [fn for fn in all_files if ".edr" in fn]
            edr_files = [fn for fn in edr_files if file_id in fn]
            edr_files.sort()
            energy_df = None
            for edr_file in edr_files:
                print("Loading:", os.path.join(sim_dir, edr_file))
                if energy_df is None:
                    energy_df = panedr.edr_to_df(os.path.join(sim_dir, edr_file))
                else:
                    energy_df = pd.concat(
                        [energy_df, panedr.edr_to_df(os.path.join(sim_dir, edr_file))]
                    )
            energy_df.sort_values(by="Time")
        edf_list.append(energy_df)
    for i, column_name in enumerate(args.obs):
        if column_name in edf_list[0].columns:
            if args.time_series:
                figure_name = column_name.replace(" ", "_").lower()
                figure_name = figure_name.replace("(", "")
                figure_name = figure_name.replace(")", "")
                figure_name = figure_name + "_" + args.output_base + ".png"
                plt.figure(dpi=300)
                for edf in edf_list:
                    time = edf["Time"] * 0.001
                    observable = edf[column_name]
                    plt.plot(time, observable)
                plt.xlabel("Time (ns)")
                if args.units:
                    plt.ylabel(column_name + " " + args.units[i])
                else:
                    plt.ylabel(column_name)
                plt.savefig(figure_name)
            if args.hist:
                figure_name = column_name.replace(" ", "_").lower()
                figure_name = figure_name.replace("(", "")
                figure_name = figure_name.replace(")", "")
                figure_name = figure_name + "_hist_" + args.output_base + ".png"
                plt.figure(dpi=300)
                for edf in edf_list:
                    observable = edf[column_name]
                    plt.hist(observable, bins=10, density=True)
                plt.ylabel("Density")
                if args.units:
                    plt.xlabel(column_name + " " + args.units[i])
                else:
                    plt.xlabel(column_name)
                plt.savefig(figure_name)


def hmr_topology():
    """
    Quick script to apply HMR to a GROMACS topology file
    """

    def parse_args():
        parser = argparse.ArgumentParser(
            description="A quick script to run HMR on  .top \
                        files"
        )

        parser.add_argument(
            "-t", "--top", type=str, help="file name of original .top file to convert"
        )

        parser.add_argument(
            "-o", "--output", type=str, help="output file name to write new .itp file"
        )

        parser.add_argument(
            "--hmr_ratio",
            type=float,
            help="ratio of the the HMR hydrogen mass to the original mass.\
                An hmr_ratio of 3 will scale the hydrogen masses by 3.",
        )

        return parser.parse_args()

    args = parse_args()
    gmx_top = pmd.load_file(args.top)

    for atom in gmx_top.atoms:
        if atom.element_name == "H":
            new_mass = atom.mass * args.hmr_ratio
            d_mass = new_mass - atom.mass
            for b_atom in atom.bond_partners:
                b_atom.mass -= d_mass
            atom.mass = new_mass

    gmx_top.write(args.output)

def parameterize_foldamer():

    def parse_args():
        parser = argparse.ArgumentParser(
            description="A script to assign parameters using OpenFF workflow"
        )

        parser.add_argument(
            "--mol", type=str, help="MOL file of the polymer to be parameterized"
        )

        parser.add_argument(
            "--pdb", type=str, help="PDB structure file of polymer to be parameterized"
        )

        parser.add_argument(
            "--output",
            type=str,
            help="Output name of files"
        )

        parser.add_argument(
            "--sdf", type=str, help="Optional SDF file if AM1-BCC charges are already assigned"
        )

        parser.add_argument(
            "--ff", type=str, help="Force field XML file to assign parameters from"
        )

        return parser.parse_args()

    args = parse_args()

    # Read MOL and PDB files and get OMM topology
    print("Reading input structures...")
    molecule = Molecule.from_file(args.mol)
    pdbfile = app.PDBFile(args.pdb)
    omm_topology = pdbfile.topology

    # Create OpenFF topology
    off_topology = Topology.from_openmm(
        omm_topology, unique_molecules=[molecule]
    )

    # Calculate or Read in partial charges
    print("Getting partial charges...")
    if args.sdf is None:
        args.sdf = args.output + ".sdf"

    if not os.path.exists(args.sdf):
        molecule.assign_partial_charges(partial_charge_method="am1bcc")
        molecule.to_file(args.sdf, file_format = 'sdf')
    else:
        molecule = Molecule.from_file(args.sdf)

    # Load OpenFF force field
    print("Assigning FF parameters...")
    if args.ff is None:
        args.ff = "openFF-2.0.0.offxml"
    forcefield = ForceField(args.ff)

    # Prepare OpenMM system
    omm_system = forcefield.create_openmm_system(off_topology, charge_from_molecule = [molecule])

    # Create Interchange object
    interchange = Interchange.from_smirnoff(
        force_filed=forcefield,
        topology=off_topology,
        charge_from_molecules = [molecule]
    )
    interchange.positions = pdbfile.positions

    # Export to Gromacs Files
    interchange.to_top(args.output + ".top")
    interchange.to_top(args.output + ".gro")
