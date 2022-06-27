import argparse
import  heteropolymer_simulations as hs

def renumber_pdb_atoms():
    
    def parse_args():
        parser = argparse.ArgumentParser(
            description = "A quick script to rename atoms in a pdb \
                        file. This can be used to give unique atom \
                        to all atoms in the pdb file."
                        )

        parser.add_argument(
            "-f", "--file",
            type = str,
            help = "file name of original file to change"
        )

        parser.add_argument(
            "-o", "--output",
            type = str,
            help = "output file name to write new pdb file"
        )

        parser.add_argument(
            "--resname",
            type = str,
            help = "output file name to write new pdb file",
            default = "TET"
        )


        return parser.parse_args()
    
    args = parse_args()
    hs.utils.renumber_pdb_atoms(args.file, args.output, args.resname)

def top_to_itp():
    def parse_args():
        parser = argparse.ArgumentParser(
            description = "A quick script to convert system .top \
                        files output from the openff-toolkit to  \
                        molecule specific .itp files."
                        )

        parser.add_argument(
            "-t", "--top",
            type = str,
            help = "file name of original .top file to convert"
        )

        parser.add_argument(
            "-o", "--output",
            type = str,
            help = "output file name to write new .itp file"
        )

        return parser.parse_args()

    args = parse_args()
    top = hs.utils.TopFileObject(args.top)
    hs.utils.write_itp_file(top, args.output)


def plot_edr_observables():
    def parse_args():
        parser = argparse.ArgumentParser(
            description = "A script to plot observables from \
                        gromacs simulations. Supply this function \
                        with which .edr obsevables you want plotted.
                        )

        parser.add_argument(
            "--obs",
            type = str,
            nargs = "+",
            required = True,
            help = "List of observables to extract from .edr files provide \
                    possible entries include: 'Time', 'Bond', 'Angle', \
                    'Proper Dih.', 'LJ-14', 'Coulomb-14', 'LJ (SR)', \
                    'Coulomb (SR)', 'COM Pull En.', 'Potential', \
                    'Kinetic En.', 'Total Energy', 'Temperature', \
                    'Pressure', 'Constr. rmsd', 'Vir-XX', 'Vir-XY', \
                    'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', \
                    'Vir-ZY', ', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', \
                    'Vir-ZY', 'Pres-YZ', 'Pres-ZX', 'Pres-ZY', 'Pres-ZZ', \
                    '#Surf*SurfTen' and 'T-System'
        )
        
        parser.add_argument(
            "--time_series",
            action = 'store_true',
            help = "flag to plot time series of provided data"
        )

        parser.add_argument(
            "--hist",
            action = 'store_true',
            help = "flag to plot histogram of provided data"
        )

        parser.add_argument(
            "--base_file_name",
            help = "basename of .edr, .trr, .log, files to search for in each directory",
            nargs='+'
        )

        parser.add_argument(
            "--directories",
            help = "directories to search for specified file names. Example sim{0..63}",
            nargs='+'
        )

        parser.add_argument(
            "-o", "--output_base",
            type = str
            help = "extension to add to output figure names",
            default = ""
        )

        return parser.parse_args()

    args = parse_args()
    for sim_dir in args.directories:
        all_files  = os.listdir(sim_dir)
        for file_id in args.base_file_name:
            edr_files = [fn for fn in all_files if ".edr" in fn]
            edr_files.sort()
            for edr_file in edr_files:
                print("Loading:", os.path.join(sim_dir,edr_file))
                if energy_df is None:
                    energy_df = panedr.edr_to_df(os.path.join(sim_dir, edr_file))
                else:
                    energy_df = pd.concat([energy_df, panedr.edr_to_df(os.path.join(sim_dir, edr_file))])
            energy_df.sort_values(by="Time")
    time = energy_df["Time"]
    for column_name in args.obs:
        if column_name in energy_df.columns():
            observable = energy_df[column_name]
            if args.time_series:
                time = energy_df["Time"]*0.002
                figure_name = column_name.replace(" ", "_").lower()
                figure_name = args.output_base + figure_name + ".png"
                plt.figre(dpi = 300)
                plt.plot(time, observable)
                plt.xlabel("Time (ns)")
                plt.ylabel(column_name)
                plt.savefig(figure_name)
            if args.hist:
                figure_name = column_name.replace(" ", "_").lower()
                figure_name = figure_name.replace("(", "")
                figure_name = figure_name.replace(")", "")
                figure_name = args.output_base + figure_name + ".png"
                plt.figre(dpi = 300)
                plt.hist(observabl, bins = 30, density = True)
                plt.xlabel(column_name)
                plt.ylabel("Density")
                plt.savefig(figure_name)



    






