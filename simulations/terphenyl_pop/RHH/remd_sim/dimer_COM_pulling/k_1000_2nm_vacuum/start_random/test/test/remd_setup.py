# Setup Script for REMD replicas at varying temperatures

import os
import shutil
import argparse
import sys
import numpy as np
def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("-N","--n_replicas",
                             type = int,
                             help = "Number of replicas",
                            )
    args_parser.add_argument("--t_range",
                             type = float,
                             nargs = 2,
                             help = "Temperature ranges over which to create for replicas",
                             required = True
                            )
    args_parser.add_argument("--mdps",
                             type = str,
                             nargs = "+",
                             help = "List of mdps to move and edit in each replica",
                             required = True
                            )
    args_parser.add_argument("-R","--common_ratio",
                             type = float,
                             help = "Common_ratio",
                            )
    args_parser.add_argument("--sim_id",
                             type = str,
                             help="base name of dirs containing simulation",
                            )
    args_parser.add_argument("--extra_files",
                             type = str,
                             nargs = "+",
                             help="Additional files requied for REMD simulations"
                            )
    args = args_parser.parse_args()

    if sum([inp is not None for inp in [args.n_replicas, args.t_range, args.common_ratio]]) < 2:
        print(sum([inp is None for inp in [args.n_replicas, args.t_range, args.common_ratio]])) 
        print("Too few parameters specified! Please specify ONLY 2 of the following inputs:", file=sys.stderr)
        print("--n_replicas --t_range --common_ratio", file = sys.stderr)
        sys.exit(1)
    
    if sum([inp is not None for inp in [args.n_replicas, args.t_range, args.common_ratio]]) > 2:
        print("Too many parameters specified! Please specify ONLY 2 of the following inputs:", file=sys.stderr)
        print("--n_replicas --t_range --common_ratio", file = sys.stderr)
        sys.exit(1)
    
    if args.n_replicas == None:
        args.n_replicas = int(np.log(args.t_range[1]/args.t_range[0])/np.log(args.common_ratio)-1)

    if args.common_ratio == None:
        args.common_ratio = np.power(args.t_range[1]/args.t_range[0], 1/(args.n_replicas-1))

    for i in range(args.n_replicas):
        t_i = args.t_range[0] * args.common_ratio ** i
        sim_path = args.sim_id + str(i)
        if os.path.isdir(sim_path):
            n = 1
            print("Backoff!")
            while os.path.isdir("#" + sim_path + "." + str(n) + "#"):
                n += 1
            os.rename(sim_path, "#" + sim_path + "." + str(n) + "#")
        os.makedirs(sim_path)
        for mdp in args.mdps:
            with open(mdp, "r") as r:
                with open(os.path.join(sim_path, mdp), "w") as w:
                    for line in r:
                        if "TEMP" in line:
                            w.write(line.replace("TEMP", str(t_i)))
                        else:
                            w.write(line)
        for extra in args.extra_files:
            shutil.copy(extra, os.path.join(sim_path, extra))
        




if __name__ == "__main__":
    main()
