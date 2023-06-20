import signac
import argparse
import glob
import shutil
from signac_init import replace_all_pattern


def parse_args():
    parser = argparse.ArgumentParser(
        description = "A script to add files to signac projects",
    )

    parser.add_argument(
        "-f", "--file_path",
        type = str,
        nargs = "+",
        help = "files to add to signac projects"
    )

    parser.add_argument(
        "-r", "--replace",
        type = str,
        nargs = "+",
        help = "List of string needed to repalce in added files. \
            Currently this only works for WALKER_DIRS and strings \
            matching the statepoint variables."
    )

    return parser.parse_args()

def main():

    args = parse_args()
    project = signac.get_project()

    for job in project.find_jobs():
        for file_path, replace in zip(args.file_path, args.replace):
            filename = file_path.split("/")[-1]
            n_walkers = len(glob.glob(job.fn("WALKER*")))
            if replace == "WALKER_DIRS":
                walker_dirs = " ".join(["WALKER" + str(walker_id) for walker_id in range(n_walkers)])
                shutil.copy(file_path, job.fn(filename))
                replace_all_pattern("WALKER_DIRS", walker_dirs, job.fn(filename))
            if replace in job.sp.keys():
                walker_dirs = ["WALKER" + str(walker_id) for walker_id in range(n_walkers)]
                for walker_dir in walker_dirs:
                    shutil.copy(file_path, job.fn(os.path.join(walker_dir, filename)))
                    replace_all_pattern(replace, str(job.sp[replace]), job.fn(os.path.join(walker_dir, filename)))

            
if __name__ == "__main__":
    main()