import argparse
import os
import subprocess
import sys

def run_name(value):
    return value[:-4] if value.endswith(".tsv") else value

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", default="dball")
    parser.add_argument("run")
    args = parser.parse_args()

    run = run_name(args.run)
    syn_tsv = run + "_syn.tsv"
    bb_dir = run + "_bb"

    if not os.path.isfile(syn_tsv):
        raise SystemExit("Missing synthons file: " + syn_tsv)

    subprocess.check_call([
        sys.executable,
        "prepare_dball_reagent_matches.py",
        syn_tsv,
        bb_dir,
        args.db,
    ])

    print("Wrote", bb_dir)

if __name__ == "__main__":
    main()
