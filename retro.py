import os
import subprocess
import sys

def run_name(value):
    return value[:-4] if value.endswith(".tsv") else value

def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 retro.py RUN_NAME")

    run = run_name(sys.argv[1])
    input_tsv = run + ".tsv"
    syn_tsv = run + "_syn.tsv"
    ret_dir = run + "_ret"

    if not os.path.isfile(input_tsv):
        raise SystemExit("Missing input file: " + input_tsv)

    subprocess.check_call([sys.executable, "run_retro_batch.py", input_tsv, syn_tsv, ret_dir])

    if os.path.exists("job_ids.tsv"):
        os.replace("job_ids.tsv", run + "_jobs.tsv")
    if os.path.exists("job_summary.txt"):
        os.replace("job_summary.txt", run + "_summary.tsv")

    print("Wrote", syn_tsv)
    print("Wrote", ret_dir)

if __name__ == "__main__":
    main()
