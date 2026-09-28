import os
import subprocess
import sys

def run_name(value):
    return value[:-4] if value.endswith(".tsv") else value

def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 enum.py RUN_NAME")

    run = run_name(sys.argv[1])
    bb_dir = run + "_bb"
    enum_dir = run + "_enum"

    if not os.path.isdir(bb_dir):
        raise SystemExit("Missing building-block folder: " + bb_dir)

    subprocess.check_call([
        sys.executable,
        "enumerate_selected_reagent_lists.py",
        bb_dir,
        enum_dir,
    ])

    print("Wrote", enum_dir)

if __name__ == "__main__":
    main()
