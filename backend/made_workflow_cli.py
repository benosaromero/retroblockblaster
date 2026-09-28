# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

from __future__ import print_function

import json
import os
import sys


# This CLI wrapper lets the FastAPI app or terminal launch the MADE workflow
# in a separate Python process. Keeping it small makes failures easier to
# isolate from the web server.

PROJECT_ROOT = "/nfs/home/bromero/retro-blockblaster_1.0"

if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from backend.made_workflow import run_made_workflow


def main():
    if len(sys.argv) != 5:
        raise SystemExit(
            "Usage: made_workflow_cli.py "
            "SMILES FILTER_MODE OUTPUT_DIRECTORY RESULT_JSON"
        )

    smiles = sys.argv[1]
    filter_mode = sys.argv[2]
    output_directory = sys.argv[3]
    result_path = sys.argv[4]

    result = run_made_workflow(
        target_smiles=smiles,
        output_directory=output_directory,
        filter_mode=filter_mode,
    )

    # Write through a temporary file so callers never read a half-written JSON.
    temporary_path = result_path + ".tmp"

    with open(temporary_path, "w") as output:
        json.dump(result, output, indent=2)

    os.replace(temporary_path, result_path)


if __name__ == "__main__":
    main()
