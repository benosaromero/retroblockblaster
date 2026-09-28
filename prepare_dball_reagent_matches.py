# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

import csv
import json
import os
import sys
import urllib.parse
import urllib.request

# This script reads the combined retrosynthesis/synthon table and creates
# editable reagent-list TSV files for each input molecule and reaction route.
#
# Input:
#   all_synthons.tsv-style file, usually created from Retro-Blockblaster output.
#
# Output:
#   cli_dball_intermediates/
#     input_id/
#       reaction_id_pair_id/
#         route_info.json
#         bb1_synthon_<id>_matches.tsv
#         bb2_synthon_<id>_matches.tsv
#         ...
#
# Each *_matches.tsv file is meant to be manually reviewable/editable before
# enumeration. Deleting rows from those files restricts which building blocks
# are used in the later enumeration step.

BASE_URL = "http://127.0.0.1:8002"

# Command line usage:
#   python3 prepare_dball_reagent_matches.py input_synthons.tsv output_folder
#
# Defaults are kept for the common workflow.
INPUT_TSV = sys.argv[1] if len(sys.argv) > 1 else "all_synthons.tsv"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "cli_dball_intermediates"
DB_SOURCE = sys.argv[3] if len(sys.argv) > 3 else "dball"

# John requested a permissive cutoff so analogs are still found when exact
# purchasable building blocks are unavailable.
TC_CUTOFF = 0.4

# Maximum number of purchasable dball matches requested per synthon.
SYNTHON_LIMIT = 10000

os.makedirs(OUT_DIR, exist_ok=True)


def safe_name(text):
    """Convert arbitrary IDs into filenames that are safe to use in shell paths."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(text))


def get_json(url):
    """Fetch JSON from the local FastAPI server."""
    with urllib.request.urlopen(url) as response:
        return json.loads(response.read().decode("utf-8"))


def write_matches_tsv(path, results, disconnected=None):
    """Write one editable reagent TSV for one synthon slot.

    The rank 0 row is the disconnected retrosynthesis synthon. It may not be
    purchasable, but it is included so we can compare:
      1. fully purchasable enumerations
      2. enumerations where one or more retrosynthetic synthons are treated as
         perfect symbolic building blocks
    """
    with open(path, "w") as out:
        out.write("rank\ttc_score\tbb_id\tsmiles\tsynthon_id\tsynthon_name\tsource\n")

        if disconnected:
            out.write("0\t1.0\t{}\t{}\t{}\t{}\tretrosynthesis\n".format(
                disconnected.get("bb_id", ""),
                disconnected.get("smiles", ""),
                disconnected.get("synthon_id", ""),
                disconnected.get("synthon_name", ""),
            ))

        for i, item in enumerate(results, start=1):
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                i,
                item.get("tc_score", ""),
                item.get("bb_id", ""),
                item.get("smiles", ""),
                item.get("synthon_id", ""),
                item.get("synthon_name", ""),
                item.get("db", ""),
            ))


with open(INPUT_TSV) as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row_index, row in enumerate(reader, start=1):
        input_id = row.get("input_id", "")
        job_id = row.get("job_id", "")
        reaction_id = row.get("reaction_id", "")
        pair_id = row.get("pair_id", "")

        # Rows without a reaction and pair ID cannot be mapped to a route folder.
        if not reaction_id or not pair_id:
            continue

        route_dir = os.path.join(
            OUT_DIR,
            safe_name(input_id),
            "{}_{}".format(safe_name(reaction_id), safe_name(pair_id)),
        )
        os.makedirs(route_dir, exist_ok=True)

        route_info = {
            "input_id": input_id,
            "job_id": job_id,
            "reaction_id": reaction_id,
            "reaction_name": row.get("reaction_name", ""),
            "pair_id": pair_id,
            "reagent_count": row.get("reagent_count", ""),
            "original_smiles": row.get("original_smiles", ""),
            "source": DB_SOURCE,
            "tc_cutoff": TC_CUTOFF,
            "synthons": [],
        }

        print("Saving matches", input_id, reaction_id, pair_id)

        for slot in ["bb1", "bb2", "bb3", "bb4"]:
            synthon_id = row.get(slot + "_synthon_id", "").strip()
            reference_smiles = row.get(slot + "_smiles", "").strip()

            if not synthon_id:
                continue

            # Ask the local API for purchasable dball matches similar to the
            # disconnected synthon from this retrosynthetic route.
            query = urllib.parse.urlencode({
                "source": DB_SOURCE,
                "limit": str(SYNTHON_LIMIT),
                "reference_smiles": reference_smiles,
                "tc_cutoff": str(TC_CUTOFF),
            })

            url = "{}/synthon/{}?{}".format(
                BASE_URL,
                urllib.parse.quote(synthon_id),
                query,
            )

            data = get_json(url)
            results = data.get("results", [])

            os.makedirs(route_dir, exist_ok=True)

            matches_path = os.path.join(
                route_dir,
                "{}_synthon_{}_matches.tsv".format(slot, safe_name(synthon_id)),
            )

            disconnected = {
                "bb_id": "DISCONNECTED_{}".format(slot.upper()),
                "smiles": reference_smiles,
                "synthon_id": synthon_id,
                "synthon_name": row.get(slot + "_synthon_name", ""),
            }

            write_matches_tsv(matches_path, results, disconnected=disconnected)

            route_info["synthons"].append({
                "slot": slot,
                "synthon_id": synthon_id,
                "reference_smiles": reference_smiles,
                "match_count": len(results),
                "matches_file": matches_path,
                "message": data.get("message", ""),
            })

        with open(os.path.join(route_dir, "route_info.json"), "w") as out:
            json.dump(route_info, out, indent=2)

print("Done. Intermediate files saved in", OUT_DIR)
