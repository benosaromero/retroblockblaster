# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

import csv
import glob
import json
import os
import sys
import urllib.request

# This script enumerates products from editable reagent-list TSV files.
#
# Typical workflow:
#   1. Run prepare_dball_reagent_matches.py to create cli_dball_intermediates/.
#   2. Manually delete unwanted rows from any *_matches.tsv files.
#   3. Run this script to enumerate only the remaining selected rows.
#
# Input:
#   cli_dball_intermediates/
#
# Output:
#   cli_dball_selected_enumerations/
#     input_id/
#       route_key/
#         enumerated_products.json
#         enumerated_products.tsv

BASE_URL = "http://127.0.0.1:8002"

# Command line usage:
#   python3 enumerate_selected_reagent_lists.py input_intermediates_dir output_dir
#
# Defaults are kept for the common workflow.
INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "cli_dball_intermediates"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "cli_dball_selected_enumerations"

# Safety limit passed to the API so accidental huge reagent lists do not explode.
MAX_COMBOS = 1000000

os.makedirs(OUT_DIR, exist_ok=True)


def safe_name(text):
    """Convert arbitrary IDs into filenames that are safe to use in shell paths."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(text))


def post_json(url, payload):
    """POST JSON to the local FastAPI server and return parsed JSON."""
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    try:
        with urllib.request.urlopen(request) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception as error:
        return {
            "error": str(error),
            "url": url,
            "payload": payload,
        }


def read_matches(path):
    """Read one edited *_matches.tsv file into API reagent block dictionaries.

    Rows removed from the TSV are not included here, which is how manual curation
    controls the enumeration.
    """
    blocks = []

    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            smiles = row.get("smiles", "").strip()
            if not smiles:
                continue

            blocks.append({
                "smiles": smiles,
                "bb_id": row.get("bb_id", "").strip(),
                "synthon_id": row.get("synthon_id", "").strip(),
                "synthon_name": row.get("synthon_name", "").strip(),
                "db": row.get("source", "").strip(),
                "tc_score": row.get("tc_score", "").strip(),
            })

    return blocks


def write_products_tsv(path, products):
    """Write a compact TSV copy of enumerated products for easy inspection."""
    with open(path, "w") as out:
        out.write("product_code\tproduct_smiles\ttranche\tbuilding_blocks\n")

        for product in products:
            out.write("{}\t{}\t{}\t{}\n".format(
                product.get("product_code", ""),
                product.get("product_smiles", ""),
                product.get("tranche", ""),
                "|".join(product.get("building_blocks", [])),
            ))


route_info_paths = sorted(glob.glob(os.path.join(INPUT_DIR, "*", "*", "route_info.json")))

print("Found route folders:", len(route_info_paths))

for index, route_info_path in enumerate(route_info_paths, start=1):
    route_dir = os.path.dirname(route_info_path)

    with open(route_info_path) as handle:
        route_info = json.load(handle)

    input_id = route_info.get("input_id", "")
    job_id = route_info.get("job_id", "")
    pair_id = route_info.get("pair_id", "")
    reaction_id = route_info.get("reaction_id", "")

    route_key = "{}_{}_{}".format(
        safe_name(input_id),
        safe_name(reaction_id),
        safe_name(pair_id),
    )

    out_route_dir = os.path.join(OUT_DIR, safe_name(input_id), route_key)
    os.makedirs(out_route_dir, exist_ok=True)

    output_json = os.path.join(out_route_dir, "enumerated_products.json")
    output_tsv = os.path.join(out_route_dir, "enumerated_products.tsv")

    # This makes the script resumable after disconnects or long cluster runs.
    if os.path.exists(output_json):
        print("Skipping existing", index, route_key)
        continue

    block_lists = {}

    for match_path in sorted(glob.glob(os.path.join(route_dir, "bb*_synthon_*_matches.tsv"))):
        filename = os.path.basename(match_path)
        slot = filename.split("_", 1)[0]

        blocks = read_matches(match_path)

        if blocks:
            block_lists[slot] = blocks

    print("Enumerating", index, route_key, reaction_id, "synthon lists:", len(block_lists))

    if len(block_lists) < 2:
        result = {
            "error": "fewer than two editable synthon match files had usable rows",
            "route_info": route_info,
            "block_lists": block_lists,
        }
    else:
        result = post_json(
            BASE_URL + "/enumerate-all",
            {
                "reaction_id": reaction_id,
                "pair_id": pair_id + "_selected",
                "block_lists": block_lists,
                "max_combos": MAX_COMBOS,
            },
        )

        result["route_info"] = route_info
        result["source_intermediate_dir"] = route_dir

    with open(output_json, "w") as out:
        json.dump(result, out, indent=2)

    write_products_tsv(output_tsv, result.get("products", []))

print("Done. Selected enumerations saved in", OUT_DIR)
