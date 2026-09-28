# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

import csv
import os
import sys
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# This script ranks enumerated analog products by similarity to the original
# input/query molecules.
#
# Input:
#   molecules.tsv
#     The original query molecules. Expected columns:
#       input_id or id
#       smiles
#
#   all_selected_dball_products.tsv
#     Combined enumerated products from the selected reagent-list workflow.
#
# Output:
#   all_selected_dball_products_sorted.tsv
#     Combined ranked products for all queries.
#
#   cli_sorted_products_by_query/
#     One ranked TSV per input molecule.
#
# Similarity:
#   ECFP4 / Morgan radius 2 fingerprints with Tanimoto coefficient.

# Command line usage:
#   python3 rank_products_by_similarity.py molecules.tsv products.tsv output.tsv output_dir
#
# Defaults are kept for the common workflow.
MOLECULES_TSV = sys.argv[1] if len(sys.argv) > 1 else "molecules.tsv"
PRODUCTS_TSV = sys.argv[2] if len(sys.argv) > 2 else "all_selected_dball_products.tsv"
OUT_TSV = sys.argv[3] if len(sys.argv) > 3 else "all_selected_dball_products_sorted.tsv"
OUT_DIR = sys.argv[4] if len(sys.argv) > 4 else "cli_sorted_products_by_query"

os.makedirs(OUT_DIR, exist_ok=True)


def fp_from_smiles(smiles):
    """Return an ECFP4 fingerprint for a SMILES string, or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def tc(query_fp, smiles):
    """Compute Tanimoto similarity between a query fingerprint and SMILES."""
    fp = fp_from_smiles(smiles)
    if query_fp is None or fp is None:
        return ""
    return round(DataStructs.TanimotoSimilarity(query_fp, fp), 4)


def add_building_block_columns(row):
    blocks = str(row.get("building_blocks", "")).split("|")
    for index in range(4):
        column = "bb{}_smiles".format(index + 1)
        row[column] = blocks[index] if index < len(blocks) else ""
    return row


queries = {}

with open(MOLECULES_TSV) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        input_id = (row.get("input_id", "") or row.get("id", "")).strip()
        smiles = row.get("smiles", "").strip()

        if input_id and smiles:
            queries[input_id] = {
                "smiles": smiles,
                "fp": fp_from_smiles(smiles),
            }

products_by_input = {}

with open(PRODUCTS_TSV) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        input_id = row.get("input_id", "").strip()
        if not input_id:
            continue
        products_by_input.setdefault(input_id, []).append(row)

all_rows = []

for input_id, query in sorted(queries.items()):
    query_smiles = query["smiles"]
    query_fp = query["fp"]

    rows = []

    # Keep the original query molecule at rank 0 as a reference row.
    rows.append(add_building_block_columns({
        "input_id": input_id,
        "rank": 0,
        "similarity_to_query": 1.0,
        "is_query": "yes",
        "route_key": "",
        "product_code": "QUERY",
        "product_smiles": query_smiles,
        "tranche": "",
        "building_blocks": "",
    }))

    seen = set([query_smiles])
    scored_products = []

    for product in products_by_input.get(input_id, []):
        smiles = product.get("product_smiles", "").strip()
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)

        row = {
            "input_id": input_id,
            "similarity_to_query": tc(query_fp, smiles),
            "is_query": "no",
            "route_key": product.get("route_key", ""),
            "product_code": product.get("product_code", ""),
            "product_smiles": smiles,
            "tranche": product.get("tranche", ""),
            "building_blocks": product.get("building_blocks", ""),
        }
        scored_products.append(add_building_block_columns(row))

    scored_products.sort(
        key=lambda row: row["similarity_to_query"] if row["similarity_to_query"] != "" else -1,
        reverse=True,
    )

    for rank, row in enumerate(scored_products, start=1):
        row["rank"] = rank
        rows.append(row)

    per_query_path = os.path.join(OUT_DIR, "{}_sorted_products.tsv".format(input_id))

    with open(per_query_path, "w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=[
            "input_id",
            "rank",
            "similarity_to_query",
            "is_query",
            "route_key",
            "product_code",
            "product_smiles",
            "tranche",
            "building_blocks",
            "bb1_smiles",
            "bb2_smiles",
            "bb3_smiles",
            "bb4_smiles",
        ])
        writer.writeheader()
        writer.writerows(rows)

    all_rows.extend(rows)

RANK_FIELDNAMES = [
    "input_id",
    "rank",
    "similarity_to_query",
    "is_query",
    "route_key",
    "product_code",
    "product_smiles",
    "tranche",
    "building_blocks",
    "bb1_smiles",
    "bb2_smiles",
    "bb3_smiles",
    "bb4_smiles",
]


def numeric_similarity(row):
    try:
        return float(row.get("similarity_to_query", ""))
    except Exception:
        return -1.0


def is_disconnected_row(row):
    text = " ".join(
        str(row.get(column, ""))
        for column in [
            "product_smiles",
            "product_code",
            "tranche",
            "building_blocks",
        ]
    ).lower()

    return "disconnected" in text


def write_rank_subset(path, rows, fieldnames, tc_cutoff, keep_disconnected):
    kept = []

    for row in rows:
        if row.get("is_query") == "yes":
            kept.append(row)
            continue

        if numeric_similarity(row) < tc_cutoff:
            continue

        if not keep_disconnected and is_disconnected_row(row):
            continue

        kept.append(row)

    with open(path, "w") as out:
        writer = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(kept)

    return len(kept)


with open(OUT_TSV, "w") as out:
    writer = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=RANK_FIELDNAMES,
    )
    writer.writeheader()
    writer.writerows(all_rows)

out_base, out_ext = os.path.splitext(OUT_TSV)

if not out_ext:
    out_ext = ".tsv"

tc06_all_path = out_base + "_tc06_all" + out_ext
tc04_connected_path = out_base + "_tc04_connected" + out_ext

tc06_count = write_rank_subset(
    tc06_all_path,
    all_rows,
    RANK_FIELDNAMES,
    tc_cutoff=0.6,
    keep_disconnected=True,
)

tc04_count = write_rank_subset(
    tc04_connected_path,
    all_rows,
    RANK_FIELDNAMES,
    tc_cutoff=0.4,
    keep_disconnected=False,
)

print("Wrote", OUT_TSV)
print("Wrote", tc06_all_path, "rows:", tc06_count)
print("Wrote", tc04_connected_path, "rows:", tc04_count)
print("Wrote per-query files in", OUT_DIR)
