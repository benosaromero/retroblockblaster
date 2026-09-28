# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

from __future__ import print_function

import hashlib
import os

from backend.made_pipeline import (
    enumerate_filtered_blocks,
    filter_made_for_synthon_strict,
    run_validated_retro_decomposition,
    write_enumerated_products,
    write_filtered_made_rows,
)


# Supported MADE filtering modes.
#
# one_step:
#   Keep MADE candidates that directly satisfy the target synthon.
#
# two_step:
#   Keep a broader set that can satisfy the target synthon through the
#   two-step filtering logic in backend.made_pipeline.
SUPPORTED_FILTERS = ("one_step", "two_step")


def safe_key(synthon_id, smiles):
    """Build a stable folder key for a synthon ID and reference SMILES.

    The hash prevents long or chemically awkward SMILES strings from becoming
    folder names, while still keeping each synthon/reference pair distinct.
    """
    digest = hashlib.sha1(
        smiles.encode("utf-8")
    ).hexdigest()[:10]

    return "{}-{}".format(
        str(synthon_id).replace(".", "-"),
        digest,
    )


def write_raw_manifest(path, rows):
    """Write a manifest describing the raw MADE candidate files used."""
    columns = [
        "synthon_id",
        "reference_smiles",
        "filename",
        "path",
        "size_bytes",
        "warning",
    ]

    with open(path, "w") as output:
        output.write("\t".join(columns) + "\n")

        for row in rows:
            output.write(
                "\t".join(
                    str(row.get(column, ""))
                    for column in columns
                )
                + "\n"
            )


def run_made_workflow(
    target_smiles,
    output_directory,
    filter_mode,
):
    """Run the MADE analog workflow for one target molecule.

    High-level steps:
      1. Run Retro-Blockblaster retrosynthesis.
      2. Keep only routes whose synthon mapping validates.
      3. For each synthon, pull/filter MADE candidates.
      4. Enumerate products for routes where every synthon has candidates.
      5. Write summary TSV files and return a JSON-serializable result.
    """
    if filter_mode not in SUPPORTED_FILTERS:
        raise ValueError(
            "filter_mode must be one_step or two_step; "
            "anon_dist is unavailable because SmallWorld "
            "does not expose incrMADE"
        )

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    retro_directory = os.path.join(
        output_directory,
        "retro",
    )
    filter_root = os.path.join(
        output_directory,
        "filters",
    )
    enumeration_root = os.path.join(
        output_directory,
        "enumeration",
    )

    retro = run_validated_retro_decomposition(
        target_smiles,
        retro_directory,
    )

    validated_routes = [
        route
        for route in retro["routes"]
        if route.get("synthon_mapping_status")
        == "validated"
    ]

    # Cache repeated synthon lookups so the same synthon/reference pair is not
    # filtered more than once across routes.
    filter_cache = {}

    raw_manifest = []
    raw_seen = set()
    filtered_rows = []
    all_products = []
    all_matches = []
    route_results = []

    for route in validated_routes:
        block_lists = {}
        route_ready = True

        for synthon in route["synthons"]:
            synthon_id = synthon["synthon_id"]
            reference_smiles = synthon["smiles"]
            cache_key = (
                str(synthon_id),
                reference_smiles,
            )

            if cache_key not in filter_cache:
                filter_directory = os.path.join(
                    filter_root,
                    safe_key(
                        synthon_id,
                        reference_smiles,
                    ),
                )

                filter_cache[cache_key] = (
                    filter_made_for_synthon_strict(
                        synthon_smiles=reference_smiles,
                        synthon_id=synthon_id,
                        output_directory=filter_directory,
                    )
                )

            filtered = filter_cache[cache_key]
            raw = filtered["raw_made_file"]
            raw_key = (
                str(synthon_id),
                raw["path"],
            )

            if raw_key not in raw_seen:
                raw_seen.add(raw_key)
                raw_manifest.append({
                    "synthon_id": synthon_id,
                    "reference_smiles": reference_smiles,
                    "filename": raw["filename"],
                    "path": raw["path"],
                    "size_bytes": raw["size_bytes"],
                    "warning": raw["warning"],
                })

            selected_key = (
                "one_step_matches"
                if filter_mode == "one_step"
                else "two_step_matches"
            )
            selected = filtered[selected_key]

            for row in selected:
                combined = dict(row)
                combined["reaction_id"] = route[
                    "reaction_id"
                ]
                combined["pair_id"] = route["pair_id"]
                filtered_rows.append(combined)

            block_lists[str(synthon_id)] = selected

            if not selected:
                route_ready = False

        if not route_ready:
            route_results.append({
                "reaction_id": route["reaction_id"],
                "status": "skipped",
                "reason": (
                    "At least one synthon had no filtered "
                    "MADE candidates"
                ),
            })
            continue

        enumeration_directory = os.path.join(
            enumeration_root,
            "{}-{}".format(
                route["reaction_id"],
                route["pair_id"],
            ),
        )

        enumeration = enumerate_filtered_blocks(
            reaction_id=route["reaction_id"],
            block_lists=block_lists,
            target_smiles=target_smiles,
            output_directory=enumeration_directory,
        )

        all_products.extend(
            enumeration.get("products", [])
        )
        all_matches.extend(
            enumeration.get("matches", [])
        )
        route_results.append({
            "reaction_id": route["reaction_id"],
            "status": enumeration["status"],
            "combination_count": enumeration.get(
                "combination_count",
                0,
            ),
            "product_count": enumeration.get(
                "product_count",
                0,
            ),
            "match_count": enumeration.get(
                "match_count",
                0,
            ),
        })

    raw_path = os.path.join(
        output_directory,
        "05_made_bb_candidates_all.txt",
    )
    filtered_path = os.path.join(
        output_directory,
        "06_made_bb_candidates_filtered.txt",
    )
    products_path = os.path.join(
        output_directory,
        "07_made_enumerated_products.txt",
    )
    matches_path = os.path.join(
        output_directory,
        "08_made_inchikey_matches.txt",
    )

    write_raw_manifest(raw_path, raw_manifest)
    write_filtered_made_rows(
        filtered_path,
        filtered_rows,
    )
    write_enumerated_products(
        products_path,
        all_products,
    )
    write_enumerated_products(
        matches_path,
        all_matches,
    )

    return {
        "filter_mode": filter_mode,
        "retro_route_count": len(retro["routes"]),
        "validated_route_count": len(validated_routes),
        "filtered_candidate_count": len(filtered_rows),
        "product_count": len(all_products),
        "match_count": len(all_matches),
        "routes": route_results,
        "files": {
            "made_all": raw_path,
            "made_filtered": filtered_path,
            "made_products": products_path,
            "made_matches": matches_path,
        },
    }
