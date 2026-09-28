# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

# MADE pipeline helpers.
#
# This module contains the lower-level functions used by backend/made_workflow.py:
#   - submit and poll ZINC26/Arthor MADE searches
#   - read precomputed incrMADE synthon files
#   - validate retrosynthetic synthons
#   - filter MADE candidates for one-step or two-step workflows
#   - enumerate filtered building blocks
#   - write workflow TSV outputs
#
# The higher-level orchestration lives in backend/made_workflow.py.

from __future__ import print_function

import json
import time
import urllib.parse
import urllib.request


ZINC26_URL = "https://zinc26.docking.org"
MADE_DATABASE = "Purchasable"
MADE_CATEGORY = "incoming"
MADE_CATALOG = "incrMADE"

POLL_INTERVAL_SECONDS = 10
SEARCH_TIMEOUT_SECONDS = 300


def get_json(url, timeout=30):
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def submit_made_search(smiles, search_type="sim"):
    parameters = urllib.parse.urlencode({
        "type": search_type,
        "database": MADE_DATABASE,
        "category": MADE_CATEGORY,
        "catalog": MADE_CATALOG,
        "query": smiles,
    })

    submission = get_json(
        ZINC26_URL + "/arthor/search?" + parameters
    )

    task_id = submission.get("task_id")

    if not task_id:
        raise RuntimeError(
            "ZINC26 returned no task_id: {}".format(submission)
        )

    return task_id


def poll_made_search(
    task_id,
    timeout_seconds=SEARCH_TIMEOUT_SECONDS,
    poll_interval=POLL_INTERVAL_SECONDS,
):
    result_url = (
        ZINC26_URL
        + "/arthor/result?"
        + urllib.parse.urlencode({"task_id": task_id})
    )
    deadline = time.time() + timeout_seconds

    while time.time() < deadline:
        result = get_json(result_url)
        status = result.get("status")

        if status == "SUCCESS":
            return result.get("results", [])

        if status == "FAILURE":
            raise RuntimeError(
                result.get("stderr")
                or "ZINC26 MADE search failed"
            )

        time.sleep(poll_interval)

    raise RuntimeError(
        "ZINC26 MADE search exceeded {} seconds".format(
            timeout_seconds
        )
    )


def search_made(smiles, search_type="sim"):
    task_id = submit_made_search(
        smiles=smiles,
        search_type=search_type,
    )
    rows = poll_made_search(task_id)

    return {
        "task_id": task_id,
        "results": rows,
    }


def exact_rows(rows):
    matches = []

    for row in rows:
        try:
            similarity = float(row.get("tc", 0))
        except (TypeError, ValueError):
            similarity = 0.0

        if similarity >= 0.999999:
            matches.append(row)

    return matches


def search_made_exact(smiles):
    search_result = search_made(
        smiles=smiles,
        search_type="sim",
    )

    return {
        "task_id": search_result["task_id"],
        "results": exact_rows(search_result["results"]),
    }


def write_rows(path, rows):
    columns = [
        "compound_id",
        "smiles",
        "tc",
        "ecfp4",
        "rdkitfp",
        "logp",
    ]

    with open(path, "w") as output:
        output.write("\t".join(columns) + "\n")

        for row in rows:
            values = [
                str(row.get(column, ""))
                for column in columns
            ]
            output.write("\t".join(values) + "\n")


MADE_SYNTHON_ROOT = (
    "/mnt/nfs/exj/Fe/synthons/incrMADE"
)


def validate_synthon_id(synthon_id):
    import re

    value = str(synthon_id).strip()

    if not re.match(r"^[0-9]+(?:\.[0-9]+)?$", value):
        raise ValueError(
            "Invalid synthon ID: {}".format(value)
        )

    return value


def find_made_synthon_file(synthon_id):
    import glob
    import os

    value = validate_synthon_id(synthon_id)
    pattern = os.path.join(
        MADE_SYNTHON_ROOT,
        "incrMADE_{}_".format(value) + "*",
    )

    matches = [
        path
        for path in glob.glob(pattern)
        if os.path.isfile(path)
    ]

    if not matches:
        raise IOError(
            "No incrMADE file found for synthon {}".format(
                value
            )
        )

    if len(matches) > 1:
        raise RuntimeError(
            "Multiple incrMADE files found for synthon {}: {}".format(
                value,
                matches,
            )
        )

    return matches[0]


def made_synthon_file_info(synthon_id):
    import os

    path = find_made_synthon_file(synthon_id)

    return {
        "synthon_id": validate_synthon_id(synthon_id),
        "path": path,
        "filename": os.path.basename(path),
        "size_bytes": os.path.getsize(path),
        "warning": (
            "This raw MADE building-block file may contain "
            "millions of rows."
        ),
    }


def iter_made_synthon_rows(synthon_id):
    path = find_made_synthon_file(synthon_id)

    with open(path, "r", errors="replace") as source:
        for line_number, line in enumerate(source, 1):
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            yield {
                "smiles": parts[0],
                "compound_id": parts[1],
                "synthon_id": parts[2],
                "line_number": line_number,
            }


RETRO_PROJECT_ROOT = (
    "/nfs/home/bromero/retro-blockblaster_1.0"
)
RETRO_JOB_PYTHON = (
    "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python"
)
RETRO_RUN_JOB = (
    RETRO_PROJECT_ROOT + "/run_job.py"
)


def parse_retro_output(output_path):
    import csv

    routes = []

    with open(
        output_path,
        "r",
        errors="replace",
    ) as source:
        reader = csv.DictReader(source, delimiter="\t")

        for row in reader:
            if not row.get("reaction_id"):
                continue

            synthons = []

            for slot_number in range(1, 5):
                prefix = "bb{}".format(slot_number)
                smiles = str(
                    row.get(prefix + "_smiles", "")
                ).strip()
                synthon_id = str(
                    row.get(prefix + "_synthon_id", "")
                ).strip()
                synthon_name = str(
                    row.get(prefix + "_synthon_name", "")
                ).strip()

                if not smiles or not synthon_id:
                    continue

                synthons.append({
                    "slot": slot_number,
                    "smiles": smiles,
                    "synthon_id": synthon_id,
                    "synthon_name": synthon_name,
                })

            routes.append({
                "reaction_id": row.get("reaction_id", ""),
                "reaction_name": row.get("reaction_name", ""),
                "pair_id": row.get("pair_id", ""),
                "original_smiles": row.get(
                    "original_smiles",
                    "",
                ),
                "synthons": synthons,
            })

    return routes


def run_retro_decomposition(smiles, job_directory):
    import os
    import subprocess

    if not os.path.isdir(job_directory):
        os.makedirs(job_directory)

    input_path = os.path.join(
        job_directory,
        "retro_input.smi",
    )
    output_path = os.path.join(
        job_directory,
        "retro_output.txt",
    )

    with open(input_path, "w") as output:
        output.write(
            "{} input_molecule\n".format(smiles)
        )

    environment = os.environ.copy()
    environment["PYTHONPATH"] = "/nfs/exj/Fe/scripts"

    subprocess.check_call(
        [
            RETRO_JOB_PYTHON,
            RETRO_RUN_JOB,
            input_path,
            output_path,
        ],
        cwd=RETRO_PROJECT_ROOT,
        env=environment,
    )

    return {
        "input_path": input_path,
        "output_path": output_path,
        "routes": parse_retro_output(output_path),
    }


AUTOMATIC_REPLACEMENTS = {
    "C": ["N", "O", "S"],
    "N": ["C", "O", "S"],
    "O": ["N", "S"],
    "S": ["N", "O"],
    "F": ["Cl", "Br", "I"],
    "Cl": ["F", "Br", "I"],
    "Br": ["F", "Cl", "I"],
    "I": ["F", "Cl", "Br"],
    "B": ["C", "N"],
    "P": ["N", "S"],
}


def build_automatic_modifications(smiles):
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError(
            "Invalid synthon SMILES: {}".format(smiles)
        )

    modifications = {}

    for atom in molecule.GetAtoms():
        replacements = AUTOMATIC_REPLACEMENTS.get(
            atom.GetSymbol(),
            [],
        )

        if replacements:
            modifications[str(atom.GetIdx())] = replacements

    if not modifications:
        raise ValueError(
            "No eligible atoms found in synthon {}".format(
                smiles
            )
        )

    return modifications


def generate_synthon_close_analogs(
    synthon_smiles,
    output_directory,
    modifications=None,
):
    import csv
    import os
    from pathlib import Path

    from backend.analog_pipeline import (
        run_close_analog_generator,
    )

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    if modifications is None:
        modifications = build_automatic_modifications(
            synthon_smiles
        )

    result = run_close_analog_generator(
        parent_smiles=synthon_smiles,
        modifications_dict=modifications,
        job_dir=Path(output_directory),
    )

    if result.get("error"):
        raise RuntimeError(
            result.get("stderr")
            or result.get("error")
        )

    csv_path = result.get("csv")

    if not csv_path or not os.path.isfile(csv_path):
        raise IOError(
            "Close-analog CSV was not generated"
        )

    one_step = []
    two_step = []

    with open(csv_path, "r", errors="replace") as source:
        for row in csv.DictReader(source):
            mutation_name = row.get("name", "")
            smiles = row.get("smiles", "").strip()

            if not smiles:
                continue

            step_count = mutation_name.count("_") + 1
            record = {
                "name": mutation_name,
                "smiles": smiles,
                "step_count": step_count,
            }

            if step_count == 1:
                one_step.append(record)
            elif step_count == 2:
                two_step.append(record)

    return {
        "parent_smiles": synthon_smiles,
        "modifications": modifications,
        "one_step": one_step,
        "two_step": two_step,
        "csv_path": csv_path,
    }


def smiles_inchikey(smiles):
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(str(smiles).strip())

    if molecule is None:
        return ""

    return Chem.MolToInchiKey(molecule)


def index_generated_analogs(generated_rows):
    indexed = {}

    for row in generated_rows:
        inchikey = smiles_inchikey(row.get("smiles", ""))

        if not inchikey:
            continue

        indexed[inchikey] = {
            "mutation_name": row.get("name", ""),
            "generated_smiles": row.get("smiles", ""),
            "step_count": row.get("step_count", ""),
        }

    return indexed


def match_made_rows_to_generated(
    made_rows,
    generated_rows,
    synthon_id,
):
    generated_by_key = index_generated_analogs(
        generated_rows
    )
    matches = []
    seen = set()

    for made_row in made_rows:
        made_smiles = made_row.get("smiles", "")
        inchikey = smiles_inchikey(made_smiles)

        if (
            not inchikey
            or inchikey not in generated_by_key
            or (inchikey, made_row.get("compound_id")) in seen
        ):
            continue

        generated = generated_by_key[inchikey]
        seen.add((inchikey, made_row.get("compound_id")))

        matches.append({
            "inchikey": inchikey,
            "smiles": made_smiles,
            "compound_id": made_row.get("compound_id", ""),
            "synthon_id": str(synthon_id),
            "step_count": generated["step_count"],
            "mutation_name": generated["mutation_name"],
            "tc": made_row.get("tc", ""),
            "ecfp4": made_row.get("ecfp4", ""),
            "rdkitfp": made_row.get("rdkitfp", ""),
            "logp": made_row.get("logp", ""),
        })

    return matches


def write_filtered_made_rows(path, rows):
    columns = [
        "inchikey",
        "smiles",
        "compound_id",
        "synthon_id",
        "step_count",
        "mutation_name",
        "tc",
        "ecfp4",
        "rdkitfp",
        "logp",
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


def filter_made_for_synthon(
    synthon_smiles,
    synthon_id,
    output_directory,
    modifications=None,
):
    import os

    safe_synthon_id = validate_synthon_id(
        synthon_id
    ).replace(".", "-")

    analog_directory = os.path.join(
        output_directory,
        "synthon-{}-analogs".format(
            safe_synthon_id
        ),
    )

    generated = generate_synthon_close_analogs(
        synthon_smiles=synthon_smiles,
        output_directory=analog_directory,
        modifications=modifications,
    )

    search_result = search_made(
        smiles=synthon_smiles,
        search_type="sim",
    )
    made_rows = search_result["results"]

    one_step_matches = match_made_rows_to_generated(
        made_rows=made_rows,
        generated_rows=generated["one_step"],
        synthon_id=synthon_id,
    )
    two_step_matches = match_made_rows_to_generated(
        made_rows=made_rows,
        generated_rows=generated["two_step"],
        synthon_id=synthon_id,
    )

    one_step_path = os.path.join(
        output_directory,
        "made_synthon_{}_one_step.txt".format(
            safe_synthon_id
        ),
    )
    two_step_path = os.path.join(
        output_directory,
        "made_synthon_{}_two_step.txt".format(
            safe_synthon_id
        ),
    )

    write_filtered_made_rows(
        one_step_path,
        one_step_matches,
    )
    write_filtered_made_rows(
        two_step_path,
        two_step_matches,
    )

    return {
        "synthon_id": str(synthon_id),
        "synthon_smiles": synthon_smiles,
        "task_id": search_result["task_id"],
        "zinc26_candidate_count": len(made_rows),
        "generated_one_step_count": len(
            generated["one_step"]
        ),
        "generated_two_step_count": len(
            generated["two_step"]
        ),
        "one_step_match_count": len(
            one_step_matches
        ),
        "two_step_match_count": len(
            two_step_matches
        ),
        "one_step_path": one_step_path,
        "two_step_path": two_step_path,
        "raw_made_file": made_synthon_file_info(
            synthon_id
        ),
    }


SYNTHONS_JSON_PATH = (
    "/nfs/exj/Fe/scripts/common_db/export/"
    "Synthons.json"
)


def load_synthon_definitions():
    with open(SYNTHONS_JSON_PATH, "r") as source:
        return json.load(source)


def get_synthon_definition(definitions, synthon_id):
    value = str(synthon_id).strip()

    if isinstance(definitions, dict):
        definition = definitions.get(value)

        if definition is None and value.isdigit():
            definition = definitions.get(int(value))

        return definition

    for definition in definitions:
        candidate_id = definition.get(
            "ID",
            definition.get("id"),
        )

        if str(candidate_id) == value:
            return definition

    return None


def molecule_matches_synthon(smiles, definition):
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None or not definition:
        return False

    inclusion_patterns = definition.get(
        "inSMARTS",
        [],
    )
    exclusion_patterns = definition.get(
        "exSMARTS",
        [],
    )

    inclusion_match = False

    for pattern_text in inclusion_patterns:
        pattern = Chem.MolFromSmarts(pattern_text)

        if (
            pattern is not None
            and molecule.HasSubstructMatch(pattern)
        ):
            inclusion_match = True
            break

    if not inclusion_match:
        return False

    for pattern_text in exclusion_patterns:
        pattern = Chem.MolFromSmarts(pattern_text)

        if (
            pattern is not None
            and molecule.HasSubstructMatch(pattern)
        ):
            return False

    return True


def reassign_route_synthons(route, definitions):
    import itertools

    original_synthons = route.get("synthons", [])
    expected_ids = [
        synthon["synthon_id"]
        for synthon in original_synthons
    ]

    if len(expected_ids) < 2:
        route["synthon_mapping_status"] = (
            "not_needed"
        )
        return route

    best_assignment = None
    best_score = -1

    for assignment in itertools.permutations(
        expected_ids,
        len(original_synthons),
    ):
        score = 0

        for synthon, assigned_id in zip(
            original_synthons,
            assignment,
        ):
            definition = get_synthon_definition(
                definitions,
                assigned_id,
            )

            if molecule_matches_synthon(
                synthon["smiles"],
                definition,
            ):
                score += 1

        if score > best_score:
            best_score = score
            best_assignment = assignment

    reassigned = []

    for synthon, assigned_id in zip(
        original_synthons,
        best_assignment,
    ):
        definition = get_synthon_definition(
            definitions,
            assigned_id,
        )
        updated = dict(synthon)
        updated["original_synthon_id"] = synthon[
            "synthon_id"
        ]
        updated["synthon_id"] = str(assigned_id)
        updated["synthon_name"] = (
            definition.get("Name", "")
            if definition
            else ""
        )
        updated["smarts_match"] = (
            molecule_matches_synthon(
                updated["smiles"],
                definition,
            )
        )
        reassigned.append(updated)

    corrected = dict(route)
    corrected["synthons"] = reassigned
    corrected["synthon_mapping_score"] = best_score
    corrected["synthon_mapping_status"] = (
        "validated"
        if best_score == len(reassigned)
        else "unresolved"
    )

    return corrected


def run_validated_retro_decomposition(
    smiles,
    job_directory,
):
    result = run_retro_decomposition(
        smiles,
        job_directory,
    )
    definitions = load_synthon_definitions()
    result["routes"] = [
        reassign_route_synthons(
            route,
            definitions,
        )
        for route in result["routes"]
    ]

    return result


def restrict_made_rows_to_synthon(
    made_rows,
    synthon_id,
):
    candidates_by_id = {
        str(row.get("compound_id", "")).strip(): row
        for row in made_rows
        if row.get("compound_id")
    }

    if not candidates_by_id:
        return []

    restricted = []
    found_ids = set()
    path = find_made_synthon_file(synthon_id)

    with open(path, "r", errors="replace") as source:
        for line in source:
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            raw_smiles = parts[0]
            compound_id = parts[1]
            row_synthon_id = parts[2]

            if (
                compound_id not in candidates_by_id
                or row_synthon_id != str(synthon_id)
                or compound_id in found_ids
            ):
                continue

            row = dict(candidates_by_id[compound_id])
            row["smiles"] = raw_smiles
            row["synthon_id"] = row_synthon_id
            restricted.append(row)
            found_ids.add(compound_id)

            if len(found_ids) == len(candidates_by_id):
                break

    return restricted


def filter_made_for_synthon_strict(
    synthon_smiles,
    synthon_id,
    output_directory,
    modifications=None,
):
    import os

    safe_synthon_id = validate_synthon_id(
        synthon_id
    ).replace(".", "-")

    analog_directory = os.path.join(
        output_directory,
        "synthon-{}-analogs".format(
            safe_synthon_id
        ),
    )

    generated = generate_synthon_close_analogs(
        synthon_smiles=synthon_smiles,
        output_directory=analog_directory,
        modifications=modifications,
    )

    search_result = search_made(
        smiles=synthon_smiles,
        search_type="sim",
    )

    restricted_rows = restrict_made_rows_to_synthon(
        search_result["results"],
        synthon_id,
    )

    one_step_matches = match_made_rows_to_generated(
        made_rows=restricted_rows,
        generated_rows=generated["one_step"],
        synthon_id=synthon_id,
    )
    two_step_matches = match_made_rows_to_generated(
        made_rows=restricted_rows,
        generated_rows=generated["two_step"],
        synthon_id=synthon_id,
    )

    one_step_path = os.path.join(
        output_directory,
        "made_synthon_{}_one_step.txt".format(
            safe_synthon_id
        ),
    )
    two_step_path = os.path.join(
        output_directory,
        "made_synthon_{}_two_step.txt".format(
            safe_synthon_id
        ),
    )

    write_filtered_made_rows(
        one_step_path,
        one_step_matches,
    )
    write_filtered_made_rows(
        two_step_path,
        two_step_matches,
    )

    return {
        "synthon_id": str(synthon_id),
        "synthon_smiles": synthon_smiles,
        "task_id": search_result["task_id"],
        "zinc26_candidate_count": len(
            search_result["results"]
        ),
        "synthon_candidate_count": len(
            restricted_rows
        ),
        "generated_one_step_count": len(
            generated["one_step"]
        ),
        "generated_two_step_count": len(
            generated["two_step"]
        ),
        "one_step_matches": one_step_matches,
        "two_step_matches": two_step_matches,
        "one_step_path": one_step_path,
        "two_step_path": two_step_path,
        "raw_made_file": made_synthon_file_info(
            synthon_id
        ),
    }


REACTIONS_JSON_PATH = (
    "/nfs/exj/Fe/scripts/common_db/export/"
    "Reactions.json"
)
ENUM_SETTINGS_PATH = (
    "/nfs/exj/Fe/scripts/testreactions2/"
    "EXAMPLE_Library_Enumeration_File_"
    "hg2_noProtected_logP350_HAC35_RotB7_V3_A-AK1.json"
)


def load_enumeration_reaction(reaction_id):
    import sys

    scripts_path = "/nfs/exj/Fe/scripts"

    if scripts_path not in sys.path:
        sys.path.insert(0, scripts_path)

    from backend import enumerator_web as enum

    with open(REACTIONS_JSON_PATH, "r") as source:
        reactions = json.load(source)

    with open(ENUM_SETTINGS_PATH, "r") as source:
        settings = json.load(source)

    for name, value in settings[
        "Product_Restrictions_Map"
    ].items():
        setattr(
            enum.Restrictions_Identifiers,
            name,
            value,
        )

    reaction = enum.Extract_Rxns(
        reactions["Reaction_IDs"],
        reactions,
        [reaction_id],
    )[0]

    return enum, reaction, settings


def order_blocks_for_reaction(blocks, reaction):
    import itertools
    from rdkit import Chem

    templates = [
        reaction.buildRxn.GetReactantTemplate(index)
        for index in range(
            reaction.buildRxn.GetNumReactantTemplates()
        )
    ]

    for ordering in itertools.permutations(blocks):
        valid = True

        for block, template in zip(
            ordering,
            templates,
        ):
            molecule = Chem.MolFromSmiles(
                block.get("smiles", "")
            )

            if (
                molecule is None
                or not molecule.HasSubstructMatch(
                    template
                )
            ):
                valid = False
                break

        if valid:
            return ordering

    return None


def write_enumerated_products(path, rows):
    columns = [
        "product_inchikey",
        "product_smiles",
        "reaction_id",
        "tranche",
        "building_block_ids",
        "building_block_smiles",
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


def enumerate_filtered_blocks(
    reaction_id,
    block_lists,
    target_smiles,
    output_directory,
    maximum_combinations=100000,
):
    import contextlib
    import io
    import itertools
    import os

    enum, reaction, settings = (
        load_enumeration_reaction(reaction_id)
    )

    ordered_lists = []

    for synthon_id in reaction.Synthons:
        blocks = block_lists.get(
            str(synthon_id),
            [],
        )

        if not blocks:
            return {
                "status": "skipped",
                "reason": (
                    "No filtered MADE blocks for synthon {}"
                    .format(synthon_id)
                ),
                "combination_count": 0,
                "products": [],
                "matches": [],
            }

        ordered_lists.append(blocks)

    combination_count = 1

    for blocks in ordered_lists:
        combination_count *= len(blocks)

    if combination_count > maximum_combinations:
        raise RuntimeError(
            "Enumeration requires {} combinations; "
            "maximum is {}".format(
                combination_count,
                maximum_combinations,
            )
        )

    target_inchikey = smiles_inchikey(
        target_smiles
    )
    products = []
    matches = []
    seen_products = set()

    quiet_output = io.StringIO()

    with contextlib.redirect_stdout(quiet_output):
        for combination in itertools.product(
            *ordered_lists
        ):
            ordered_blocks = order_blocks_for_reaction(
                combination,
                reaction,
            )

            if ordered_blocks is None:
                continue

            normalized_blocks = []

            for block in ordered_blocks:
                normalized = dict(block)
                normalized["bb_id"] = (
                    normalized.get("bb_id")
                    or normalized.get("compound_id")
                    or "BB"
                )
                normalized_blocks.append(normalized)

            ordered_blocks = tuple(normalized_blocks)

            product = enum.Determine_and_Run_Reactions(
                ordered_blocks,
                reaction,
                settings,
            )

            if (
                not product
                or isinstance(product, bool)
            ):
                continue

            product_inchikey = smiles_inchikey(
                product.smiles
            )
            block_ids = [
                block.get("compound_id", "")
                for block in ordered_blocks
            ]
            block_smiles = [
                block.get("smiles", "")
                for block in ordered_blocks
            ]
            identity = (
                product_inchikey,
                tuple(block_ids),
            )

            if identity in seen_products:
                continue

            seen_products.add(identity)

            row = {
                "product_inchikey": product_inchikey,
                "product_smiles": product.smiles,
                "reaction_id": reaction_id,
                "tranche": product.tranche,
                "building_block_ids": ",".join(
                    block_ids
                ),
                "building_block_smiles": ",".join(
                    block_smiles
                ),
            }
            products.append(row)

            if (
                product_inchikey
                and product_inchikey
                == target_inchikey
            ):
                matches.append(row)

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    products_path = os.path.join(
        output_directory,
        "made_enumerated_products.txt",
    )
    matches_path = os.path.join(
        output_directory,
        "made_inchikey_matches.txt",
    )

    write_enumerated_products(
        products_path,
        products,
    )
    write_enumerated_products(
        matches_path,
        matches,
    )

    return {
        "status": "complete",
        "combination_count": combination_count,
        "product_count": len(products),
        "match_count": len(matches),
        "products": products,
        "matches": matches,
        "products_path": products_path,
        "matches_path": matches_path,
    }
