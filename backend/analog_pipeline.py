import csv
import json
import subprocess
import uuid
from pathlib import Path


ANALOG_PYTHON = "python"
ANALOG_ENV_SOURCE = "/nfs/home/omailhot/pyenv_source.sh"

def canonical_smiles(smiles):
    smiles = str(smiles).strip()

    if not smiles:
        return ""

    command = [
        "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python",
        "backend/canonicalize_smiles_cli.py",
        smiles,
    ]

    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
    except Exception:
        return ""

    return completed.stdout.strip()


def classify_analogs(all_analogs, buyable_analogs, makeable_analogs):
    all_set = set(canonical_smiles(smiles) for smiles in all_analogs)
    buyable_set = set(canonical_smiles(smiles) for smiles in buyable_analogs)
    makeable_set = set(canonical_smiles(smiles) for smiles in makeable_analogs)

    all_set.discard("")
    buyable_set.discard("")
    makeable_set.discard("")

    rows = []

    for smiles in sorted(all_set):
        buyable = smiles in buyable_set
        makeable = smiles in makeable_set

        if buyable and makeable:
            status = "both"
        elif buyable:
            status = "buyable_only"
        elif makeable:
            status = "makeable_only"
        else:
            status = "neither"

        rows.append({
            "canonical_smiles": smiles,
            "buyable": buyable,
            "makeable": makeable,
            "status": status,
        })

    summary = {
        "total_close_analogs": len(all_set),
        "buyable": len(all_set & buyable_set),
        "makeable": len(all_set & makeable_set),
        "both": len(all_set & buyable_set & makeable_set),
        "buyable_only": len((all_set & buyable_set) - makeable_set),
        "makeable_only": len((all_set & makeable_set) - buyable_set),
        "neither": len(all_set - buyable_set - makeable_set),
    }

    return summary, rows


def write_coverage_csv(path, rows):
    with open(path, "w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["canonical_smiles", "buyable", "makeable", "status"],
        )
        writer.writeheader()
        writer.writerows(rows)


def run_close_analog_generator(parent_smiles, modifications_dict, job_dir):
    modifications_path = job_dir / "modifications.json"
    output_prefix = str(job_dir / "close_analogs")

    with open(modifications_path, "w") as handle:
        json.dump(modifications_dict, handle)

    command = (
        "source {} && python backend/analog_generator_cli.py '{}' '{}' '{}'"
        .format(
            ANALOG_ENV_SOURCE,
            parent_smiles.replace("'", ""),
            str(modifications_path),
            output_prefix,
        )
    )

    completed = subprocess.run(
        ["bash", "-lc", command],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    if completed.returncode != 0:
        return {
            "error": "analog generator failed",
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "command": command,
            "analogs": [],
        }

    for line in reversed(completed.stdout.splitlines()):
        line = line.strip()

        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)

    return {
        "error": "analog generator did not return JSON",
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "command": command,
        "analogs": [],
    }


def run_coverage_pipeline(parent_smiles, modifications_dict, buyable_analogs, makeable_analogs):
    job_id = uuid.uuid4().hex[:12]
    job_dir = Path("jobs") / job_id / "coverage"
    job_dir.mkdir(parents=True, exist_ok=True)

    analog_result = run_close_analog_generator(
        parent_smiles=parent_smiles,
        modifications_dict=modifications_dict,
        job_dir=job_dir,
    )

    if analog_result.get("error"):
        return {
            "job_id": job_id,
            "parent_smiles": parent_smiles,
            "analog_generation": analog_result,
            "summary": {},
            "rows": [],
            "coverage_csv": "",
        }

    all_analogs = analog_result.get("analogs", [])

    summary, rows = classify_analogs(
        all_analogs=all_analogs,
        buyable_analogs=buyable_analogs,
        makeable_analogs=makeable_analogs,
    )

    coverage_csv = job_dir / "coverage.csv"
    write_coverage_csv(coverage_csv, rows)

    return {
        "job_id": job_id,
        "parent_smiles": parent_smiles,
        "analog_generation": analog_result,
        "summary": summary,
        "rows": rows[:500],
        "coverage_csv": str(coverage_csv),
        "downloads": {
            "coverage_csv": "/coverage-file/{}/coverage.csv".format(job_id),
            "close_analogs_csv": "/coverage-file/{}/close_analogs.csv".format(job_id),
            "close_analogs_png": "/coverage-file/{}/close_analogs.png".format(job_id),
            "parent_numbered_png": "/coverage-file/{}/close_analogs_parent_numbered.png".format(job_id),
        },
    }
