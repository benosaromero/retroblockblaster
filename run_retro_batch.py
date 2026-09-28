# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

import csv
import os
import subprocess
import sys
import uuid

RETRO_PYTHON = "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python"

# Batch driver for the full terminal workflow.
#
# This script starts from a company-provided TSV containing molecule IDs and
# SMILES strings. For each molecule, it runs run_job.py, saves the per-molecule
# retrosynthesis output, and combines all successful route/synthon rows into
# all_synthons.tsv.
#
# Expected input columns:
#   id or input_id
#   smiles or SMILES
#
# Usage:
#   python3 run_retro_batch.py molecules.tsv
#   python3 run_retro_batch.py molecules.tsv all_synthons.tsv cli_jobs

INPUT_TSV = sys.argv[1] if len(sys.argv) > 1 else "molecules.tsv"
OUT_TSV = sys.argv[2] if len(sys.argv) > 2 else "all_synthons.tsv"
JOB_DIR = sys.argv[3] if len(sys.argv) > 3 else "cli_jobs"

JOB_IDS_TSV = "job_ids.tsv"
JOB_SUMMARY_TSV = "job_summary.txt"

os.makedirs(JOB_DIR, exist_ok=True)


def safe_name(text):
    """Convert arbitrary IDs into filenames that are safe to use in shell paths."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(text))


def get_first_value(row, names):
    """Return the first non-empty value from a list of possible column names."""
    for name in names:
        value = row.get(name, "")
        if value:
            return value.strip()
    return ""


def read_input_molecules(path):
    """Read query molecules from a company TSV."""
    molecules = []

    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for index, row in enumerate(reader, start=1):
            input_id = get_first_value(row, ["input_id", "id", "ID", "name", "compound_id"])
            smiles = get_first_value(row, ["smiles", "SMILES", "canonical_smiles"])

            if not input_id:
                input_id = "mol_{:04d}".format(index)

            if not smiles:
                molecules.append({
                    "input_id": input_id,
                    "smiles": "",
                    "status": "failed",
                    "message": "missing SMILES",
                })
                continue

            molecules.append({
                "input_id": input_id,
                "smiles": smiles,
                "status": "pending",
                "message": "",
            })

    return molecules


def run_one_molecule(molecule):
    """Run run_job.py for one molecule and return job metadata."""
    input_id = molecule["input_id"]
    smiles = molecule["smiles"]

    job_id = uuid.uuid4().hex[:12]
    job_path = os.path.join(JOB_DIR, safe_name(input_id) + "_" + job_id)
    os.makedirs(job_path, exist_ok=True)

    input_path = os.path.join(job_path, "input.smi")
    output_path = os.path.join(job_path, "output.txt")

    with open(input_path, "w") as out:
        out.write("{} {}\n".format(smiles, input_id))

    try:
        subprocess.check_call([
            RETRO_PYTHON,
            "run_job.py",
            input_path,
            output_path,
        ])

        route_count = count_output_routes(output_path)
        status = "complete"
        message = "ok"

        if route_count == 0:
            status = "no_routes"
            message = "no matching retro reactions"

    except Exception as error:
        route_count = 0
        status = "failed"
        message = str(error)

    return {
        "input_id": input_id,
        "smiles": smiles,
        "job_id": job_id,
        "job_path": job_path,
        "input_path": input_path,
        "output_path": output_path,
        "status": status,
        "message": message,
        "route_count": route_count,
    }


def count_output_routes(output_path):
    """Count non-comment data rows in one run_job.py output file."""
    if not os.path.isfile(output_path):
        return 0

    with open(output_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)

    if len(rows) < 2:
        return 0

    return sum(1 for row in rows[1:] if row and not row[0].startswith("#"))


def write_job_tables(job_rows):
    """Write job ID and summary files for the batch run."""
    with open(JOB_IDS_TSV, "w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=[
            "input_id",
            "job_id",
            "job_path",
            "output_path",
        ])
        writer.writeheader()

        for row in job_rows:
            writer.writerow({
                "input_id": row["input_id"],
                "job_id": row["job_id"],
                "job_path": row["job_path"],
                "output_path": row["output_path"],
            })

    with open(JOB_SUMMARY_TSV, "w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=[
            "input_id",
            "smiles",
            "job_id",
            "status",
            "route_count",
            "message",
        ])
        writer.writeheader()

        for row in job_rows:
            writer.writerow({
                "input_id": row["input_id"],
                "smiles": row["smiles"],
                "job_id": row["job_id"],
                "status": row["status"],
                "route_count": row["route_count"],
                "message": row["message"],
            })


def write_all_synthons(job_rows):
    """Combine all per-molecule retrosynthesis outputs into one TSV."""
    wrote_header = False

    with open(OUT_TSV, "w") as combined:
        writer = None

        for job in job_rows:
            output_path = job["output_path"]

            if not os.path.isfile(output_path):
                continue

            with open(output_path) as handle:
                reader = csv.DictReader(
                    (line for line in handle if not line.startswith("#")),
                    delimiter="\t",
                )

                if not reader.fieldnames:
                    continue

                fieldnames = ["input_id", "job_id"] + reader.fieldnames

                if not wrote_header:
                    writer = csv.DictWriter(combined, delimiter="\t", fieldnames=fieldnames)
                    writer.writeheader()
                    wrote_header = True

                for row in reader:
                    if not row.get("reaction_id"):
                        continue

                    output_row = {
                        "input_id": job["input_id"],
                        "job_id": job["job_id"],
                    }
                    output_row.update(row)
                    writer.writerow(output_row)


def main():
    molecules = read_input_molecules(INPUT_TSV)

    if not molecules:
        raise SystemExit("No molecules found in {}".format(INPUT_TSV))

    job_rows = []

    for index, molecule in enumerate(molecules, start=1):
        if not molecule["smiles"]:
            print("Skipping", molecule["input_id"], molecule["message"])
            job_rows.append({
                "input_id": molecule["input_id"],
                "smiles": molecule["smiles"],
                "job_id": "",
                "job_path": "",
                "input_path": "",
                "output_path": "",
                "status": molecule["status"],
                "message": molecule["message"],
                "route_count": 0,
            })
            continue

        print("Running", index, "of", len(molecules), molecule["input_id"])
        job_rows.append(run_one_molecule(molecule))

    write_job_tables(job_rows)
    write_all_synthons(job_rows)

    print("Wrote", JOB_IDS_TSV)
    print("Wrote", JOB_SUMMARY_TSV)
    print("Wrote", OUT_TSV)
    print("Done.")


if __name__ == "__main__":
    main()
