import os
import subprocess
import sys
from pathlib import Path

RDKIT_PYTHON = "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python"

def run_name(value):
    return value[:-4] if value.endswith(".tsv") else value

def combine_products(run):
    enum_dir = run + "_enum"
    paths = sorted(Path(enum_dir).glob("*/*/enumerated_products.tsv"))
    all_tsv = Path(run + "_all.tsv")

    if not paths:
        raise SystemExit("No enumerated product files found in " + enum_dir)

    wrote_header = False

    with all_tsv.open("w") as out:
        for path in paths:
            input_id = path.parts[1]
            route_key = path.parts[2]

            with path.open() as handle:
                header = handle.readline().rstrip("\n").split("\t")

                if not wrote_header:
                    out.write("input_id\troute_key\t" + "\t".join(header) + "\n")
                    wrote_header = True

                for line in handle:
                    if line.strip():
                        out.write(input_id + "\t" + route_key + "\t" + line)

    print("Wrote", all_tsv, "from", len(paths), "route files")
    return str(all_tsv)

def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 rank.py RUN_NAME")

    run = run_name(sys.argv[1])
    input_tsv = run + ".tsv"
    all_tsv = combine_products(run)
    rank_tsv = run + "_rank.tsv"
    ranked_dir = run + "_ranked"

    if not os.path.isfile(input_tsv):
        raise SystemExit("Missing input file: " + input_tsv)

    subprocess.check_call([
        RDKIT_PYTHON,
        "rank_products_by_similarity.py",
        input_tsv,
        all_tsv,
        rank_tsv,
        ranked_dir,
    ])

    print("Wrote", rank_tsv)
    print("Wrote", ranked_dir)

if __name__ == "__main__":
    main()
