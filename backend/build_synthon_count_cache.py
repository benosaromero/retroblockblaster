import json
import re
import subprocess
import sys
from pathlib import Path

from backend.main import COUNT_ROWS, find_synthon_file, read_job_rows

CACHE_PATH = Path("cache/synthon_counts.json")


def load_cache():
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text())
    return {}


def save_cache(cache):
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(json.dumps(cache, indent=2, sort_keys=True))

def count_from_log(path):
    log_path = Path(str(path) + ".smi_log")

    if not log_path.exists():
        return None

    text = log_path.read_text(errors="replace")

    match = re.search(r"Encoded\s+([\d,]+)\s+molecules", text)
    if match:
        return int(match.group(1).replace(",", ""))

    match = re.search(r"Encoded\s+([\d,]+)\s+records", text)
    if match:
        return int(match.group(1).replace(",", ""))

    return None


def line_count(path):
    result = subprocess.check_output(["wc", "-l", str(path)])
    return int(result.decode("utf-8").strip().split()[0])


def sources_from_rows():
    sources = []
    for row in COUNT_ROWS:
        for key in ["cumulative", "incremental"]:
            source = row.get(key)
            if source and source not in sources:
                sources.append(source)
    return sources


def synthon_ids_for_job(job_id):
    ids = set()
    for row in read_job_rows(job_id):
        for key in ["bb1_synthon_id", "bb2_synthon_id", "bb3_synthon_id", "bb4_synthon_id"]:
            value = str(row.get(key, "")).strip()
            if value:
                ids.add(value)
    return sorted(ids)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 backend/build_synthon_count_cache.py JOB_ID")
        sys.exit(1)

    job_id = sys.argv[1]
    cache = load_cache()
    sources = sources_from_rows()

    for synthon_id in synthon_ids_for_job(job_id):
        cache.setdefault(str(synthon_id), {})

        for source in sources:

            if source in cache[str(synthon_id)]:
                continue

            path = find_synthon_file(synthon_id, source)
            if path is None:
                cache[str(synthon_id)][source] = 0
                print("{} {} missing".format(synthon_id, source))
                continue

            count = count_from_log(path)

            if count is None:
                count = line_count(path)

            cache[str(synthon_id)][source] = count
            print("{} {} {}".format(synthon_id, source, count))
            save_cache(cache)

    save_cache(cache)


if __name__ == "__main__":
    main()
