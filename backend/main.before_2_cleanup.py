# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

from urllib.parse import urlencode
from urllib.request import Request, urlopen
from pathlib import Path
import subprocess
import uuid
import zipfile
import os
import smtplib
from email.message import EmailMessage
import json
import fnmatch

from fastapi.staticfiles import StaticFiles
from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel
from backend.analog_pipeline import run_coverage_pipeline
from backend.made_api import router as made_router

load_dotenv(override=True)

app = FastAPI()
app.include_router(made_router)

app.mount("/static", StaticFiles(directory="static"), name="static")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

COUNT_ROWS = [
    {
        "label": "$50",
        "new_source": "incr050",
        "total_source": "db050",
    },
    {
        "label": "$100",
        "new_source": "incr100",
        "total_source": "db100",
    },
    {
        "label": "$250",
        "new_source": "incr250",
        "total_source": "db250",
    },
    {
        "label": "$500",
        "new_source": "incr500",
        "total_source": "db500",
    },
    {
        "label": "$999",
        "new_source": "incr999",
        "total_source": "db999",
    },
    {
        "label": "FAST",
        "new_source": "bbAFAS",
        "total_source": "bbAFAS",
    },
    {
        "label": "MADE",
        "new_source": "dball",
        "total_source": "dball",
    },
]


def send_finished_email(to_email: str, job_id: str):
    download_url = f"http://127.0.0.1:8000/download/{job_id}"

    message = EmailMessage()
    message["Subject"] = "Your BenZinc26 job is finished"
    message["From"] = os.environ["FROM_EMAIL"]
    message["To"] = to_email

    message.set_content(
        f"""Your BenZinc26 job is finished.

Job ID: {job_id}

Download your output.zip here:
{download_url}
"""
    )

    with smtplib.SMTP(os.environ["SMTP_HOST"], int(os.environ["SMTP_PORT"])) as smtp:
        smtp.starttls()
        smtp.login(os.environ["SMTP_USER"], os.environ["SMTP_PASSWORD"])
        smtp.send_message(message)


class HelloJobRequest(BaseModel):
    user_input: str
    email: str
    source: str = "db050"
    tc_cutoff: float = 0.6
    sort_mode: str = "balanced"

class WorkflowRequest(BaseModel):
    run: str
    db: str = "dball"
    input_tsv: str = ""

class CoveragePipelineRequest(BaseModel):
    parent_smiles: str
    modifications: dict
    buyable_analogs: list = []
    makeable_analogs: list = []

def lookup_synthon_building_blocks(synthon_id, limit=50, reference_smiles=""):
    import arthor

    synthons_path = "/nfs/exj/Fe/scripts/common_db/export/Synthons.json"
    db_path = "/mnt/nfs/exj/arthor_db_backup/public_arthor/bb_instock.smi"

    with open(synthons_path, "r") as f:
        synthons = json.load(f)

    synthon = None

    if isinstance(synthons, dict):
        synthon = synthons.get(str(synthon_id))
        if synthon is None and str(synthon_id).isdigit():
            synthon = synthons.get(int(synthon_id))
    elif isinstance(synthons, list):
        for item in synthons:
            if isinstance(item, dict) and str(item.get("id")) == str(synthon_id):
                synthon = item
                break

    if synthon is None:
        return {"error": "synthon not found", "synthon_id": str(synthon_id), "results": []}

    smarts_list = synthon.get("inSMARTS") or synthon.get("smarts") or synthon.get("SMARTS") or synthon.get("query")
    smarts = smarts_list[0] if isinstance(smarts_list, list) else smarts_list

    if not smarts:
        return {
            "error": "synthon has no SMARTS query",
            "synthon_id": str(synthon_id),
            "synthon_name": synthon.get("Name", ""),
            "in_smarts": smarts_list,
            "results": [],
        }

    db = arthor.SubDb("{}.atdb".format(db_path), smifile=db_path)
    qry = arthor.Query(smarts, "Sma:QI")

    ref_fp = None
    if reference_smiles:
        try:
            from rdkit import Chem, DataStructs
            from rdkit.Chem import AllChem

            ref_mol = Chem.MolFromSmiles(reference_smiles)
            if ref_mol is not None:
                ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
        except Exception:
            ref_fp = None

    results = []
    for hit in db.search(qry):
        line = hit.decode("utf-8") if hasattr(hit, "decode") else str(hit)
        parts = line.strip().split()

        smiles = parts[0] if len(parts) > 0 else ""
        bb_id = parts[1] if len(parts) > 1 else ""
        tc_score = ""

        if ref_fp is not None and smiles:
            try:
                hit_mol = Chem.MolFromSmiles(smiles)
                if hit_mol is not None:
                    hit_fp = AllChem.GetMorganFingerprintAsBitVect(hit_mol, 2, nBits=2048)
                    tc_score = round(DataStructs.TanimotoSimilarity(ref_fp, hit_fp), 3)
            except Exception:
                tc_score = ""

        if ref_fp is not None and (tc_score == "" or tc_score < 0.2):
            continue

        results.append({
            "smiles": smiles,
            "bb_id": bb_id,
            "synthon_id": str(synthon_id),
            "synthon_name": synthon.get("Name", ""),
            "db": "bb-stock",
            "tc_score": tc_score,
        })

        if len(results) >= limit:
            break

    if ref_fp is not None:
        results.sort(key=lambda item: item.get("tc_score", 0), reverse=True)

    return {
        "synthon_id": str(synthon_id),
        "synthon_name": synthon.get("Name", ""),
        "in_smarts": smarts_list,
        "db": "bb-stock",
        "results": results,
    }

def run_enumeration(data, mode="manual"):
    reaction_id = data.get("reaction_id", "")
    pair_id = data.get("pair_id", uuid.uuid4().hex[:12])

    if not reaction_id:
        return {"error": "missing reaction_id"}

    if mode == "all":
        block_lists = data.get("block_lists", {})
        if not block_lists:
            return {"error": "missing block_lists"}

        enum_dir = Path("jobs") / pair_id / "enumerate_all"
        enum_dir.mkdir(parents=True, exist_ok=True)

        input_files = []
        total_combos = 1

        for slot in sorted(block_lists.keys()):
            blocks = block_lists.get(slot, [])
            if not blocks:
                return {"error": "missing building blocks for {}".format(slot)}

            total_combos *= len(blocks)

            first_block = blocks[0]
            synthon_id = str(first_block.get("synthon_id", "")).strip()
            synthon_name = str(first_block.get("synthon_name", "Selected")).replace(" ", "-")

            if not synthon_id:
                return {"error": "missing synthon_id for {}".format(slot)}

            path = enum_dir / "web_synthon_{}_{}".format(synthon_id, synthon_name)

            with path.open("w") as f:
                for block in blocks:
                    smiles = str(block.get("smiles", "")).strip()
                    bb_id = str(block.get("bb_id", "BB")).strip()
                    if smiles:
                        f.write("{} {}\n".format(smiles, bb_id))

            input_files.append(str(path))

    else:
        reagents = data.get("reagents", {})
        if not reagents:
            return {"error": "missing reagents"}

        enum_dir = Path("jobs") / pair_id / "enumeration"
        enum_dir.mkdir(parents=True, exist_ok=True)

        input_files = []

        for slot in sorted(reagents.keys()):
            reagent = reagents[slot]
            synthon_id = str(reagent.get("synthon_id", "")).strip()
            synthon_name = str(reagent.get("synthon_name", "Selected")).replace(" ", "-")
            smiles = str(reagent.get("smiles", "")).strip()
            bb_id = str(reagent.get("bb_id", "BB")).strip()

            if not synthon_id or not smiles:
                continue

            path = enum_dir / "web_synthon_{}_{}".format(synthon_id, synthon_name)

            with path.open("w") as f:
                f.write("{} {}\n".format(smiles, bb_id))

            input_files.append(str(path))

        total_combos = len(input_files)

    if len(input_files) < 2:
        return {"error": "need at least two selected building blocks"}

    output_prefix = str(enum_dir / "enum_output")

    cmd = [
        "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python",
        "backend/enumerator_web.py",
        "-i",
    ] + input_files + [
        "-r", reaction_id,
        "-o", output_prefix,
        "-j", "/nfs/exj/Fe/scripts/common_db/export/Synthons.json",
        "-e", "/nfs/exj/Fe/scripts/testreactions2/EXAMPLE_Library_Enumeration_File_hg2_noProtected_logP350_HAC35_RotB7_V3_A-AK1.json",
        "-a", "/nfs/exj/Fe/scripts/common_db/export/Reactions.json",
        "-p", "1",
    ]

    env = os.environ.copy()
    env["PYTHONPATH"] = "/nfs/exj/Fe/scripts"

    subprocess.check_call(cmd, env=env)

    output_path = "{}_{}.out".format(output_prefix, reaction_id)
    products = []

    if os.path.exists(output_path):
        with open(output_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    products.append({
                        "product_code": parts[0],
                        "product_smiles": parts[1],
                        "tranche": parts[2] if len(parts) > 2 else "",
                        "building_blocks": parts[3:] if len(parts) > 3 else [],
                    })

    return {
        "message": "Enumeration finished" if mode != "all" else "Enumerate all finished",
        "reaction_id": reaction_id,
        "pair_id": pair_id,
        "total_combos": total_combos,
        "products": products,
    }

def read_job_rows(job_id):
    import csv

    output_path = Path("jobs") / job_id / "output.txt"

    if not output_path.exists():
        return []

    with output_path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def synthon_ids_for_job(job_id):
    rows = read_job_rows(job_id)
    synthon_ids = set()

    for row in rows:
        for key in ["bb1_synthon_id", "bb2_synthon_id", "bb3_synthon_id", "bb4_synthon_id"]:
            value = str(row.get(key, "")).strip()
            if value:
                synthon_ids.add(value)

    return sorted(synthon_ids)

COUNT_SOURCES = [
    "db050",
    "db100",
    "db250",
    "db500",
    "db999",
    "bbAFAS",
    "dball",
    "incrMADE",
]

def count_file_lines(path):
    if path is None or not path.exists():
        return 0

    output = subprocess.check_output(["wc", "-l", str(path)])
    return int(output.decode("utf-8").strip().split()[0])

def synthon_source_counts(synthon_id):
    source_counts = {}

    for row in COUNT_ROWS:
        for source_key in ["new_source", "total_source"]:
            source = row.get(source_key)

            if not source:
                continue

            if source in source_counts:
                continue

            path = find_synthon_file(synthon_id, source)
            source_counts[source] = count_file_lines(path) if path else 0

    rows = []

    for row in COUNT_ROWS:
        new_source = row.get("new_source")
        total_source = row.get("total_source")

        rows.append({
            "label": row["label"],
            "new_source": new_source,
            "new_count": source_counts.get(new_source, 0),
            "total_source": total_source,
            "total_count": source_counts.get(total_source, 0),
        })

    return rows

SYNTHON_ROOT = Path("/mnt/nfs/exj/Fe/synthons")


def normalize_synthon_id_for_file(synthon_id):
    return str(synthon_id).strip().replace(".", "-")


def synthon_source_dir(source):
    allowed = {
        "db050": SYNTHON_ROOT / "db050",
        "db100": SYNTHON_ROOT / "db100",
        "db250": SYNTHON_ROOT / "db250",
        "db500": SYNTHON_ROOT / "db500",
        "db999": SYNTHON_ROOT / "db999",
        "incr050": SYNTHON_ROOT / "incr050",
        "incr100": SYNTHON_ROOT / "incr100",
        "incr250": SYNTHON_ROOT / "incr250",
        "incr500": SYNTHON_ROOT / "incr500",
        "incr999": SYNTHON_ROOT / "incr999",
        "bbA050": SYNTHON_ROOT / "all" / "bbA050",
        "bbA100": SYNTHON_ROOT / "all" / "bbA100",
        "bbA250": SYNTHON_ROOT / "all" / "bbA250",
        "bbA500": SYNTHON_ROOT / "all" / "bbA500",
        "bbA999": SYNTHON_ROOT / "all" / "bbA999",
        "bbAFAS": SYNTHON_ROOT / "incrFAST",
        "bbAall": SYNTHON_ROOT / "dball",
        "dball": SYNTHON_ROOT / "dball",
        "incrFAST": SYNTHON_ROOT / "incrFAST",
        "incrMADE": SYNTHON_ROOT / "incrMADE",
    }

    return allowed.get(source, allowed["db050"])


def find_synthon_file(synthon_id, source):
    directory = synthon_source_dir(source)

    if not directory.exists():
        return None

    file_synthon_id = normalize_synthon_id_for_file(synthon_id)
    file_source = "incrFAST" if source == "bbAFAS" else source

    prefixes = [
        "{}_{}_".format(file_source, file_synthon_id),
        "{}_synthon_{}_".format(file_source, file_synthon_id),
    ]

    for path in directory.iterdir():
        if not path.is_file():
            continue

        name = path.name

        if name.endswith((".smi", ".smi_log", ".atdb", ".atfp", ".fp")):
            continue

        if any(name.startswith(prefix) for prefix in prefixes):
            return path

    return None


def clean_synthon_name_from_file(path_name, synthon_id):
    name = path_name

    for suffix in [".smi", ".smi_log", ".atdb", ".atfp", ".fp"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]

    file_synthon_id = normalize_synthon_id_for_file(synthon_id)

    markers = [
        "_{}_".format(file_synthon_id),
        "_synthon_{}_".format(file_synthon_id),
    ]

    for marker in markers:
        if marker in name:
            return name.split(marker, 1)[1]

    return name


def read_synthon_file(synthon_id, source="db050", limit=50):
    path = find_synthon_file(synthon_id, source)

    if path is None:
        return {
            "error": "synthon file not found",
            "synthon_id": str(synthon_id),
            "source": source,
            "results": [],
        }

    results = []
    clean_synthon_name = clean_synthon_name_from_file(path.name, synthon_id)

    with path.open("r", errors="replace") as handle:
        for line in handle:
            parts = line.strip().split()

            if len(parts) < 1:
                continue

            smiles = parts[0]
            bb_id = parts[1] if len(parts) > 1 else ""

            results.append({
                "smiles": smiles,
                "bb_id": bb_id,
                "synthon_id": str(synthon_id),
                "synthon_name": clean_synthon_name,
                "db": source,
                "tc_score": "",
            })

            if len(results) >= limit:
                break

    return {
        "synthon_id": str(synthon_id),
        "synthon_name": clean_synthon_name,
        "source": source,
        "file": str(path),
        "results": results,
    }


@app.get("/", response_class=HTMLResponse)
def home():
    return Path("index.html").read_text()

@app.post("/hello-job")
def submit_hello_job(data: HelloJobRequest):
    job_id = uuid.uuid4().hex[:12]
    job_dir = Path("jobs") / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    settings_path = job_dir / "settings.json"
    settings_path.write_text(json.dumps({
        "source": data.source,
        "tc_cutoff": data.tc_cutoff,
        "sort_mode": data.sort_mode,
    }, indent=2))

    input_path = job_dir / "input.smi"
    output_path = job_dir / "output.txt"
    zip_path = job_dir / "output.zip"

    input_path.write_text(data.user_input)

    subprocess.run(
        [
            "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python",
            "run_job.py",
            str(input_path),
            str(output_path),
        ],
        check=True,
    )

    with zipfile.ZipFile(zip_path, "w") as zip_file:
        zip_file.write(output_path, arcname="output.txt")

    email_status = "sent"

    try:
        send_finished_email(data.email, job_id)
    except Exception as error:
        email_status = f"failed: {error}"

    return {
        "job_id": job_id,
        "message": "Job finished",
        "output": output_path.read_text(),
        "email_status": email_status,
    }


@app.post("/reaction")
def reaction_job(data: HelloJobRequest):
    return submit_hello_job(data)


@app.post("/retro")
def retro_job(data: HelloJobRequest):
    return submit_hello_job(data)

@app.get("/health")
def health():
    return {
        "status": "running",
        "service": "Retro-Blockblaster 2.0",
        "port": 8003,
        "main_page": "/workflow",
        "workflow_api": "/api/workflow/full",
    }

@app.get("/synthons/{job_id}", response_class=HTMLResponse)
def synthons(job_id: str):
    output_path = Path("jobs") / job_id / "output.txt"

    if not output_path.exists():
        return Path("index.html").read_text()

    return f"""
<!DOCTYPE html>
<html>
<head>
  <title>Retro-Blockblaster Synthons</title>
  <link rel="stylesheet" href="/static/styles.css">
</head>
<body data-job-id="{job_id}">
  <div class="page">
    <div class="container">
      <header class="header">
        <h1>Retro-Blockblaster 1.0</h1>
        <p>Synthons and retrosynthetic route results | Job ID: {job_id}</p>
      </header>

      <section class="panel">
        <h2>Synthon Results</h2>
        <p id="result">Loading job results...</p>
      </section>

      <footer class="footer">
        <p>
          <a class="button-link" href="/">Submit another job</a>
          <a class="button-link" href="/download/{job_id}">Download output.zip</a>
        </p>
      </footer>
    </div>
  </div>

  <script src="/static/app.js"></script>
</body>
</html>
"""

@app.get("/api/synthons/{job_id}")
def synthons_api(job_id: str):
    output_path = Path("jobs") / job_id / "output.txt"

    if not output_path.exists():
        return {"error": "Job output not found", "job_id": job_id}

    return {
        "job_id": job_id,
        "output": output_path.read_text(),
        "download": f"/download/{job_id}",
    }


@app.get("/image/{job_id}/{filename}")
def image(job_id: str, filename: str):
    image_path = Path("jobs") / job_id / "images" / filename

    if not image_path.exists():
        return {"error": "image not found"}

    return FileResponse(image_path, media_type="image/svg+xml")

BETA_BACKEND = "http://127.0.0.1:8001"


def proxy_beta_get(path):
    with urlopen(BETA_BACKEND + path, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def proxy_beta_post(path, payload):
    from urllib.error import HTTPError, URLError

    url = "http://127.0.0.1:8001{}".format(path)
    body = json.dumps(payload).encode("utf-8")

    request = Request(
        url,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    try:
        with urlopen(request, timeout=300) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        error_body = error.read().decode("utf-8", errors="replace")
        return {
            "error": "beta endpoint returned HTTP {}".format(error.code),
            "beta_path": path,
            "beta_body": error_body,
        }
    except URLError as error:
        return {
            "error": "could not connect to beta endpoint",
            "beta_path": path,
            "details": str(error),
        }

def add_tc_scores(reference_smiles, results):
    if not reference_smiles or not results:
        return results

    payload = json.dumps({
        "reference_smiles": reference_smiles,
        "candidates": results,
    }).encode("utf-8")

    process = subprocess.Popen(
        ["/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python", "backend/tc_score_cli.py"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    stdout, stderr = process.communicate(payload)

    if process.returncode != 0:
        return results

    try:
        return json.loads(stdout.decode("utf-8"))
    except Exception:
        return results


def safe_workflow_name(value):
    run = str(value).strip()

    if not run:
        raise HTTPException(status_code=400, detail="missing run name")

    if run.endswith(".tsv"):
        run = run[:-4]

    allowed = []
    for char in run:
        if char.isalnum() or char in "-_":
            allowed.append(char)

    cleaned = "".join(allowed)

    if cleaned != run:
        raise HTTPException(
            status_code=400,
            detail="Run name can only contain letters, numbers, dash, and underscore",
        )

    return cleaned


def run_workflow_command(command):
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    if completed.returncode != 0:
        raise HTTPException(
            status_code=500,
            detail=completed.stderr or completed.stdout or "Workflow command failed",
        )

    return completed.stdout.strip()


@app.post("/api/workflow/input")
def workflow_input(data: WorkflowRequest):
    run = safe_workflow_name(data.run)

    if not data.input_tsv.strip():
        raise HTTPException(status_code=400, detail="missing input TSV")

    Path(run + ".tsv").write_text(data.input_tsv)

    return {
        "run": run,
        "step": "input",
        "message": "Saved {}.tsv".format(run),
    }


@app.post("/api/workflow/retro")
def workflow_retro(data: WorkflowRequest):
    run = safe_workflow_name(data.run)
    output = run_workflow_command(["python3", "retro.py", run])
    return {"run": run, "step": "retro", "message": output}


@app.post("/api/workflow/bbs")
def workflow_bbs(data: WorkflowRequest):
    run = safe_workflow_name(data.run)
    output = run_workflow_command(["python3", "bbs.py", "-db", data.db, run])
    return {"run": run, "step": "bbs", "message": output}


@app.post("/api/workflow/make")
def workflow_make(data: WorkflowRequest):
    run = safe_workflow_name(data.run)
    output = run_workflow_command(["python3", "make.py", run])
    return {"run": run, "step": "make", "message": output}


@app.post("/api/workflow/rank")
def workflow_rank(data: WorkflowRequest):
    run = safe_workflow_name(data.run)
    output = run_workflow_command(["python3", "rank.py", run])
    return {"run": run, "step": "rank", "message": output}


@app.post("/api/workflow/full")
def workflow_full(data: WorkflowRequest):
    run = safe_workflow_name(data.run)
    output = run_workflow_command(["python3", "workflow.py", "-db", data.db, run])
    return {"run": run, "step": "full", "message": output}

@app.get("/synthon/{synthon_id}")
def synthon_lookup(synthon_id: str, source: str = "db050", limit: int = 50, reference_smiles: str = "", tc_cutoff: float = 0.6, fallback_top_match: bool = True):
    result = read_synthon_file(
        synthon_id=synthon_id,
        source=source,
        limit=limit,
    )

    scored_results = add_tc_scores(
        reference_smiles=reference_smiles,
        results=result.get("results", []),
    )

    if reference_smiles:
        scored_results.sort(
            key=lambda item: item.get("tc_score") if item.get("tc_score") != "" else -1,
            reverse=True,
        )

        passing_results = [
            item for item in scored_results
            if item.get("tc_score") != "" and item.get("tc_score", 0) >= tc_cutoff
        ]

        if passing_results:
            result["results"] = passing_results
            result["message"] = "{} building blocks passed Tc cutoff {}.".format(
                len(passing_results),
                tc_cutoff,
            )
        else:
            if fallback_top_match:
                result["results"] = scored_results[:1]
                result["message"] = "No building blocks were over Tc cutoff {}. Showing highest Tc match.".format(
                    tc_cutoff,
                )
            else:
                result["results"] = []
                result["message"] = "No building blocks were over Tc cutoff {}. Lower Tc cutoff or allow top-match fallback.".format(
                    tc_cutoff,
                )

    else:
        result["results"] = scored_results
        result["message"] = "No reference SMILES provided; Tc cutoff was not applied."

    result["tc_cutoff"] = tc_cutoff
    result["shown_count"] = len(result["results"])

    return result

@app.post("/enumerate")
def enumerate_selected(data: dict):
    return proxy_beta_post("/enumerate", data)


@app.post("/enumerate-all")
def enumerate_all(data: dict):
    return proxy_beta_post("/enumerate-all", data)

@app.post("/api/preload-synthons/{job_id}")
def preload_synthons(job_id: str):
    synthon_ids = synthon_ids_for_job(job_id)
    loaded = []
    failed = []

    for synthon_id in synthon_ids:
        try:
            result = proxy_beta_get("/synthon/{}?reference_smiles=".format(synthon_id))

            loaded.append({
                "synthon_id": synthon_id,
                "count": len(result.get("results", [])),
            })
        except Exception as error:
            failed.append({
                "synthon_id": synthon_id,
                "error": str(error),
            })

    return {
        "job_id": job_id,
        "synthon_count": len(synthon_ids),
        "loaded": loaded,
        "failed": failed,
    }

@app.get("/api/reaction/{job_id}/{reaction_id}")
def reaction_detail_api(job_id: str, reaction_id: str):
    rows = [
        row for row in read_job_rows(job_id)
        if row.get("reaction_id") == reaction_id
    ]

    if not rows:
        return {
            "error": "Reaction not found",
            "job_id": job_id,
            "reaction_id": reaction_id,
        }

    return {
        "job_id": job_id,
        "reaction_id": reaction_id,
        "rows": rows,
    }


@app.get("/reaction/{job_id}/{reaction_id}", response_class=HTMLResponse)
def reaction_detail_page(job_id: str, reaction_id: str):
    return f"""
<!DOCTYPE html>
<html>
<head>
  <title>Retro-Blockblaster {reaction_id}</title>
  <link rel="stylesheet" href="/static/styles.css">
</head>
<body data-job-id="{job_id}" data-reaction-id="{reaction_id}">
  <div class="page">
    <div class="container">
      <header class="header">
        <h1>{reaction_id}</h1>
        <p>Reaction detail | Job ID: {job_id}</p>
      </header>

      <section class="panel">
        <h2>Reaction Detail</h2>
        <p id="result">Loading reaction detail...</p>
      </section>

      <footer class="footer">
        <p>
          <a class="button-link" href="/synthons/{job_id}">Back to synthons</a>
          <a class="button-link" href="/download/{job_id}">Download output.zip</a>
        </p>
      </footer>
    </div>
  </div>

  <script src="/static/app.js"></script>
</body>
</html>
"""

@app.get("/draw")
def draw_molecule(smiles: str = ""):
    if not smiles:
        return {"error": "missing smiles"}

    result = subprocess.run(
        [
            "/nfs/exj/Fe/scripts/common_db/arthor-env/bin/python",
            "backend/draw_molecule_cli.py",
            smiles,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    if result.returncode != 0:
        raise HTTPException(status_code=500, detail=result.stderr)

    return HTMLResponse(result.stdout, media_type="image/svg+xml")

def load_synthon_count_cache():
    path = Path("cache/synthon_counts.json")

    if not path.exists():
        return {}

    return json.loads(path.read_text())


def synthon_ids_for_job(job_id):
    ids = set()

    for row in read_job_rows(job_id):
        for key in ["bb1_synthon_id", "bb2_synthon_id", "bb3_synthon_id", "bb4_synthon_id"]:
            synthon_id = str(row.get(key, "")).strip()

            if synthon_id:
                ids.add(synthon_id)

    return sorted(ids)

@app.get("/coverage", response_class=HTMLResponse)
def coverage_page():
    return """
<!DOCTYPE html>
<html>
<head>
  <title>Retro-Blockblaster 1.0 | Analog Coverage</title>
  <link rel="stylesheet" href="/static/styles.css">
</head>
<body>
  <div class="page">
    <div class="container">
      <header class="header">
        <h1>Retro-Blockblaster 1.0</h1>
        <p>Analog coverage across buyable, makeable, overlapping, and inaccessible chemical space.</p>
      </header>

      <section class="panel coverage-panel">
        <h2>Analog Coverage</h2>

        <div class="coverage-instructions">
          <p>
            Enter a parent SMILES and define atom-level modifications. The tool generates close analogs, then compares them against optional buyable and makeable analog lists.
          </p>
          <p>
            Use the submitted Retro-Blockblaster job to find makeable analogs, then paste product SMILES into the makeable list if needed.
          </p>
        </div>

        <form id="coverage-form">
          <label for="coverage-parent-smiles">Parent SMILES</label>
          <textarea id="coverage-parent-smiles" class="coverage-textarea is-compact" required>CCOC(=O)C1=CN=CN1[C@H](C)C2=CC=CC=C2</textarea>

          <div class="coverage-action-row">
            <button type="button" class="coverage-secondary-button" onclick="openRetroForCoverage()">
              Open Submitted Retro-Blockblaster Job
            </button>
          </div>

          <label for="coverage-modifications">Modifications JSON</label>
          <textarea id="coverage-modifications" class="coverage-textarea" required>{
  "0": ["O", "F", "C"],
  "6": ["C"],
  "14": ["C"],
  "15": ["N", "O"],
  "16": ["O"]
}</textarea>

          <div class="coverage-optional-controls">
            <button type="button" class="coverage-secondary-button" onclick="toggleCoverageOptional('coverage-buyable-wrap', this)">
              Add buyable analog SMILES
            </button>

            <button type="button" class="coverage-secondary-button" onclick="toggleCoverageOptional('coverage-makeable-wrap', this)">
              Add makeable analog SMILES
            </button>
          </div>

          <div id="coverage-buyable-wrap" class="coverage-optional-section">
            <label for="coverage-buyable">Buyable analog SMILES</label>
            <textarea id="coverage-buyable" class="coverage-textarea" placeholder="One SMILES per line"></textarea>
          </div>

          <div id="coverage-makeable-wrap" class="coverage-optional-section">
            <label for="coverage-makeable">Makeable analog SMILES</label>
            <textarea id="coverage-makeable" class="coverage-textarea" placeholder="One SMILES per line"></textarea>
          </div>

          <button type="submit">Run Coverage</button>
        </form>
      </section>

      <section class="panel">
        <h2>Coverage Result</h2>
        <div id="coverage-result">No coverage job submitted yet.</div>
      </section>

      <footer class="footer">
        <p>
          Irwin Lab, UCSF | Retro-Blockblaster 1.0
        </p>
      </footer>
    </div>
  </div>

  <script src="/static/coverage.js"></script>
</body>
</html>
"""

@app.get("/api/synthon-counts/{job_id}")
def synthon_counts_api(job_id: str):
    cache = load_synthon_count_cache()
    counts = {}

    for synthon_id in synthon_ids_for_job(job_id):
        source_counts = cache.get(str(synthon_id), {})
        rows = []

        for row in COUNT_ROWS:
            new_source = row.get("new_source")
            total_source = row.get("total_source")

            rows.append({
                "label": row["label"],
                "new_source": new_source,
                "new_count": source_counts.get(new_source) if new_source else 0,
                "total_source": total_source,
                "total_count": source_counts.get(total_source) if total_source else 0,
            })

        counts[str(synthon_id)] = rows

    return {
        "job_id": job_id,
        "counts": counts,
        "rows": COUNT_ROWS,
    }

@app.get("/api/job-settings/{job_id}")
def job_settings(job_id: str):
    settings_path = Path("jobs") / job_id / "settings.json"

    if not settings_path.exists():
        return {
            "source": "db050",
            "tc_cutoff": 0.6,
            "sort_mode": "balanced",
        }

    return json.loads(settings_path.read_text())

@app.post("/api/coverage-pipeline")
def coverage_pipeline(data: CoveragePipelineRequest):
    return run_coverage_pipeline(
        parent_smiles=data.parent_smiles,
        modifications_dict=data.modifications,
        buyable_analogs=data.buyable_analogs,
        makeable_analogs=data.makeable_analogs,
    )

@app.get("/coverage-file/{job_id}/{filename}")
def coverage_file(job_id: str, filename: str):
    allowed = {
        "coverage.csv",
        "close_analogs.csv",
        "close_analogs.png",
        "close_analogs_parent_numbered.png",
    }

    if filename not in allowed:
        return {"error": "file not allowed"}

    path = Path("jobs") / job_id / "coverage" / filename

    if not path.exists():
        return {"error": "file not found"}

    if filename.endswith(".png"):
        return FileResponse(path, media_type="image/png", filename=filename)

    return FileResponse(path, media_type="text/csv", filename=filename)

@app.get("/workflow", response_class=HTMLResponse)
def workflow_dashboard():
    return Path("workflow.html").read_text()

@app.get("/download/{job_id}")
def download_output(job_id: str):
    zip_path = Path("jobs") / job_id / "output.zip"

    if not zip_path.exists():
        return {"error": "output.zip not found"}

    return FileResponse(
        zip_path,
        media_type="application/zip",
        filename="output.zip",
    )
