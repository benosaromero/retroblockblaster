# Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
# Developed in the Irwin Lab at UCSF.
# Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.

from __future__ import print_function

import json
import os
import re
import subprocess
import uuid

from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel


# FastAPI router for running the MADE analog workflow as background jobs.
#
# The API keeps job state in jobs/made/<job_id>/state.json and exposes:
#   POST /api/made/jobs
#   GET  /api/made/jobs/{job_id}
#   GET  /api/made/jobs/{job_id}/files/{file_key}
#
# The actual chemistry workflow runs in backend/made_workflow_cli.py under the
# Arthor/RDKit Python environment so the web server does not need RDKit loaded.

router = APIRouter(prefix="/api/made", tags=["MADE"])

PROJECT_ROOT = "/nfs/home/bromero/retro-blockblaster_1.0"
JOB_ROOT = os.path.join(
    PROJECT_ROOT,
    "jobs",
    "made",
)

RESULT_FILES = {
    "made_all": "05_made_bb_candidates_all.txt",
    "made_filtered": "06_made_bb_candidates_filtered.txt",
    "made_products": "07_made_enumerated_products.txt",
    "made_matches": "08_made_inchikey_matches.txt",
}


class MadeJobRequest(BaseModel):
    smiles: str
    filter_mode: str = "one_step"


def validated_job_id(job_id):
    """Accept only UUID-hex job IDs before building filesystem paths."""
    if not re.match(r"^[a-f0-9]{32}$", job_id):
        raise HTTPException(
            status_code=404,
            detail="MADE job not found",
        )

    return job_id


def job_directory(job_id):
    """Return the on-disk folder for a MADE job."""
    return os.path.join(
        JOB_ROOT,
        validated_job_id(job_id),
    )


def state_path(job_id):
    """Return the JSON state file path for a MADE job."""
    return os.path.join(
        job_directory(job_id),
        "state.json",
    )


def write_state(job_id, state):
    """Atomically write job state so status reads never see partial JSON."""
    directory = job_directory(job_id)

    if not os.path.isdir(directory):
        os.makedirs(directory)

    path = state_path(job_id)
    temporary_path = path + ".tmp"

    with open(temporary_path, "w") as output:
        json.dump(state, output, indent=2)

    os.replace(temporary_path, path)


def read_state(job_id):
    """Read current job state or return a FastAPI 404."""
    path = state_path(job_id)

    if not os.path.isfile(path):
        raise HTTPException(
            status_code=404,
            detail="MADE job not found",
        )

    with open(path, "r") as source:
        return json.load(source)


def execute_made_job(
    job_id,
    smiles,
    filter_mode,
):
    """Background task that runs the MADE workflow subprocess."""
    state = read_state(job_id)
    state["status"] = "running"
    state["error"] = None
    write_state(job_id, state)

    try:
        directory = job_directory(job_id)
        result_path = os.path.join(
            directory,
            "workflow_result.json",
        )

        subprocess.check_call(
            [
                "/nfs/exj/Fe/scripts/common_db/"
                "arthor-env/bin/python",
                os.path.join(
                    PROJECT_ROOT,
                    "backend",
                    "made_workflow_cli.py",
                ),
                smiles,
                filter_mode,
                directory,
                result_path,
            ],
            cwd=PROJECT_ROOT,
        )

        with open(result_path, "r") as source:
            result = json.load(source)

        state["status"] = "complete"
        state["result"] = result
        state["files"] = {
            key: (
                "/api/made/jobs/{}/files/{}"
                .format(job_id, key)
            )
            for key in RESULT_FILES
        }
    except Exception as error:
        state["status"] = "failed"
        state["error"] = str(error)

    write_state(job_id, state)


@router.post("/jobs")
def submit_made_job(
    request: MadeJobRequest,
    background_tasks: BackgroundTasks,
):
    """Queue a MADE workflow job and return a polling URL."""
    smiles = request.smiles.strip()
    filter_mode = request.filter_mode.strip()

    if not smiles:
        raise HTTPException(
            status_code=422,
            detail="SMILES is required",
        )

    if filter_mode not in ("one_step", "two_step"):
        raise HTTPException(
            status_code=422,
            detail=(
                "filter_mode must be one_step or two_step; "
                "anon_dist is unavailable because SmallWorld "
                "does not expose incrMADE"
            ),
        )

    job_id = uuid.uuid4().hex

    write_state(job_id, {
        "job_id": job_id,
        "status": "queued",
        "smiles": smiles,
        "filter_mode": filter_mode,
        "result": None,
        "files": {},
        "error": None,
    })

    background_tasks.add_task(
        execute_made_job,
        job_id,
        smiles,
        filter_mode,
    )

    return {
        "job_id": job_id,
        "status": "queued",
        "status_url": (
            "/api/made/jobs/{}".format(job_id)
        ),
    }


@router.get("/jobs/{job_id}")
def get_made_job(job_id: str):
    """Return queued/running/complete/failed job state."""
    return read_state(job_id)


@router.get("/jobs/{job_id}/files/{file_key}")
def download_made_file(
    job_id: str,
    file_key: str,
):
    """Download one result file from a completed MADE job."""
    read_state(job_id)

    filename = RESULT_FILES.get(file_key)

    if filename is None:
        raise HTTPException(
            status_code=404,
            detail="Unknown MADE result file",
        )

    path = os.path.join(
        job_directory(job_id),
        filename,
    )

    if not os.path.isfile(path):
        raise HTTPException(
            status_code=404,
            detail="MADE result file is not available",
        )

    return FileResponse(
        path,
        media_type="text/plain",
        filename=filename,
    )
