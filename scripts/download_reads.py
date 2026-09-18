#!/usr/bin/env python3
"""Sequential SRA download; optional validated FASTQ conversion. Run from repo root."""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.validation.reads import validate_reads
from src.validation.validate import sha256


def archive_directory(base, run):
    candidates = [base / run, base / run / run]
    occupied = [p for p in candidates if any(p.glob(f"{run}.sra*"))]
    if len(occupied) > 1:
        raise ValueError(f"Ambiguous archive locations for {run}: {occupied}")
    return occupied[0] if occupied else candidates[0]


def invoke(command, log):
    print(shlex.join(map(str, command)), flush=True)
    with log.open("a") as handle:
        handle.write(f"\n{datetime.now(timezone.utc).isoformat()} {shlex.join(map(str, command))}\n")
        handle.flush()
        child = subprocess.Popen(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        try:
            code = child.wait()
        except KeyboardInterrupt:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            raise
    if code:
        raise RuntimeError(f"Command exited {code}; stopped entire batch. See {log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accessions", default="metadata/source/SRR_Acc_List.txt")
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--run", help="Process only this run from the accession list")
    parser.add_argument("--sra-dir", default="data/sra")
    parser.add_argument("--log-dir", default="logs/download")
    parser.add_argument("--convert", action="store_true", help="After download, validate and convert to gzip FASTQ")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true", help="Print plan without touching data or running tools")
    args = parser.parse_args()
    runs = Path(args.accessions).read_text().splitlines()
    if not runs or len(runs) != len(set(runs)) or any(not re.fullmatch(r"SRR\d+", r) for r in runs):
        raise ValueError("Accession list must contain unique SRR IDs, one per line")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.run:
        if args.run not in runs:
            raise ValueError("--run must be in the accession list")
        runs = [args.run]
    with open(args.samples) as handle:
        samples = {r["run_accession"]: r for r in csv.DictReader(handle, delimiter="\t")}
    base, logs = Path(args.sra_dir), Path(args.log_dir)
    plan = [(run, archive_directory(base, run)) for run in runs]
    if args.dry_run:
        for run, directory in plan:
            print(f"{run}: prefetch -O {directory.parent}; archive={directory}; convert={args.convert}")
        return 0
    for tool in (["prefetch", "vdb-validate", "fasterq-dump"] if args.convert else ["prefetch"]):
        if not shutil.which(tool):
            raise ValueError(f"Missing command: {tool}")
    base.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    with (base / ".e2t-download.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("Another E2T download script is running") from error
        for run, directory in plan:
            completed = logs / f"{run}.completed.json"
            if args.convert and completed.is_file():
                previous = json.loads(completed.read_text())
                expected = [samples[run][k] for k in ("fastq_1", "fastq_2") if samples[run][k]]
                artifacts = previous["reads"]["files"]
                if previous.get("status") != "converted" or [a["path"] for a in artifacts] != expected:
                    raise ValueError(f"Completed conversion manifest conflicts with samples.tsv: {run}")
                for artifact in artifacts:
                    path = Path(artifact["path"])
                    if not path.is_file() or sha256(path) != artifact["sha256"]:
                        raise ValueError(f"Completed FASTQ missing or modified: {path}; inspect before resuming")
                print(f"Skipping verified completed conversion: {run}", flush=True)
                continue
            if (directory / f"{run}.sra.lock").exists():
                raise ValueError(f"Existing SRA lock: {directory}. Check existing download first; do not force/delete it.")
            print(f"Processing {run}", flush=True)
            event = {"run": run, "started_at": datetime.now(timezone.utc).isoformat(),
                     "archive_directory": str(directory), "status": "started"}
            status = logs / f"{run}.status.json"
            status.write_text(json.dumps(event, indent=2) + "\n")
            try:
                invoke(["prefetch", run, "--resume", "yes", "--progress", "-O", directory.parent],
                       logs / f"{run}.prefetch.log")
                archive = directory / f"{run}.sra"
                if not archive.is_file():
                    raise ValueError(f"prefetch returned success but expected complete archive absent: {archive}")
                event["status"] = "downloaded"
                if args.convert:
                    sample = samples[run]
                    outputs = [Path(sample[k]) for k in ("fastq_1", "fastq_2") if sample[k]]
                    if any(p.exists() for p in outputs):
                        raise ValueError(f"FASTQ already exists for {run}; inspect it before conversion. No files overwritten.")
                    invoke(["vdb-validate", directory], logs / f"{run}.validate.log")
                    Path("data/tmp").mkdir(parents=True, exist_ok=True)
                    # Keep failed conversion files for diagnosis; never delete user's archives.
                    staging = Path(tempfile.mkdtemp(prefix=f"{run}-", dir="data/tmp"))
                    invoke(["fasterq-dump", directory.resolve(), "--split-3", "--threads", str(args.threads),
                            "--outdir", staging, "--temp", staging], logs / f"{run}.fasterq.log")
                    generated = []
                    names = [f"{run}_1", f"{run}_2"] if sample["layout"] == "PAIRED" else [run]
                    for name in names:
                        source = staging / f"{name}.fastq"
                        compressed = staging / f"{name}.fastq.gz"
                        with source.open("rb") as reader, gzip.open(compressed, "wb") as writer:
                            shutil.copyfileobj(reader, writer)
                        generated.append(compressed)
                    event["reads"] = validate_reads(generated)
                    for source, target in zip(generated, outputs):
                        target.parent.mkdir(parents=True, exist_ok=True)
                        source.replace(target)
                    for item, target in zip(event["reads"]["files"], outputs):
                        item["path"] = str(target)
                    # Only delete raw temporary files generated successfully by this invocation.
                    for name in names:
                        (staging / f"{name}.fastq").unlink()
                    event["staging_directory"] = str(staging)
                    event["status"] = "converted"
            except BaseException as error:
                event.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed", error=str(error))
                raise
            finally:
                event["updated_at"] = datetime.now(timezone.utc).isoformat()
                status.write_text(json.dumps(event, indent=2) + "\n")
                if event["status"] == "converted":
                    completed.write_text(json.dumps(event, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Interrupted: entire batch stopped; downloaded partial files preserved.", file=sys.stderr)
        sys.exit(130)
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        print(f"Download stopped: {error}", file=sys.stderr)
        sys.exit(1)
