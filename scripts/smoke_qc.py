#!/usr/bin/env python3
"""Create isolated synthetic SE/PE data and exercise the real QC tools and DAG."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", default="results/synthetic-qc-smoke")
    parser.add_argument("--dry-run-only", action="store_true")
    parser.add_argument("--use-conda", action="store_true")
    parser.add_argument("--mode", choices=["raw", "trim"], default="trim")
    args = parser.parse_args()
    destination = Path(args.directory).resolve()
    if destination.exists():
        raise ValueError(f"Smoke directory exists; choose a new --directory: {destination}")
    destination.mkdir(parents=True)
    for directory in ("config", "metadata", "src", "workflow", "envs"):
        shutil.copytree(ROOT / directory, destination / directory, ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
    config = yaml.safe_load((destination / "config/config.yaml").read_text())
    config["stage"] = "qc"
    (destination / "config/config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    qc = yaml.safe_load((destination / "config/qc.yaml").read_text())
    qc.update(mode=args.mode, synthetic=True, sample_ids=["young_rna_1", "young_ribo_1"])
    adapter = "AGATCGGAAGAGCACACGTCT"
    for assay, policy in qc["policies"].items():
        policy.update(confirmed=True, evidence="Synthetic test: 4bp prefix, known adapter, Q40; not GEO experimental evidence",
                      adapter_mode="sequence", adapter_r1=adapter, adapter_r2=adapter if assay == "rnaseq" else None,
                      barcode_mode="fixed_front", trim_front1=4, trim_front2=4 if assay == "rnaseq" else 0, min_length=20)
    if args.mode == "raw":
        qc["policies"] = yaml.safe_load((ROOT / "config/qc.yaml").read_text())["policies"]
    (destination / "config/qc.yaml").write_text(yaml.safe_dump(qc, sort_keys=False))
    randomizer = random.Random(203147)
    with (destination / "config/samples.tsv").open() as handle:
        samples = [r for r in csv.DictReader(handle, delimiter="\t") if r["sample_id"] in qc["sample_ids"]]
    for sample in samples:
        for mate, key in enumerate(("fastq_1", "fastq_2"), 1):
            if not sample[key]:
                continue
            path = destination / sample[key]
            path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(path, "wt") as handle:
                for index in range(200):
                    length = 80 if sample["assay"] == "rnaseq" else 28
                    seq = "AAAA" + "".join(randomizer.choices("ACGT", k=length)) + adapter
                    handle.write(f"@synthetic_{index}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n")
    executable = shutil.which("snakemake") or str(ROOT / ".venv/bin/snakemake")
    command = [executable, "--snakefile", "workflow/Snakefile", "--cores", "2"]
    if args.use_conda:
        command += ["--use-conda"]
    for step, extra in [("dry-run", ["--dry-run"])] + ([] if args.dry_run_only else [("run", []), ("rerun", [])]):
        result = subprocess.run(command + extra, cwd=destination, text=True, capture_output=True)
        (destination / f"{step}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"Smoke {step} failed; see {destination / (step + '.log')}")
        if step == "rerun" and "Nothing to be done" not in result.stdout + result.stderr:
            raise RuntimeError("Unchanged smoke DAG unexpectedly reran")
    if not args.dry_run_only:
        for sample in samples:
            suffix = ".trimmed" if args.mode == "trim" else ""
            validation = json.loads((destination / f"results/qc/{args.mode}/validation/{sample['sample_id']}{suffix}.json").read_text())
            length = 80 if sample["assay"] == "rnaseq" else 28
            if args.mode == "raw":
                length += 4 + len(adapter)
            for item in validation["files"]:
                assert item["reads"] == 200 and item["bases"] == 200 * length, item
        report = destination / f"results/qc/{args.mode}/multiqc/multiqc_report.html"
        assert report.is_file() and report.stat().st_size > 1000
        summary = destination / f"results/qc/{args.mode}/multiqc/multiqc_report_data/multiqc_data.json"
        parsed = json.loads(summary.read_text())["report_saved_raw_data"]
        assert any("fastqc" in key for key in parsed), parsed.keys()
        assert any("fastp" in key for key in parsed) == (args.mode == "trim"), parsed.keys()
    print(f"Synthetic QC {'DAG' if args.dry_run_only else 'real-tool smoke and no-op rerun'} passed: {destination}")


if __name__ == "__main__":
    main()
