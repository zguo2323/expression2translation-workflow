"""Validate local source evidence and the curated E2T sample sheet.

Metadata-only by design: never downloads data or opens planned FASTQ paths.
"""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform
import re
import sys

from jsonschema import Draft202012Validator
import yaml


class ValidationError(ValueError):
    """An actionable metadata/configuration failure."""


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_path(root, value):
    require(isinstance(value, str) and value.strip() == value and bool(value),
            f"Invalid relative path: {value!r}")
    path = Path(value)
    require(not path.is_absolute() and ".." not in path.parts and value != ".",
            f"Path must be repository-relative without '..': {value}")
    resolved = (root / path).resolve()
    require(resolved.is_relative_to(root), f"Path escapes repository: {value}")
    return resolved


def check_schema(data, path, label):
    schema = json.loads(path.read_text())
    errors = sorted(Draft202012Validator(schema).iter_errors(data), key=lambda e: str(e.path))
    if errors:
        details = "; ".join(f"{'.'.join(map(str, e.path)) or '<root>'}: {e.message}" for e in errors)
        raise ValidationError(f"{label}: {details}")


def read_table(path, delimiter):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing or duplicate header in {path.name}")
        rows = list(reader)
    require(rows, f"Empty table: {path.name}")
    for index, row in enumerate(rows, 2):
        require(None not in row and None not in row.values(),
                f"{path.name}:{index}: malformed row / missing column")
        require(all(v == v.strip() for v in row.values()),
                f"{path.name}:{index}: leading/trailing whitespace")
    return rows


def read_soft(path):
    samples = {}
    current = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                current = line.split(" = ", 1)[1]
                require(current not in samples, f"Duplicate GSM in SOFT: {current}")
                samples[current] = {}
            elif line.startswith("^"):
                current = None
            elif current and line.startswith("!Sample_") and " = " in line:
                key, value = line.split(" = ", 1)
                samples[current].setdefault(key, []).append(value)
    return samples


def unique(rows, field):
    values = [r[field] for r in rows]
    duplicates = sorted(v for v, count in Counter(values).items() if count > 1)
    require(not duplicates, f"Duplicate {field}: {duplicates}")


def validate_sample_rows(rows, root, design, raw_path, schema_path=None):
    """One row per biological library; no accession or file existence required."""
    schema_path = schema_path or root / "workflow/schemas/samples.schema.json"
    for row in rows:
        check_schema(row, schema_path, f"sample {row.get('sample_id')}")
    unique(rows, "sample_id")
    for key in ("run_accession", "geo_accession"):
        unique([r for r in rows if r.get(key)], key)
    identities = [(r["condition"], r["assay"], r["replicate"]) for r in rows]
    require(len(identities) == len(set(identities)),
            "Duplicate biological identity: condition × assay × replicate; merge technical runs/lanes upstream")
    # The same biological specimen cannot be counted twice within an assay.
    specimens = [(r["assay"], r["biosample"]) for r in rows if r.get("biosample")]
    require(len(specimens) == len(set(specimens)), "Duplicate BioSample within assay; not independent biological replicates")
    fastq_paths = set()
    raw_dir = relative_path(root, raw_path) if raw_path is not None else None
    for row in rows:
        sid = row["sample_id"]
        require(row["condition"] in design["conditions"], f"{sid}: condition absent from design.conditions")
        require(row["layout"] == design["assays"][row["assay"]], f"{sid}: assay/layout mismatch")
        require(bool(row["fastq_2"]) == (row["layout"] == "PAIRED"),
                f"{sid}: PAIRED requires fastq_2; SINGLE requires empty fastq_2")
        for value in (row["fastq_1"], row["fastq_2"]):
            if not value or raw_dir is None:
                continue
            path = relative_path(root, value)
            require(path.is_relative_to(raw_dir) and value.endswith((".fastq.gz", ".fq.gz")),
                    f"{sid}: FASTQ must be .fastq.gz/.fq.gz under paths.raw")
            require(path not in fastq_paths, f"{sid}: duplicate FASTQ path: {value}")
            fastq_paths.add(path)
    return fastq_paths


def select_samples(ids, samples, assay, layout):
    require(isinstance(ids, list) and all(isinstance(s, str) for s in ids)
            and len(ids) == len(set(ids)), f"Invalid {assay} sample_ids: expected unique list")
    if not ids:
        ids = [s["sample_id"] for s in samples if s["assay"] == assay]
    by_id = {s["sample_id"]: s for s in samples}
    require(ids and set(ids) <= set(by_id), f"{assay} selection contains unknown sample or no {assay} libraries")
    chosen = {sid: by_id[sid] for sid in ids}
    label = "paired-end RNA" if assay == "rnaseq" else "single-end riboseq"
    require(all(s["assay"] == assay and s["layout"] == layout for s in chosen.values()),
            f"Selection requires known {label} samples")
    identities = [(s["condition"], s["replicate"]) for s in chosen.values()]
    require(len(identities) == len(set(identities)), f"{assay}: duplicate biological replicate identity")
    return ids, chosen


def validate_geo_sra(config, root, rows):
    manifest_path = relative_path(root, config["source_manifest"])
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get("schema_version") == 1, "Unsupported source manifest version")
    entries = manifest.get("files")
    require(isinstance(entries, list) and bool(entries), "Empty source manifest")
    seen = set()
    for entry in entries:
        require(isinstance(entry, dict), "Invalid source manifest entry")
        path = relative_path(root, entry["path"])
        require(path not in seen, f"Duplicate source manifest path: {entry['path']}")
        seen.add(path)
        require(path.is_file(), f"Missing source file: {entry['path']}")
        require(isinstance(entry["sha256"], str) and re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]),
                f"Invalid SHA-256: {entry['path']}")
        require(path.stat().st_size == entry["bytes"] and sha256(path) == entry["sha256"],
                f"Source checksum/size mismatch: {entry['path']}")
    sources = {key: relative_path(root, value) for key, value in config["sources"].items()}
    require(all(path in seen for path in sources.values()), "Required source absent from manifest")
    accessions = sources["accessions"].read_text().splitlines()
    require(all(re.fullmatch(r"(?:SRR|ERR|DRR)[0-9]+", item) for item in accessions), "Invalid accession list")
    require(len(accessions) == len(set(accessions)), "Duplicate SRR in accession list")
    require({r["run_accession"] for r in rows} <= set(accessions),
            "Sample runs do not match selected accession list")
    sra_rows = read_table(sources["sra"], ",")
    unique(sra_rows, "Run")
    sra = {r["Run"]: r for r in sra_rows}
    geo = read_soft(sources["geo_soft"])
    design = config["design"]
    evidence = []
    for row in rows:
        sid, run, gsm = row["sample_id"], row["run_accession"], row["geo_accession"]
        require(run in sra, f"{sid}: SRR missing from SRA CSV")
        source = sra[run]
        checks = {"BioSample": row["biosample"], "Sample Name": gsm,
                  "Library Name": gsm, "AGE": row["condition"], "LibraryLayout": row["layout"],
                  "BioProject": config["dataset"]["bioproject"],
                  "SRA Study": config["dataset"]["sra_study"],
                  "Organism": config["dataset"]["organism"]}
        for key, expected_value in checks.items():
            require(source.get(key) == expected_value, f"{sid}: SRA {key} mismatch")
        require(row["layout"] == design["assays"][row["assay"]], f"{sid}: assay/layout mismatch")
        require(gsm in geo, f"{sid}: GSM missing from GEO SOFT")
        sample_geo = geo[gsm]
        label = {"rnaseq": "RNA", "riboseq": "Ribo"}[row["assay"]]
        title = f"{label}, {row['condition']}, {row['replicate']}"
        require(sample_geo.get("!Sample_title") == [title], f"{sid}: GEO title / assay / replicate mismatch")
        require(sample_geo.get("!Sample_series_id") == [config["dataset"]["geo"]],
                f"{sid}: GEO series mismatch")
        require(sample_geo.get("!Sample_organism_ch1") == [config["dataset"]["organism"]],
                f"{sid}: GEO organism mismatch")
        relations = sample_geo.get("!Sample_relation", [])
        require(any(v.startswith("BioSample:") and v.rstrip("/").endswith("/" + row["biosample"])
                    for v in relations), f"{sid}: GEO BioSample mismatch")
        require(any(v.startswith("SRA:") and v.endswith("term=" + source["Experiment"])
                    for v in relations), f"{sid}: GEO SRA experiment mismatch")
        if row["assay"] == "rnaseq":
            require(source["Assay Type"] == "RNA-Seq", f"{sid}: unexpected SRA RNA assay")
        archive_bytes, bases = int(source["Bytes"]), int(source["Bases"])
        require(archive_bytes > 0 and bases > 0, f"{sid}: invalid SRA size/bases")
        evidence.append({**row, "sra_assay_type": source["Assay Type"], "geo_title": title,
                         "archive_bytes": archive_bytes, "bases": bases})

    return evidence, manifest, [config["source_manifest"]] + [e["path"] for e in entries]


def _validate_metadata(config, root):
    check_schema(config, root / "workflow/schemas/config.schema.json", "config")
    for value in config["paths"].values():
        relative_path(root, value)
    for key in ("genome", "annotation", "transcriptome", "rrna"):
        if config["reference"][key] is not None:
            relative_path(root, config["reference"][key])
    if config["preprocessing"]["ribo_psite_offsets"] is not None:
        relative_path(root, config["preprocessing"]["ribo_psite_offsets"])
    sample_path = relative_path(root, config["samples"])
    rows = read_table(sample_path, "\t")
    design = config["design"]
    fastq_paths = validate_sample_rows(rows, root, design, config["paths"]["raw"])
    manifest, source_files = None, []
    if config["source_validation"] == "geo_sra":
        study_path = relative_path(root, config["study_samples"])
        study = read_table(study_path, "\t")
        # Historical study paths/conditions are independent of the active input location/subset.
        study_design = {**design, "conditions": sorted({r["condition"] for r in study})}
        validate_sample_rows(study, root, study_design, None)
        study_evidence, manifest, source_files = validate_geo_sra(config, root, study)
        by_id = {r["sample_id"]: r for r in study_evidence}
        evidence = []
        for row in rows:
            require(row["sample_id"] in by_id, f"Unknown study sample: {row['sample_id']}")
            original = by_id[row["sample_id"]]
            for key in ("run_accession", "geo_accession", "biosample", "condition", "assay", "replicate", "layout"):
                require(row.get(key) == original[key], f"{row['sample_id']}: study {key} mismatch")
            evidence.append({**original, **row})
        source_files.append(config["study_samples"])
    else:
        evidence = [{**row, "archive_bytes": None, "bases": None} for row in rows]
    total_bytes = sum(r["archive_bytes"] for r in evidence) if manifest else None
    total_bases = sum(r["bases"] for r in evidence) if manifest else None
    # Prevent configurable report/log locations from overwriting an input or future reads.
    protected = {sample_path} | {relative_path(root, p) for p in source_files}
    for key in ("results", "logs"):
        directory = relative_path(root, config["paths"][key])
        for protected_dir in ("config", "metadata", "src", "workflow", "envs", "docs", "tests", "resources", ".git"):
            other = root / protected_dir
            require(not directory.is_relative_to(other) and not other.is_relative_to(directory),
                    f"paths.{key} overlaps protected directory {protected_dir}")
        require(not any(p.is_relative_to(directory) for p in protected),
                f"paths.{key} overlaps metadata input")
        for other_key in ("raw", "sra", "reference"):
            other = relative_path(root, config["paths"][other_key])
            require(not directory.is_relative_to(other) and not other.is_relative_to(directory),
                    f"paths.{key} overlaps paths.{other_key}")

    unresolved = [f"{section}.{key}" for section in ("reference", "preprocessing", "integration")
                  for key, value in config[section].items() if value is None]
    return {
        "status": "passed", "stage": "metadata", "analysis_ready": False,
        "source_validation": config["source_validation"],
        "sample_count": len(rows), "expected_fastq_count": len(fastq_paths),
        "groups": [{"condition": c, "assay": a,
                    "count": sum(r["condition"] == c and r["assay"] == a for r in rows)}
                   for c in design["conditions"] for a in design["assays"]],
        "archive_bytes": total_bytes, "archive_gb": total_bytes / 1e9 if total_bytes is not None else None,
        "archive_gib": total_bytes / 1024 ** 3 if total_bytes is not None else None, "bases": total_bases,
        "source_files": source_files, "source_manifest": manifest,
        "samples": evidence, "unresolved_parameters": unresolved,
        "limitations": ["FASTQ and reference existence/content have not been checked.",
                        "Biological independence requires method evidence; identifiers alone cannot prove it.",
                        "Cross-assay biological pairing is unconfirmed.",
                        "Metadata success does not authorize or establish analysis readiness."],
    }


def validate_metadata(config, root=Path(".")):
    """Return evidence on success; raise ValidationError without writing files."""
    try:
        return _validate_metadata(config, Path(root).resolve())
    except ValidationError:
        raise
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        raise ValidationError(f"Cannot validate metadata: {error}") from error


def write_report(config, output):
    report = validate_metadata(config)
    expected_output = Path(config["paths"]["results"]) / "validation/metadata.json"
    require(Path(output) == expected_output, f"Output must be {expected_output}")
    root = Path.cwd()
    code_paths = sorted(p for directory in ("src/validation", "workflow", "envs")
                        for p in Path(directory).rglob("*")
                        if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    report["provenance"] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "effective_config": config,
        "effective_config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
        "samples_sha256": sha256(root / config["samples"]),
        "source_manifest_sha256": (sha256(root / config["source_manifest"])
                                   if report["source_manifest"] is not None else None),
        "source_files_sha256": {p: sha256(root / p) for p in report["source_files"]},
        "code_and_environment_sha256": {p.as_posix(): sha256(p) for p in code_paths},
        "software": {"python": platform.python_version()},
    }
    for package in ("snakemake", "PyYAML", "jsonschema"):
        try:
            report["provenance"]["software"][package] = version(package)
        except PackageNotFoundError:
            report["provenance"]["software"][package] = "not installed"
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(target)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--config", default="config/config.yaml")
    source.add_argument("--config-json", help="Effective Snakemake configuration including overrides")
    parser.add_argument("--output", help="Write results/validation/metadata.json; otherwise read-only")
    args = parser.parse_args()
    try:
        config = json.loads(args.config_json) if args.config_json else yaml.safe_load(Path(args.config).read_text())
        report = write_report(config, args.output) if args.output else validate_metadata(config)
    except (ValidationError, OSError, ValueError, yaml.YAMLError) as error:
        print(f"Validation failed: {error}", file=sys.stderr)
        return 1
    print(f"Metadata validation passed: {report['sample_count']} samples, "
          f"{report['expected_fastq_count']} planned FASTQ files. Analysis ready: false.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
