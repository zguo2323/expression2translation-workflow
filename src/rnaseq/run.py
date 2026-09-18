"""Salmon adapters with output/mapping gates and explicit reference identity."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from src.validation.reference import validate_rna_reference
from src.validation.sequences import fasta, write_fasta
from src.validation.validate import require, sha256


def run(task):
    out = Path(task["directory"])
    out.mkdir(parents=True, exist_ok=True)
    command = ["salmon"]
    version = subprocess.check_output(["salmon", "--version"], text=True).strip()
    if task["action"] == "index":
        ref = validate_rna_reference(task["manifest"])
        sequences = fasta(ref["files"]["transcriptome"]["path"])
        genome = fasta(ref["files"]["genome"]["path"])
        require(not set(sequences) & set(genome), "Decoy/transcript names overlap")
        write_fasta(out / "gentrome.fa", {**sequences, **genome})
        (out / "decoys.txt").write_text("\n".join(genome) + "\n")
        command += ["index", "-t", str(out / "gentrome.fa"), "-d", str(out / "decoys.txt"),
                    "-i", str(out / "index"), "-k", str(task["k"]), "-p", str(task["threads"]), "--keepDuplicates"]
        subprocess.run(command, check=True)
        require((out / "index/versionInfo.json").is_file(), "Salmon index missing versionInfo.json")
        metrics = {"reference_id": ref["reference_id"], "transcripts": len(sequences), "decoys": len(genome)}
    else:
        command += ["quant", "-i", task["index"], "-l", task["library_type"],
                    "-1", task["reads"][0], "-2", task["reads"][1], "-p", str(task["threads"]),
                    "--validateMappings", "-o", str(out)]
        subprocess.run(command, check=True)
        metrics = json.loads((out / "aux_info/meta_info.json").read_text())
        formats = json.loads((out / "lib_format_counts.json").read_text())
        require(metrics["num_processed"] > 0 and metrics["num_mapped"] > 0, "Salmon assigned no fragments")
        require(metrics["percent_mapped"] / 100 >= task["min_mapping_rate"], "Salmon mapping rate below configured gate")
        require(formats["compatible_fragment_ratio"] >= task["min_library_compatibility"], "Salmon library compatibility below gate; review library_type")
        metrics["library_format_evidence"] = formats
    provenance = dict(tool_version=version, command=command, parameters=task, metrics=metrics,
                      reference_manifest_sha256=sha256(Path(task["manifest"])),
                      code_sha256=sha256(Path(__file__)), synthetic=task["synthetic"],
                      generated_at=datetime.now(timezone.utc).isoformat())
    (out / "e2t.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    run(json.loads(parser.parse_args().task))
