#!/usr/bin/env python3
"""Derive consistent transcript/rRNA FASTA and tx2gene from a single genome/GTF."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.validation.sequences import reconstruct, write_fasta
from src.validation.validate import require, sha256, relative_path


def prepare(genome_path, annotation_path, output, reference_id, provider, release, sources):
    root = Path.cwd()
    genome_path = relative_path(root, str(genome_path))
    annotation_path = relative_path(root, str(annotation_path))
    output = relative_path(root, str(output))
    require(not output.exists(), f"Reference output exists; choose a new directory: {output}")
    genome, transcripts = reconstruct(genome_path, annotation_path)
    rrna = {tid: t["sequence"] for tid, t in transcripts.items() if t["biotype"].lower() == "rrna"}
    require(rrna, "No rRNA exons found: verify GTF biotypes before preparing reference")
    output.mkdir(parents=True)
    write_fasta(output / "transcripts.fa", {tid: t["sequence"] for tid, t in transcripts.items()})
    write_fasta(output / "rrna.fa", rrna)
    with (output / "tx2gene.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "gene_id"])
        writer.writerows((tid, t["gene_id"]) for tid, t in transcripts.items())
    timestamp = datetime.now(timezone.utc).isoformat()
    entries = {}
    for role, path in dict(genome=genome_path, annotation=annotation_path, transcriptome=output / "transcripts.fa", rrna=output / "rrna.fa", tx2gene=output / "tx2gene.tsv").items():
        entries[role] = dict(path=str(path.relative_to(root)), sha256=sha256(path), source=sources.get(role, "derived from manifest genome + annotation"),
                             retrieved_at=timestamp, derivation=None if role in sources else "scripts/prepare_reference.py: sorted exons, 1-based inclusive GTF, reverse-complement negative strand")
    manifest = dict(reference_id=reference_id, provider=provider, release=release, files=entries)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "preparation.json").write_text(json.dumps(dict(genome_contigs=len(genome), transcripts=len(transcripts), genes=len({t['gene_id'] for t in transcripts.values()}),
        rrna_transcripts=len(rrna), prepared_at=timestamp, sources_sha256={k: entries[k]["sha256"] for k in sources},
        code_sha256={p: sha256(Path(p)) for p in ("scripts/prepare_reference.py", "src/validation/sequences.py")}), indent=2) + "\n")
    return output / "manifest.json"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("genome", "annotation", "output", "reference-id", "provider", "release", "genome-source", "annotation-source"):
        p.add_argument("--" + name, required=True)
    args = p.parse_args()
    print(prepare(args.genome, args.annotation, args.output, args.reference_id, args.provider, args.release,
                  {"genome": args.genome_source, "annotation": args.annotation_source}))
