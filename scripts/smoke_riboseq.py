#!/usr/bin/env python3
"""Run four synthetic Ribo libraries through real QC, Bowtie2, SAMtools and P-site counting."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.validation.sequences import reverse_complement, write_fasta


def smoke(args):
    target = Path(args.directory).resolve()
    if target.exists():
        raise ValueError(f"Choose a new smoke directory: {target}")
    target.mkdir(parents=True)
    for folder in ("config", "metadata", "src", "workflow", "envs", "scripts"):
        shutil.copytree(ROOT / folder, target / folder, ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))

    rng = random.Random(375616)
    genome_parts, gtf, transcripts = [], [], {}
    cursor = 0
    for index in range(50):
        tid, gid = f"tx{index:03}", f"gene{index:03}"
        transcript = "".join(rng.choices("ACGT", k=720))
        strand = "+" if index % 2 == 0 else "-"
        genomic = transcript if strand == "+" else reverse_complement(transcript)
        start, end = cursor + 1, cursor + len(genomic)
        attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "protein_coding";'
        gtf.append(f"chrSynthetic\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
        # The symmetric 60 nt UTRs produce transcript-coordinate CDS [60, 660) on both strands.
        gtf.append(f"chrSynthetic\ttest\tCDS\t{start + 60}\t{start + 659}\t.\t{strand}\t0\t{attrs}\n")
        transcripts[tid] = transcript
        genome_parts.extend([genomic, "N" * 80])
        cursor += 800
    rrna = "".join(rng.choices("ACGT", k=900))
    rrna_start = cursor + 1
    rrna_attrs = 'gene_id "rrna_gene"; transcript_id "rrna_tx"; gene_biotype "rRNA";'
    gtf.append(f"chrSynthetic\ttest\texon\t{rrna_start}\t{rrna_start + 899}\t.\t+\t.\t{rrna_attrs}\n")
    transcripts["rrna_tx"] = rrna
    genome_parts.append(rrna)

    source = target / "resources/reference/synthetic/source"
    source.mkdir(parents=True)
    write_fasta(source / "genome.fa", {"chrSynthetic": "".join(genome_parts)})
    (source / "annotation.gtf").write_text("".join(gtf))
    subprocess.run([sys.executable, "scripts/prepare_reference.py", "--genome", "resources/reference/synthetic/source/genome.fa",
                    "--annotation", "resources/reference/synthetic/source/annotation.gtf", "--output", "resources/reference/synthetic/derived",
                    "--reference-id", "synthetic-ribo-v1", "--provider", "synthetic", "--release", "test-v1",
                    "--genome-source", "generated synthetic", "--annotation-source", "generated synthetic"], cwd=target, check=True)
    manifest_path = "resources/reference/synthetic/derived/manifest.json"
    manifest = json.loads((target / manifest_path).read_text())

    config = yaml.safe_load((target / "config/config.yaml").read_text())
    config["stage"] = "riboseq"
    config["reference"] = {key: manifest[key] for key in ("reference_id", "provider", "release")}
    config["reference"].update({role: manifest["files"][role]["path"] for role in ("genome", "annotation", "transcriptome", "rrna")})
    (target / "config/config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

    offsets = target / "config/synthetic_psite_offsets.tsv"
    offsets.write_text("read_length\tpsite_offset\n28\t12\n")
    ribo = yaml.safe_load((target / "config/riboseq.yaml").read_text())
    ribo.update(reference_manifest=manifest_path, synthetic=True, strand_mode="forward",
                strand_evidence="synthetic reads generated in transcript orientation",
                offset_table="config/synthetic_psite_offsets.tsv", offset_evidence="synthetic P-sites planted at frame 0",
                min_mapping_rate=0.9, max_rrna_fraction=0.2, min_frame0_fraction=0.95)
    (target / "config/riboseq.yaml").write_text(yaml.safe_dump(ribo, sort_keys=False))
    qc = yaml.safe_load((target / "config/qc.yaml").read_text())
    qc.update(mode="trim", synthetic=True)
    qc["policies"]["riboseq"].update(confirmed=True, evidence="synthetic fixed 28 nt reads without adapters/barcodes",
        adapter_mode="none", barcode_mode="none", min_length=20)
    (target / "config/qc.yaml").write_text(yaml.safe_dump(qc, sort_keys=False))

    with (target / "config/samples.tsv").open() as handle:
        samples = [row for row in csv.DictReader(handle, delimiter="\t") if row["assay"] == "riboseq"]
    for sample_index, sample in enumerate(samples):
        path = target / sample["fastq_1"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as handle:
            read_number = 0
            for gene_index in range(50):
                sequence = transcripts[f"tx{gene_index:03}"]
                observations = 6 + (gene_index + sample_index) % 5
                for observation in range(observations):
                    site = 60 + 30 + 3 * ((gene_index * 7 + observation * 11 + sample_index) % 160)
                    start = site - 12
                    read = sequence[start:start + 28]
                    handle.write(f"@{sample['sample_id']}_{read_number}\n{read}\n+\n{'I' * 28}\n")
                    read_number += 1
            for contamination in range(20):
                start = (contamination * 31 + sample_index * 7) % (len(rrna) - 28)
                read = rrna[start:start + 28]
                handle.write(f"@{sample['sample_id']}_rrna_{contamination}\n{read}\n+\n{'I' * 28}\n")

    environment = dict(os.environ)
    if not args.use_conda:
        environment["PATH"] = os.pathsep.join([str(ROOT / ".conda/riboseq/bin"), str(ROOT / ".conda/qc/bin"), environment["PATH"]])
    command = [shutil.which("snakemake") or str(ROOT / ".venv/bin/snakemake"),
               "--cores", "2", "--snakefile", "workflow/Snakefile"]
    if args.use_conda:
        command += ["--use-conda"]
    for step, extra in (("dry-run", ["--dry-run"]), ("run", []), ("rerun", [])):
        result = subprocess.run(command + extra, cwd=target, env=environment, text=True, capture_output=True)
        (target / f"{step}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"Ribo smoke failed: {target / (step + '.log')}")
        if step == "rerun":
            assert "Nothing to be done" in result.stdout + result.stderr

    summary = target / "results/riboseq/summary"
    with (summary / "qc_metrics.tsv").open() as handle:
        qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (summary / "gene_counts.tsv").open() as handle:
        count_rows = list(csv.DictReader(handle, delimiter="\t"))
    provenance = json.loads((summary / "provenance.json").read_text())
    assert len(qc_rows) == 4 and len(count_rows) == 200
    assert all(float(row["frame0_fraction"]) >= 0.95 for row in qc_rows)
    assert all(float(row["mapping_rate"]) >= 0.9 for row in qc_rows)
    assert provenance["synthetic"] and provenance["reference_id"] == "synthetic-ribo-v1"
    assert all((target / f"results/riboseq/sites/{sample['sample_id']}/periodicity.svg").is_file() for sample in samples)

    ribo["min_mapq"] = 21
    (target / "config/riboseq.yaml").write_text(yaml.safe_dump(ribo, sort_keys=False))
    changed = subprocess.run(command + ["--dry-run"], cwd=target, env=environment, text=True, capture_output=True)
    (target / "parameter-change.log").write_text(changed.stdout + changed.stderr)
    assert changed.returncode == 0 and "ribo_count_sites" in changed.stdout + changed.stderr
    ribo["min_mapq"] = 20
    (target / "config/riboseq.yaml").write_text(yaml.safe_dump(ribo, sort_keys=False))
    restored = subprocess.run(command, cwd=target, env=environment, text=True, capture_output=True)
    (target / "restored-config.log").write_text(restored.stdout + restored.stderr)
    assert restored.returncode == 0
    print(f"Ribo real-tool smoke passed; 4 libraries, 50 genes, frame-0 P-sites and parameter invalidation: {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", default="results/synthetic-ribo-smoke")
    parser.add_argument("--use-conda", action="store_true")
    smoke(parser.parse_args())
