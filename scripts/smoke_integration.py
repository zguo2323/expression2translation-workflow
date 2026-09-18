#!/usr/bin/env python3
"""Run a synthetic RNA+Ribo dataset through the complete E2T MVP."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import sqlite3
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
    rng = random.Random(203147375616)
    genome_parts, gtf, transcripts = [], [], {}
    cursor = 0
    for index in range(50):
        tid, gid = f"tx{index:03}", f"gene{index:03}"
        transcript = "".join(rng.choices("ACGT", k=720))
        strand = "+" if index % 2 == 0 else "-"
        genomic = transcript if strand == "+" else reverse_complement(transcript)
        start, end = cursor + 1, cursor + 720
        attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "protein_coding";'
        gtf.extend([f"chrSynthetic\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n",
                    f"chrSynthetic\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n",
                    f"chrSynthetic\ttest\tCDS\t{start+60}\t{start+659}\t.\t{strand}\t0\t{attrs}\n"])
        transcripts[tid] = transcript
        genome_parts.extend([genomic, "N" * 80])
        cursor += 800
    rrna = "".join(rng.choices("ACGT", k=900))
    start = cursor + 1
    attrs = 'gene_id "rrna_gene"; transcript_id "rrna_tx"; gene_biotype "rRNA";'
    gtf.extend([f"chrSynthetic\ttest\tgene\t{start}\t{start+899}\t.\t+\t.\t{attrs}\n",
                f"chrSynthetic\ttest\texon\t{start}\t{start+899}\t.\t+\t.\t{attrs}\n"])
    transcripts["rrna_tx"] = rrna
    genome_parts.append(rrna)
    source = target / "resources/reference/synthetic/source"
    source.mkdir(parents=True)
    write_fasta(source / "genome.fa", {"chrSynthetic": "".join(genome_parts)})
    (source / "annotation.gtf").write_text("".join(gtf))
    subprocess.run([sys.executable, "scripts/prepare_reference.py", "--genome", "resources/reference/synthetic/source/genome.fa",
                    "--annotation", "resources/reference/synthetic/source/annotation.gtf", "--output", "resources/reference/synthetic/derived",
                    "--reference-id", "synthetic-e2e-v1", "--provider", "synthetic", "--release", "test-v1",
                    "--genome-source", "generated synthetic", "--annotation-source", "generated synthetic"], cwd=target, check=True)
    manifest_path = "resources/reference/synthetic/derived/manifest.json"
    manifest = json.loads((target / manifest_path).read_text())
    config = yaml.safe_load((target / "config/config.yaml").read_text())
    config["stage"] = "integration"
    config["reference"] = {key: manifest[key] for key in ("reference_id", "provider", "release")}
    config["reference"].update({role: manifest["files"][role]["path"] for role in ("genome", "annotation", "transcriptome", "rrna")})
    (target / "config/config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    rna = yaml.safe_load((target / "config/rnaseq.yaml").read_text())
    rna.update(reference_manifest=manifest_path, synthetic=True, library_type="ISF",
               library_evidence="synthetic inward-facing transcript-oriented fragments", run_deseq2=False,
               design_confirmed=False, design_evidence=None, min_samples=2)
    (target / "config/rnaseq.yaml").write_text(yaml.safe_dump(rna, sort_keys=False))
    (target / "config/synthetic_psite_offsets.tsv").write_text("read_length\tpsite_offset\n28\t12\n")
    ribo = yaml.safe_load((target / "config/riboseq.yaml").read_text())
    ribo.update(reference_manifest=manifest_path, synthetic=True, strand_mode="forward",
                strand_evidence="synthetic transcript-oriented reads", offset_table="config/synthetic_psite_offsets.tsv",
                offset_evidence="synthetic frame-0 P-sites", min_mapping_rate=0.9,
                max_rrna_fraction=0.2, min_frame0_fraction=0.95)
    (target / "config/riboseq.yaml").write_text(yaml.safe_dump(ribo, sort_keys=False))
    integration = yaml.safe_load((target / "config/integration.yaml").read_text())
    integration.update(synthetic=True, min_rna_tpm=1.0, min_ribo_tpm=1.0, report_title="Synthetic E2T end-to-end report")
    (target / "config/integration.yaml").write_text(yaml.safe_dump(integration, sort_keys=False))
    qc = yaml.safe_load((target / "config/qc.yaml").read_text())
    qc.update(mode="trim", synthetic=True)
    qc["policies"]["rnaseq"].update(confirmed=True, evidence="synthetic reads without adapters/barcodes",
        adapter_mode="none", barcode_mode="none", min_length=50)
    qc["policies"]["riboseq"].update(confirmed=True, evidence="synthetic 28 nt reads without adapters/barcodes",
        adapter_mode="none", barcode_mode="none", min_length=20)
    (target / "config/qc.yaml").write_text(yaml.safe_dump(qc, sort_keys=False))
    with (target / "config/samples.tsv").open() as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    for sample_index, sample in enumerate(samples):
        first = target / sample["fastq_1"]
        first.parent.mkdir(parents=True, exist_ok=True)
        if sample["assay"] == "rnaseq":
            second = target / sample["fastq_2"]
            with gzip.open(first, "wt") as r1, gzip.open(second, "wt") as r2:
                number = 0
                for gene_index in range(50):
                    sequence = transcripts[f"tx{gene_index:03}"]
                    for observation in range(12):
                        fragment_start = 80 + ((gene_index * 19 + observation * 13 + sample_index * 7) % 380)
                        fragment = sequence[fragment_start:fragment_start + 200]
                        r1.write(f"@{sample['sample_id']}_{number}/1\n{fragment[:75]}\n+\n{'I'*75}\n")
                        r2.write(f"@{sample['sample_id']}_{number}/2\n{reverse_complement(fragment[-75:])}\n+\n{'I'*75}\n")
                        number += 1
        else:
            with gzip.open(first, "wt") as handle:
                number = 0
                for gene_index in range(50):
                    sequence = transcripts[f"tx{gene_index:03}"]
                    observations = 24 if sample["condition"] == "Middle" and gene_index < 5 else 8
                    for observation in range(observations):
                        site = 90 + 3 * ((gene_index * 7 + observation * 11 + sample_index) % 160)
                        read = sequence[site - 12:site - 12 + 28]
                        handle.write(f"@{sample['sample_id']}_{number}\n{read}\n+\n{'I'*28}\n")
                        number += 1
                for contamination in range(20):
                    offset = (contamination * 31 + sample_index * 7) % (len(rrna) - 28)
                    handle.write(f"@{sample['sample_id']}_rrna_{contamination}\n{rrna[offset:offset+28]}\n+\n{'I'*28}\n")
    environment = dict(os.environ)
    if not args.use_conda:
        environment["PATH"] = os.pathsep.join([str(ROOT / ".conda/riboseq/bin"), str(ROOT / ".conda/qc/bin"),
            str(ROOT / ".conda/salmon/bin"), str(ROOT / ".conda/rnaseq-stats/bin"), environment["PATH"]])
    command = [str(ROOT / ".venv/bin/snakemake"), "--cores", "2", "--snakefile", "workflow/Snakefile"]
    if args.use_conda:
        command += ["--use-conda"]
    for step, extra in (("dry-run", ["--dry-run"]), ("run", []), ("rerun", [])):
        result = subprocess.run(command + extra, cwd=target, env=environment, text=True, capture_output=True)
        (target / f"{step}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"E2E smoke failed: {target / (step + '.log')}")
        if step == "rerun":
            assert "Nothing to be done" in result.stdout + result.stderr
    database = target / "results/catalog.db"
    connection = sqlite3.connect(database)
    try:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute("SELECT count(*) FROM run_samples").fetchone()[0] == 8
        assert connection.execute("SELECT count(*) FROM integration_metrics").fetchone()[0] == 100
        positive = connection.execute("SELECT gene_id,te_log2_change FROM integration_contrasts WHERE eligible=1 ORDER BY te_log2_change DESC LIMIT 5").fetchall()
        assert {gene for gene, _ in positive} == {f"gene{i:03}" for i in range(5)} and all(value > 1 for _, value in positive)
    finally:
        connection.close()
    assert (target / "results/report/report.html").is_file()
    assert "SYNTHETIC TEST" in (target / "results/report/report.html").read_text()
    integration["pseudocount"] = 2.0
    (target / "config/integration.yaml").write_text(yaml.safe_dump(integration, sort_keys=False))
    changed = subprocess.run(command + ["--dry-run"], cwd=target, env=environment, text=True, capture_output=True)
    (target / "parameter-change.log").write_text(changed.stdout + changed.stderr)
    text = changed.stdout + changed.stderr
    assert changed.returncode == 0 and "integrate_expression_translation" in text and "build_catalog" in text and "render_report" in text
    assert "salmon_quant" not in text and "ribo_count_sites" not in text
    integration["pseudocount"] = 1.0
    (target / "config/integration.yaml").write_text(yaml.safe_dump(integration, sort_keys=False))
    restored = subprocess.run(command, cwd=target, env=environment, text=True, capture_output=True)
    (target / "restored-config.log").write_text(restored.stdout + restored.stderr)
    assert restored.returncode == 0
    print(f"E2T complete synthetic smoke passed; RNA+Ribo, descriptive TE, SQLite and reports: {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", default="results/synthetic-e2e-smoke")
    parser.add_argument("--use-conda", action="store_true")
    smoke(parser.parse_args())
