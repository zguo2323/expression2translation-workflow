#!/usr/bin/env python3
"""Isolated synthetic reference + 4 PE libraries through QC/Salmon/tximport/DESeq2."""
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
    rng = random.Random(203147)
    transcripts = {}
    genome = []
    gtf = []
    offset = 0
    # Introns exercise decoy-aware mapping with spliced reference sequences.
    for i in range(81):
        tid, gid = f"tx{i:03}", f"gene{i:03}"
        exon1 = "".join(rng.choices("ACGT", k=450))
        intron = "".join(rng.choices("ACGT", k=80))
        exon2 = "".join(rng.choices("ACGT", k=450))
        strand = "+" if i % 2 == 0 else "-"
        seq = exon1 + exon2
        transcripts[tid] = seq if strand == "+" else reverse_complement(seq)
        biotype = "rRNA" if i == 80 else "protein_coding"
        attr = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "{biotype}";'
        for start, end in ((offset + 1, offset + 450), (offset + 531, offset + 980)):
            gtf.append(f"chrSynthetic\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attr}\n")
        genome.append(exon1 + intron + exon2 + "N" * 50)
        offset += 1030
    source = target / "resources/reference/synthetic/source"
    source.mkdir(parents=True)
    write_fasta(source / "genome.fa", {"chrSynthetic": "".join(genome)})
    (source / "annotation.gtf").write_text("".join(gtf))
    subprocess.run([sys.executable, "scripts/prepare_reference.py", "--genome", "resources/reference/synthetic/source/genome.fa",
                    "--annotation", "resources/reference/synthetic/source/annotation.gtf", "--output", "resources/reference/synthetic/derived",
                    "--reference-id", "synthetic-rna-v1", "--provider", "synthetic", "--release", "test-v1",
                    "--genome-source", "generated synthetic", "--annotation-source", "generated synthetic"], cwd=target, check=True)
    config = yaml.safe_load((target / "config/config.yaml").read_text())
    config["stage"] = "rnaseq"
    manifest_path = "resources/reference/synthetic/derived/manifest.json"
    manifest = json.loads((target / manifest_path).read_text())
    config["reference"] = {k: manifest[k] for k in ("reference_id", "provider", "release")}
    config["reference"].update({role: manifest["files"][role]["path"] for role in ("genome", "annotation", "transcriptome", "rrna")})
    (target / "config/config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    rna = yaml.safe_load((target / "config/rnaseq.yaml").read_text())
    rna.update(reference_manifest=manifest_path, synthetic=True, library_type="ISF", library_evidence="synthetic first mate forward",
               design_confirmed=True, design_evidence="synthetic independent libraries with planted expression changes", fit_type="mean")
    (target / "config/rnaseq.yaml").write_text(yaml.safe_dump(rna, sort_keys=False))
    qc = yaml.safe_load((target / "config/qc.yaml").read_text())
    qc.update(mode="trim", synthetic=True)
    qc["policies"]["rnaseq"].update(confirmed=True, evidence="Synthetic Q40 fragments without adapters/barcodes", adapter_mode="none",
        barcode_mode="none", min_length=50)
    (target / "config/qc.yaml").write_text(yaml.safe_dump(qc, sort_keys=False))
    with (target / "config/samples.tsv").open() as handle:
        samples = [s for s in csv.DictReader(handle, delimiter="\t") if s["assay"] == "rnaseq"]
    for sample in samples:
        p1, p2 = (target / sample[k] for k in ("fastq_1", "fastq_2"))
        p1.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(p1, "wt") as r1, gzip.open(p2, "wt") as r2:
            n = 0
            for i, (tid, seq) in enumerate(transcripts.items()):
                if i == 80:
                    continue  # absent rRNA must survive import as zero-count then be filtered
                count = 80 + (i * 17) % 120
                factor = (4 if i < 8 else 0.25 if i < 16 else 1) if sample["condition"] == "Middle" else 1
                count = max(10, int(count * factor * rng.uniform(0.85, 1.15)))
                for _ in range(count):
                    fragment_length = rng.randint(220, 300)
                    start = rng.randint(0, len(seq) - fragment_length)
                    fragment = seq[start:start + fragment_length]
                    for mate, read, handle in ((1, fragment[:75], r1), (2, reverse_complement(fragment[-75:]), r2)):
                        handle.write(f"@synthetic_{n}/{mate}\n{read}\n+\n{'I'*75}\n")
                    n += 1
    env = dict(os.environ)
    if not args.use_conda:
        env["PATH"] = os.pathsep.join([str(ROOT / ".conda/qc/bin"), str(ROOT / ".conda/salmon/bin"),
                                      str(ROOT / ".conda/rnaseq-stats/bin"), env["PATH"]])
    cmd = [shutil.which("snakemake") or str(ROOT / ".venv/bin/snakemake"), "--cores", "2", "--snakefile", "workflow/Snakefile"]
    if args.use_conda:
        cmd += ["--use-conda"]
    for step, extra in (("dry-run", ["--dry-run"]), ("run", []), ("rerun", [])):
        result = subprocess.run(cmd + extra, cwd=target, env=env, text=True, capture_output=True)
        (target / f"{step}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"RNA smoke failed: {target / (step + '.log')}")
        if step == "rerun":
            assert "Nothing to be done" in result.stdout + result.stderr
    out = target / "results/rnaseq/gene_analysis"
    report = json.loads((out / "analysis.json").read_text())
    hashes = json.loads((out / "provenance.json").read_text())
    assert len(hashes["inputs_sha256"]) == 7 and any("differential_expression.tsv" in p for p in hashes["artifacts_sha256"])
    assert report["synthetic"] and report["deseq2_run"] and report["genes_imported"] == 81
    with (out / "differential_expression.tsv").open() as handle:
        results = {r["gene_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert all(float(results[f"gene{i:03}"]["log2FoldChange"]) > 1 for i in range(8))
    assert all(float(results[f"gene{i:03}"]["log2FoldChange"]) < -1 for i in range(8, 16))
    assert "gene080" not in results
    # A real parameter change must invalidate quantification and downstream statistics.
    rna["library_type"] = "A"
    (target / "config/rnaseq.yaml").write_text(yaml.safe_dump(rna, sort_keys=False))
    changed = subprocess.run(cmd + ["--dry-run"], cwd=target, env=env, text=True, capture_output=True)
    (target / "parameter-change.log").write_text(changed.stdout + changed.stderr)
    assert changed.returncode == 0 and "salmon_quant" in changed.stdout + changed.stderr
    rna["library_type"] = "ISF"
    (target / "config/rnaseq.yaml").write_text(yaml.safe_dump(rna, sort_keys=False))
    restored = subprocess.run(cmd, cwd=target, env=env, text=True, capture_output=True)
    (target / "restored-config.log").write_text(restored.stdout + restored.stderr)
    assert restored.returncode == 0
    print(f"RNA real-tool smoke passed; 81 genes, planted up/down directions, no-op rerun and parameter invalidation: {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", default="results/synthetic-rna-smoke")
    parser.add_argument("--use-conda", action="store_true")
    smoke(parser.parse_args())
