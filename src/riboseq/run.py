"""Command adapter for reproducible Ribo-seq workflow rules."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import re
import subprocess
import sys

from src.riboseq.config import load_offsets
from src.riboseq.features import derive_cds_features, read_features, write_features
from src.validation.sequences import fasta
from src.validation.validate import require, sha256


def _run(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, **kwargs)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    require(result.returncode == 0, f"Command failed ({result.returncode}): {' '.join(map(str, command))}")
    return result


def _version(command):
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    return (result.stdout or result.stderr).splitlines()[0].strip()


def _fastq_count(path):
    lines = 0
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for lines, _ in enumerate(handle, 1):
            pass
    require(lines and lines % 4 == 0, f"Invalid/empty FASTQ after decontamination: {path}")
    return lines // 4


def build_index(task):
    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / "index"
    _run(["bowtie2-build", "--threads", str(task["threads"]), task["reference"], str(prefix)])
    report = dict(role=task["role"], reference=task["reference"], reference_sha256=sha256(Path(task["reference"])),
                  bowtie2_version=_version(["bowtie2", "--version"]), synthetic=task["synthetic"])
    (directory / "index.json").write_text(json.dumps(report, indent=2) + "\n")


def make_features(task):
    lengths = {name: len(seq) for name, seq in fasta(task["transcriptome"]).items()}
    rows = derive_cds_features(task["annotation"], lengths)
    write_features(task["output"], rows)
    report = dict(transcripts=len(rows), genes=len({row["gene_id"] for row in rows}),
                  annotation_sha256=sha256(Path(task["annotation"])), transcriptome_sha256=sha256(Path(task["transcriptome"])))
    Path(task["report"]).write_text(json.dumps(report, indent=2) + "\n")


def decontaminate(task):
    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    clean = directory / "clean.fastq.gz"
    input_reads = _fastq_count(task["reads"])
    command = ["bowtie2", "--very-sensitive", "--threads", str(task["threads"]), "-x", task["index"],
               "-U", task["reads"], "--un-gz", str(clean), "-S", "/dev/null"]
    _run(command)
    retained = _fastq_count(clean)
    require(retained <= input_reads, "rRNA removal produced more reads than its input")
    fraction = (input_reads - retained) / input_reads
    require(fraction <= task["max_rrna_fraction"],
            f"rRNA fraction {fraction:.4f} exceeds configured maximum {task['max_rrna_fraction']}")
    report = dict(sample_id=task["sample_id"], input_reads=input_reads, rrna_reads=input_reads - retained,
                  retained_reads=retained, rrna_fraction=fraction, max_rrna_fraction=task["max_rrna_fraction"],
                  input_sha256=sha256(Path(task["reads"])), clean_sha256=sha256(clean),
                  bowtie2_version=_version(["bowtie2", "--version"]), synthetic=task["synthetic"])
    (directory / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")


def align(task):
    import pysam
    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    bam = directory / "alignment.bam"
    bowtie = ["bowtie2", "--very-sensitive", "--no-unal", "--threads", str(task["threads"]),
              "-x", task["index"], "-U", task["reads"]]
    sort = ["samtools", "sort", "-@", str(task["threads"]), "-o", str(bam), "-"]
    first = subprocess.Popen(bowtie, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    second = subprocess.Popen(sort, stdin=first.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    first.stdout.close()
    _, sort_error = second.communicate()
    _, bowtie_error = first.communicate()
    if bowtie_error:
        print(bowtie_error.decode(errors="replace"), end="", file=sys.stderr)
    if sort_error:
        print(sort_error.decode(errors="replace"), end="", file=sys.stderr)
    require(first.returncode == 0 and second.returncode == 0, "Bowtie2/SAMtools alignment pipeline failed")
    _run(["samtools", "index", str(bam)])
    flagstat = _run(["samtools", "flagstat", str(bam)]).stdout
    (directory / "flagstat.txt").write_text(flagstat)
    input_reads = _fastq_count(task["reads"])
    primary = 0
    with pysam.AlignmentFile(bam, "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                primary += 1
    rate = primary / input_reads
    require(rate >= task["min_mapping_rate"],
            f"Transcriptome mapping rate {rate:.4f} below configured minimum {task['min_mapping_rate']}")
    report = dict(sample_id=task["sample_id"], input_reads=input_reads, primary_mapped_reads=primary,
                  mapping_rate=rate, min_mapping_rate=task["min_mapping_rate"], bam_sha256=sha256(bam),
                  bowtie2_version=_version(["bowtie2", "--version"]),
                  samtools_version=_version(["samtools", "--version"]), synthetic=task["synthetic"])
    (directory / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")


def psite_coordinate(reference_start, reference_end, reverse, read_length, offset):
    require(0 <= offset < read_length, f"Offset {offset} outside read length {read_length}")
    return reference_end - 1 - offset if reverse else reference_start + offset


def _svg_plot(path, points, title, bars=False):
    width, height, left, right, top, bottom = 640, 360, 65, 20, 35, 50
    xs, ys = [point[0] for point in points], [point[1] for point in points]
    xmin, xmax, ymax = min(xs), max(xs), max(max(ys), 1)
    xspan = max(xmax - xmin, 1)
    sx = lambda value: left + (value - xmin) * (width - left - right) / xspan
    sy = lambda value: height - bottom - value * (height - top - bottom) / ymax
    shapes = []
    if bars:
        bar_width = (width - left - right) / max(len(points), 1) * 0.65
        shapes = [f'<rect x="{sx(x) - bar_width / 2:.2f}" y="{sy(y):.2f}" width="{bar_width:.2f}" '
                  f'height="{height - bottom - sy(y):.2f}" fill="#3b82f6"/>' for x, y in points]
    else:
        coordinates = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        shapes = [f'<polyline points="{coordinates}" fill="none" stroke="#2563eb" stroke-width="2"/>']
    labels = [f'<text x="{sx(x):.2f}" y="{height - bottom + 20}" text-anchor="middle" font-size="11">{x}</text>'
              for x in sorted(set((xmin, xmax) if len(points) > 8 else xs))]
    content = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
               '<rect width="100%" height="100%" fill="white"/>',
               f'<text x="{width / 2}" y="22" text-anchor="middle" font-size="16">{title}</text>',
               f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black"/>',
               f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black"/>',
               *shapes, *labels,
               f'<text x="{left - 8}" y="{sy(ymax) + 4:.2f}" text-anchor="end" font-size="11">{ymax}</text>',
               f'<text x="{left - 8}" y="{height-bottom + 4}" text-anchor="end" font-size="11">0</text>', '</svg>']
    Path(path).write_text("\n".join(content) + "\n")


def count_sites(task):
    import pysam
    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    features = read_features(task["features"])
    offsets = load_offsets(task["offset_table"])
    filters = Counter()
    lengths = defaultdict(Counter)
    counts = Counter()
    frames = Counter()
    metagene = Counter()
    missing_lengths = set()
    with pysam.AlignmentFile(task["bam"], "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped:
                filters["unmapped"] += 1
                continue
            if read.is_secondary or read.is_supplementary:
                filters["secondary_or_supplementary"] += 1
                continue
            length = read.query_length
            lengths[length]["primary_mapped"] += 1
            if read.mapping_quality < task["min_mapq"]:
                filters["low_mapq"] += 1
                continue
            if task["require_unique"] and read.has_tag("XS"):
                filters["non_unique"] += 1
                continue
            if ((task["strand_mode"] == "forward" and read.is_reverse) or
                    (task["strand_mode"] == "reverse" and not read.is_reverse)):
                filters["strand_mismatch"] += 1
                continue
            if not read.cigartuples or any(operation not in (0, 7, 8) for operation, _ in read.cigartuples):
                filters["gapped_or_clipped"] += 1
                continue
            if length not in offsets:
                missing_lengths.add(length)
                continue
            lengths[length]["offset_eligible"] += 1
            feature = features.get(read.reference_name)
            if feature is None:
                filters["noncoding_transcript"] += 1
                continue
            site = psite_coordinate(read.reference_start, read.reference_end, read.is_reverse, length, offsets[length])
            if not 0 <= site < handle.get_reference_length(read.reference_name):
                filters["psite_out_of_bounds"] += 1
                continue
            relative = site - feature["cds_start"]
            if -task["metagene_upstream"] <= relative <= task["metagene_downstream"]:
                metagene[relative] += 1
            if not feature["cds_start"] <= site < feature["cds_end"]:
                filters["outside_cds"] += 1
                continue
            counts[feature["gene_id"]] += 1
            frames[relative % 3] += 1
            lengths[length]["cds_psites"] += 1
    require(not missing_lengths, "P-site offset table does not cover eligible read lengths: " +
            ", ".join(map(str, sorted(missing_lengths))))
    assigned = sum(counts.values())
    require(assigned > 0, "No P-sites assigned inside CDS")
    frame0 = frames[0] / assigned
    require(frame0 >= task["min_frame0_fraction"],
            f"Frame-0 P-site fraction {frame0:.4f} below configured minimum {task['min_frame0_fraction']}")
    gene_lengths = {row["gene_id"]: row["cds_length"] for row in features.values()}
    with (directory / "gene_counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "p_site_count", "cds_length_nt", "p_sites_per_kb", "rpm"])
        for gene in sorted(gene_lengths):
            count = counts[gene]
            writer.writerow([gene, count, gene_lengths[gene], f"{count * 1000 / gene_lengths[gene]:.8f}",
                             f"{count * 1_000_000 / assigned:.8f}"])
    with (directory / "read_lengths.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["read_length", "psite_offset", "primary_mapped", "offset_eligible", "cds_psites"])
        for length in sorted(lengths):
            writer.writerow([length, offsets.get(length, ""), lengths[length]["primary_mapped"],
                             lengths[length]["offset_eligible"], lengths[length]["cds_psites"]])
    with (directory / "start_metagene.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["position_from_cds_start", "p_site_count"])
        for position in range(-task["metagene_upstream"], task["metagene_downstream"] + 1):
            writer.writerow([position, metagene[position]])
    with (directory / "periodicity.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["frame", "p_site_count", "fraction"])
        for frame in range(3):
            writer.writerow([frame, frames[frame], f"{frames[frame] / assigned:.8f}"])
    _svg_plot(directory / "start_metagene.svg",
              [(position, metagene[position]) for position in range(-task["metagene_upstream"], task["metagene_downstream"] + 1)],
              "P-site metagene around CDS start")
    _svg_plot(directory / "read_lengths.svg", [(length, lengths[length]["primary_mapped"]) for length in sorted(lengths)],
              "Primary mapped read lengths", bars=True)
    _svg_plot(directory / "periodicity.svg", [(frame, frames[frame]) for frame in range(3)],
              "P-site frame periodicity", bars=True)
    artifact_names = ("gene_counts.tsv", "read_lengths.tsv", "start_metagene.tsv", "periodicity.tsv",
                      "start_metagene.svg", "read_lengths.svg", "periodicity.svg")
    report = dict(sample_id=task["sample_id"], assigned_cds_psites=assigned,
                  frame_counts={str(i): frames[i] for i in range(3)}, frame0_fraction=frame0,
                  min_frame0_fraction=task["min_frame0_fraction"], filters=dict(sorted(filters.items())),
                  parameters={key: task[key] for key in ("strand_mode", "strand_evidence", "min_mapq", "require_unique",
                                                        "offset_evidence", "metagene_upstream", "metagene_downstream")},
                  inputs_sha256={key: sha256(Path(task[key])) for key in ("bam", "features", "offset_table")},
                  artifacts_sha256={name: sha256(directory / name) for name in artifact_names},
                  pysam_version=pysam.__version__, synthetic=task["synthetic"])
    (directory / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")


def summarize(task):
    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    qc_rows, count_rows, inputs = [], [], {}
    for sample in task["samples"]:
        decontam = json.loads(Path(sample["decontam_metrics"]).read_text())
        alignment = json.loads(Path(sample["alignment_metrics"]).read_text())
        sites = json.loads(Path(sample["site_metrics"]).read_text())
        require(all(item["sample_id"] == sample["sample_id"] for item in (decontam, alignment, sites)),
                f"Sample metrics disagree: {sample['sample_id']}")
        qc_rows.append([sample["sample_id"], sample["condition"], decontam["input_reads"], decontam["rrna_fraction"],
                        alignment["mapping_rate"], sites["assigned_cds_psites"], sites["frame0_fraction"]])
        with Path(sample["gene_counts"]).open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                count_rows.append([sample["sample_id"], sample["condition"], row["gene_id"], row["p_site_count"],
                                   row["cds_length_nt"], row["p_sites_per_kb"], row["rpm"]])
        for key in ("decontam_metrics", "alignment_metrics", "site_metrics", "gene_counts"):
            inputs[sample[key]] = sha256(Path(sample[key]))
    with (directory / "qc_metrics.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "condition", "input_reads", "rrna_fraction", "mapping_rate",
                         "assigned_cds_psites", "frame0_fraction"])
        writer.writerows(qc_rows)
    with (directory / "gene_counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "condition", "gene_id", "p_site_count", "cds_length_nt", "p_sites_per_kb", "rpm"])
        writer.writerows(count_rows)
    artifacts = {str(path): sha256(path) for path in (directory / "qc_metrics.tsv", directory / "gene_counts.tsv")}
    provenance = dict(reference_id=task["reference_id"], synthetic=task["synthetic"], effective_config=task["config"],
                      reference_manifest_sha256=sha256(Path(task["reference_manifest"])),
                      config_sha256=sha256(Path(task["config_file"])),
                      inputs_sha256=inputs, artifacts_sha256=artifacts)
    (directory / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


ACTIONS = {"index": build_index, "features": make_features, "decontaminate": decontaminate,
           "align": align, "count": count_sites, "summarize": summarize}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    task = json.loads(parser.parse_args().task)
    require(task.get("action") in ACTIONS, "Unknown Ribo action")
    ACTIONS[task["action"]](task)
