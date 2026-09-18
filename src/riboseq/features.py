"""Derive coding intervals in transcript coordinates from the reference GTF."""
import csv
from collections import defaultdict
from pathlib import Path
import re

from src.validation.sequences import open_text
from src.validation.validate import require


FIELDS = ["transcript_id", "gene_id", "cds_start", "cds_end", "cds_length", "strand"]


def _attributes(text):
    return dict(re.findall(r'(\w+) "([^"]*)"', text))


def derive_cds_features(annotation_path, transcript_lengths):
    exons, coding, identities = defaultdict(list), defaultdict(list), {}
    with open_text(annotation_path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            require(len(fields) == 9, f"GTF line {number}: expected 9 columns")
            contig, _, feature, start, end, _, strand, _, attributes = fields
            if feature not in ("exon", "CDS"):
                continue
            attrs = _attributes(attributes)
            require("transcript_id" in attrs and "gene_id" in attrs, f"{feature} missing gene/transcript ID: {number}")
            require(strand in ("+", "-"), f"Invalid strand at GTF line {number}")
            tid, gid = attrs["transcript_id"], attrs["gene_id"]
            require(tid in transcript_lengths, f"GTF transcript absent from transcriptome: {tid}")
            identity = (gid, contig, strand)
            require(tid not in identities or identities[tid] == identity, f"Conflicting coding identity: {tid}")
            identities[tid] = identity
            interval = (int(start) - 1, int(end))
            require(0 <= interval[0] < interval[1], f"Invalid GTF coordinates at line {number}")
            (exons if feature == "exon" else coding)[tid].append(interval)

    rows = []
    genes = set()
    for tid in sorted(coding):
        require(exons[tid], f"CDS transcript has no exon: {tid}")
        ordered_exons = sorted(exons[tid])
        require(sum(end - start for start, end in ordered_exons) == transcript_lengths[tid],
                f"Transcript length disagrees with GTF exons: {tid}")
        exon_offsets, cursor = [], 0
        for start, end in ordered_exons:
            exon_offsets.append((start, end, cursor))
            cursor += end - start
        mapped = []
        for cds_start, cds_end in sorted(coding[tid]):
            hits = []
            for exon_start, exon_end, offset in exon_offsets:
                left, right = max(cds_start, exon_start), min(cds_end, exon_end)
                if left < right:
                    hits.append((offset + left - exon_start, offset + right - exon_start))
            require(sum(b - a for a, b in hits) == cds_end - cds_start,
                    f"CDS is not fully contained in exons: {tid}")
            mapped.extend(hits)
        gid, _, strand = identities[tid]
        if strand == "-":
            length = transcript_lengths[tid]
            mapped = [(length - end, length - start) for start, end in mapped]
        mapped.sort()
        require(all(a[1] == b[0] for a, b in zip(mapped, mapped[1:])),
                f"CDS is not contiguous in spliced transcript: {tid}")
        cds_start, cds_end = mapped[0][0], mapped[-1][1]
        require(0 <= cds_start < cds_end <= transcript_lengths[tid], f"CDS outside transcript: {tid}")
        require(gid not in genes, f"Multiple coding transcripts per gene are unsupported: {gid}")
        genes.add(gid)
        rows.append(dict(transcript_id=tid, gene_id=gid, cds_start=cds_start, cds_end=cds_end,
                         cds_length=cds_end - cds_start, strand=strand))
    require(rows, "No coding transcripts found in GTF")
    return rows


def write_features(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_features(path):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == FIELDS, "Invalid transcript feature header")
        rows = list(reader)
    result = {}
    for row in rows:
        tid = row["transcript_id"]
        require(tid not in result, f"Duplicate transcript feature: {tid}")
        converted = dict(row)
        for key in ("cds_start", "cds_end", "cds_length"):
            converted[key] = int(converted[key])
        require(converted["cds_end"] - converted["cds_start"] == converted["cds_length"],
                f"Invalid CDS length: {tid}")
        result[tid] = converted
    require(result, "Empty transcript features")
    return result
