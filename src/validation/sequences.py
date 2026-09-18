"""Strict yeast-sized genome/GTF parser and deterministic exon reconstruction."""
from collections import defaultdict
import gzip
from pathlib import Path
import re
from src.validation.validate import require


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta(path):
    sequences = {}
    name = None
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                require(name not in sequences, f"Duplicate FASTA ID: {name}")
                sequences[name] = []
            else:
                require(name is not None and re.fullmatch("[ACGTNacgtn]+", line), f"Invalid FASTA sequence: {path}")
                sequences[name].append(line.upper())
    result = {k: "".join(v) for k, v in sequences.items()}
    require(result and all(result.values()), f"Empty FASTA record/file: {path}")
    return result


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def reconstruct(genome_path, gtf_path):
    genome = fasta(genome_path)
    transcripts = {}
    with open_text(gtf_path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            require(len(fields) == 9, f"GTF line {number}: expected 9 columns")
            contig, _, feature, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)
            require(contig in genome, f"GTF contig absent from genome: {contig}")
            require(1 <= start <= end <= len(genome[contig]), f"GTF coordinates out of range: line {number}")
            if feature != "exon":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', attributes))
            require("transcript_id" in attrs and "gene_id" in attrs, f"Exon missing gene/transcript ID: {number}")
            require(strand in ("+", "-"), f"Invalid exon strand: {number}")
            tid, gid = attrs["transcript_id"], attrs["gene_id"]
            require(not any(c.isspace() for c in tid + gid), "Whitespace in gene/transcript ID")
            biotype = attrs.get("transcript_biotype", attrs.get("gene_biotype", attrs.get("gene_type", "unknown")))
            identity = (gid, contig, strand, biotype)
            entry = transcripts.setdefault(tid, {"identity": identity, "exons": []})
            require(entry["identity"] == identity, f"Conflicting transcript mapping: {tid}")
            entry["exons"].append((start - 1, end))
    require(transcripts, "No exons found in GTF")
    result = {}
    for tid, item in sorted(transcripts.items()):
        gid, contig, strand, biotype = item["identity"]
        exons = sorted(item["exons"])
        require(all(left[1] <= right[0] for left, right in zip(exons, exons[1:])), f"Overlapping/duplicate exons: {tid}")
        seq = "".join(genome[contig][start:end] for start, end in exons)
        if strand == "-":
            seq = reverse_complement(seq)
        result[tid] = dict(gene_id=gid, contig=contig, strand=strand, biotype=biotype, sequence=seq)
    require(not set(genome).intersection(result), "Transcript IDs collide with genome decoy IDs")
    return genome, result


def write_fasta(path, sequences):
    with Path(path).open("w") as handle:
        for name, seq in sequences.items():
            handle.write(f">{name}\n")
            for offset in range(0, len(seq), 80):
                handle.write(seq[offset:offset+80] + "\n")
