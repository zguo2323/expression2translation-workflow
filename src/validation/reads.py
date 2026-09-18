"""Streaming validation of gzip FASTQ, including paired-end order and IDs."""
import argparse
import gzip
from itertools import zip_longest
import json
from pathlib import Path
import re

from src.validation.validate import ValidationError, require, sha256


def records(path):
    with gzip.open(path, "rt", encoding="ascii") as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            number += 1
            seq, plus, qual = (handle.readline().rstrip("\r\n") for _ in range(3))
            require(header.startswith("@") and len(header.strip()) > 1 and plus.startswith("+"),
                    f"{path}: invalid FASTQ record {number}")
            require(bool(seq) and len(seq) == len(qual), f"{path}: sequence/quality length mismatch at {number}")
            require(re.fullmatch("[ACGTNacgtn]+", seq), f"{path}: invalid sequence alphabet at {number}")
            require(all(33 <= ord(c) <= 126 for c in qual), f"{path}: invalid quality at {number}")
            identifier = header[1:].split()[0]
            identifier = re.sub(r"/[12]$", "", identifier)
            yield identifier, len(seq)


def validate_reads(paths):
    paths = list(map(Path, paths))
    require(len(paths) in (1, 2), "Expected one SE or two PE FASTQ files")
    require(len(set(p.resolve() for p in paths)) == len(paths), "Duplicate FASTQ inputs")
    counts = [0] * len(paths)
    bases = [0] * len(paths)
    try:
        for pair in zip_longest(*(records(p) for p in paths)):
            require(all(item is not None for item in pair), "Paired FASTQ read counts differ")
            require(len(paths) == 1 or pair[0][0] == pair[1][0], "Paired FASTQ IDs/order differ")
            for index, (_, length) in enumerate(pair):
                counts[index] += 1
                bases[index] += length
        require(all(counts), "Empty FASTQ file")
        return {"status": "passed", "layout": "PAIRED" if len(paths) == 2 else "SINGLE",
                "files": [{"path": str(p), "reads": c, "bases": b, "sha256": sha256(p), "bytes": p.stat().st_size}
                          for p, c, b in zip(paths, counts, bases)]}
    except (OSError, EOFError, UnicodeError) as error:
        raise ValidationError(f"FASTQ/gzip validation failed: {error}") from error


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = validate_reads(args.reads)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + "\n")
