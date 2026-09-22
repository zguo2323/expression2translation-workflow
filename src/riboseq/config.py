"""Strict configuration loader for the Ribo-seq branch."""
import csv
from pathlib import Path
import yaml

from src.validation.reference import validate_rna_reference
from src.validation.validate import require, relative_path, select_samples


KEYS = {"reference_manifest", "sample_ids", "synthetic", "strand_mode", "strand_evidence",
        "footprint_lengths", "offset_table", "offset_evidence", "min_mapq", "require_unique", "min_mapping_rate",
        "max_rrna_fraction", "min_frame0_fraction", "metagene_upstream", "metagene_downstream"}


def load_offsets(path):
    offsets = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["read_length", "psite_offset"], "Offset table header must be read_length, psite_offset")
        for row in reader:
            try:
                length, offset = int(row["read_length"]), int(row["psite_offset"])
            except (TypeError, ValueError) as error:
                raise ValueError("Offset table values must be integers") from error
            require(length > 0 and 0 <= offset < length, f"Invalid P-site offset for length {length}")
            require(length not in offsets, f"Duplicate P-site read length: {length}")
            offsets[length] = offset
    require(offsets, "Empty P-site offset table")
    return offsets


def load_ribo(config, samples):
    path = relative_path(Path.cwd(), config.get("ribo_config", "config/riboseq.yaml"))
    ribo = yaml.safe_load(path.read_text())
    require(isinstance(ribo, dict) and set(ribo) == KEYS, "Invalid Ribo config keys")
    require(type(ribo["synthetic"]) is bool, "Ribo synthetic must be boolean")
    ribo["sample_ids"], selected = select_samples(ribo["sample_ids"], samples, "riboseq", "SINGLE")
    require(ribo["strand_mode"] in ("forward", "reverse", "unstranded") and isinstance(ribo["strand_evidence"], str)
            and ribo["strand_evidence"].strip(), "Ribo strand_mode/evidence unresolved")
    require(isinstance(ribo["offset_evidence"], str) and ribo["offset_evidence"].strip(), "Ribo offset evidence unresolved")
    require(isinstance(ribo["offset_table"], str) and ribo["offset_table"].strip(), "Ribo offset table unresolved")
    offset_path = relative_path(Path.cwd(), ribo["offset_table"])
    offsets = load_offsets(offset_path)
    lengths = ribo["footprint_lengths"]
    require(isinstance(lengths, list) and lengths and all(type(length) is int and length > 0 for length in lengths)
            and len(lengths) == len(set(lengths)), "Ribo footprint_lengths must be a non-empty list of unique positive integers")
    require(set(lengths) == set(offsets), "Ribo footprint_lengths must exactly match offset table read lengths")
    for key in ("min_mapq", "metagene_upstream", "metagene_downstream"):
        require(type(ribo[key]) is int and ribo[key] >= 0, f"Invalid Ribo {key}")
    require(type(ribo["require_unique"]) is bool, "require_unique must be boolean")
    for key in ("min_mapping_rate", "max_rrna_fraction", "min_frame0_fraction"):
        require(type(ribo[key]) in (int, float) and 0 <= ribo[key] <= 1, f"Invalid Ribo {key}")
    manifest = validate_rna_reference(ribo["reference_manifest"], Path.cwd())
    return ribo, selected, manifest, offsets
