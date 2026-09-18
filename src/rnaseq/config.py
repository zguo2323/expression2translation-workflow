from collections import Counter
from pathlib import Path
import yaml
from src.validation.validate import require, relative_path
from src.validation.reference import validate_rna_reference


def load_rna(config, samples):
    path = relative_path(Path.cwd(), config.get("rna_config", "config/rnaseq.yaml"))
    rna = yaml.safe_load(path.read_text())
    expected = {"reference_manifest", "sample_ids", "synthetic", "library_type", "library_evidence", "index_k", "min_mapping_rate", "min_library_compatibility",
                "run_deseq2", "design_confirmed", "design_evidence", "contrast", "min_count", "min_samples", "alpha", "fit_type"}
    require(isinstance(rna, dict) and set(rna) == expected, "Invalid RNA config keys")
    ids = rna["sample_ids"]
    require(isinstance(ids, list) and ids and all(isinstance(s, str) for s in ids) and len(ids) == len(set(ids)), "Invalid RNA sample_ids")
    chosen = {s["sample_id"]: s for s in samples if s["sample_id"] in ids}
    require(set(chosen) == set(ids) and all(s["assay"] == "rnaseq" and s["layout"] == "PAIRED" for s in chosen.values()), "RNA selection requires known paired-end RNA samples")
    require(rna["library_type"] in ("A", "IU", "ISF", "ISR") and isinstance(rna["library_evidence"], str) and bool(rna["library_evidence"].strip()), "RNA library_type/evidence unresolved; choose documented auto-detection or verified orientation")
    require(type(rna["index_k"]) is int and 3 <= rna["index_k"] <= 31 and rna["index_k"] % 2 == 1, "Salmon index_k must be odd, 3..31")
    for key in ("synthetic", "run_deseq2", "design_confirmed"):
        require(type(rna[key]) is bool, f"RNA {key} must be boolean")
    for key in ("min_count", "min_samples"):
        require(type(rna[key]) is int and rna[key] >= 1, f"RNA {key} must be positive integer")
    require(0 < rna["alpha"] < 1 and 0 <= rna["min_mapping_rate"] <= 1, "Invalid RNA alpha/mapping rate")
    require(0 <= rna["min_library_compatibility"] <= 1, "Invalid library compatibility threshold")
    require(rna["fit_type"] in ("parametric", "local", "mean"), "Invalid DESeq2 fit_type")
    contrast = rna["contrast"]
    require(isinstance(contrast, list) and len(contrast) == 3 and contrast[0] == "condition" and contrast[1] != contrast[2], "RNA contrast must be [condition, numerator, denominator]")
    groups = Counter(s["condition"] for s in chosen.values())
    require(set(groups) == set(contrast[1:]), "RNA selected conditions must match contrast")
    if rna["run_deseq2"]:
        require(rna["design_confirmed"] and isinstance(rna["design_evidence"], str) and rna["design_evidence"].strip(), "DESeq2 design not confirmed")
        require(all(n >= 2 for n in groups.values()), "DESeq2 requires >=2 libraries per condition")
    manifest = validate_rna_reference(rna["reference_manifest"])
    reference_id = config["reference"]["reference_id"]
    require(reference_id == manifest["reference_id"], "RNA reference_id disagrees with project reference")
    for role in ("genome", "annotation", "transcriptome", "rrna"):
        require(config["reference"][role] == manifest["files"][role]["path"], f"Project reference {role} disagrees with manifest")
    require(config["reference"]["provider"] == manifest["provider"] and config["reference"]["release"] == manifest["release"], "Project reference provider/release mismatch")
    return rna, chosen, manifest
