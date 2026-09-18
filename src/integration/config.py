"""Configuration and design checks for descriptive integration."""
from collections import Counter
from pathlib import Path
import yaml

from src.validation.validate import require, relative_path


KEYS = {"synthetic", "normalization", "aggregation", "pseudocount", "min_rna_tpm", "min_ribo_tpm",
        "contrast", "differential_te", "report_title"}


def load_integration(config, rna, rna_samples, ribo, ribo_samples, rna_reference, ribo_reference):
    path = relative_path(Path.cwd(), config.get("integration_config", "config/integration.yaml"))
    settings = yaml.safe_load(path.read_text())
    require(isinstance(settings, dict) and set(settings) == KEYS, "Invalid integration config keys")
    require(type(settings["synthetic"]) is bool and settings["synthetic"] == rna["synthetic"] == ribo["synthetic"],
            "Integration/RNA/Ribo synthetic flags disagree")
    require(settings["normalization"] == "rna_gene_tpm__ribo_cds_tpm", "Unsupported integration normalization")
    require(settings["aggregation"] == "arithmetic_mean_by_condition", "Unsupported integration aggregation")
    for key in ("pseudocount", "min_rna_tpm", "min_ribo_tpm"):
        require(type(settings[key]) in (int, float) and settings[key] >= 0, f"Invalid integration {key}")
    require(settings["pseudocount"] > 0, "Integration pseudocount must be positive")
    require(settings["differential_te"] is False, "MVP integration is descriptive; differential_te must be false")
    require(isinstance(settings["report_title"], str) and settings["report_title"].strip(), "Missing report title")
    contrast = settings["contrast"]
    require(isinstance(contrast, list) and len(contrast) == 2 and contrast[0] != contrast[1],
            "Integration contrast must be [numerator, denominator]")
    rna_groups = Counter(row["condition"] for row in rna_samples.values())
    ribo_groups = Counter(row["condition"] for row in ribo_samples.values())
    require(set(rna_groups) == set(ribo_groups) == set(contrast),
            "RNA and Ribo selections must cover the same configured conditions")
    require(rna_reference["reference_id"] == ribo_reference["reference_id"], "RNA/Ribo reference IDs disagree")
    require(rna["reference_manifest"] == ribo["reference_manifest"], "RNA/Ribo reference manifests disagree")
    return settings
