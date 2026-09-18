from pathlib import Path
import re
import yaml

from src.validation.validate import require, relative_path


def load_qc(config, samples, sample_ids=None):
    path = relative_path(Path.cwd(), config.get("qc_config", "config/qc.yaml"))
    qc = yaml.safe_load(path.read_text())
    require(isinstance(qc, dict) and set(qc) == {"mode", "sample_ids", "synthetic", "policies"}, "Invalid QC config keys")
    require(qc["mode"] in ("raw", "trim"), "QC mode must be raw or trim")
    require(type(qc["synthetic"]) is bool, "synthetic must be boolean")
    ids = qc["sample_ids"] if sample_ids is None else sample_ids
    require(isinstance(ids, list) and all(isinstance(s, str) for s in ids) and len(ids) == len(set(ids)), "Invalid QC sample_ids")
    selected = {r["sample_id"]: r for r in samples if not ids or r["sample_id"] in ids}
    require(not ids or set(ids) == set(selected), "QC sample_ids contains unknown sample")
    require(set(qc["policies"]) == {"rnaseq", "riboseq"}, "Expected policies for both assays")
    if qc["mode"] == "trim":
        for assay in {r["assay"] for r in selected.values()}:
            validate_policy(qc["policies"][assay], assay)
    for row in selected.values():
        for key in ("fastq_1", "fastq_2"):
            if row[key]:
                require(Path(row[key]).is_file(), f"Missing FASTQ for {row['sample_id']}: {row[key]}; choose downloaded sample_ids or keep stage=metadata")
    return qc, selected


def validate_policy(policy, assay):
    keys = {"confirmed", "evidence", "adapter_mode", "adapter_r1", "adapter_r2", "barcode_mode",
            "trim_front1", "trim_front2", "min_length", "quality_phred", "unqualified_percent"}
    require(isinstance(policy, dict) and set(policy) == keys, f"{assay}: invalid preprocessing policy keys")
    require(policy["confirmed"] is True and isinstance(policy["evidence"], str) and policy["evidence"].strip(),
            f"{assay}: preprocessing policy needs confirmed=true and evidence")
    require(policy["adapter_mode"] in ("sequence", "none"), f"{assay}: adapter strategy unresolved")
    for key in ("adapter_r1", "adapter_r2"):
        value = policy[key]
        require(value is None or (isinstance(value, str) and re.fullmatch("[ACGT]+", value)), f"{assay}: invalid {key}")
    if policy["adapter_mode"] == "sequence":
        require(bool(policy["adapter_r1"]) and (assay != "rnaseq" or bool(policy["adapter_r2"])), f"{assay}: explicit adapter sequences required")
    else:
        require(policy["adapter_r1"] is None and policy["adapter_r2"] is None, f"{assay}: adapter_mode=none conflicts with sequences")
    require(policy["barcode_mode"] in ("none", "fixed_front"), f"{assay}: barcode sorting/strategy unsupported or unresolved")
    for key in ("trim_front1", "trim_front2"):
        require(type(policy[key]) is int and policy[key] >= 0, f"{assay}: invalid {key}")
    require(policy["barcode_mode"] != "none" or policy["trim_front1"] == policy["trim_front2"] == 0,
            f"{assay}: barcode_mode=none conflicts with trimming")
    require(assay != "riboseq" or (policy["trim_front2"] == 0 and policy["adapter_r2"] is None), "Ribo SINGLE cannot configure read2")
    for key, low, high in (("min_length", 1, 10000), ("quality_phred", 0, 93), ("unqualified_percent", 0, 100)):
        require(type(policy[key]) is int and low <= policy[key] <= high, f"{assay}: invalid/unresolved {key}")
