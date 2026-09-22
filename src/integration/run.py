"""Create normalized condition-level descriptive TE tables."""
import argparse
from collections import defaultdict
import csv
import itertools
import json
import math
from pathlib import Path

from src.validation.validate import require, sha256


def read_matrix(path, expected_samples):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames and reader.fieldnames[0] == "gene_id", f"Invalid gene matrix header: {path}")
        require(reader.fieldnames[1:] == expected_samples, f"Matrix samples/order disagree: {path}")
        result = {}
        for row in reader:
            gene = row["gene_id"]
            require(gene and gene not in result, f"Duplicate/empty gene in {path}: {gene}")
            values = {}
            for sample in expected_samples:
                value = float(row[sample])
                require(math.isfinite(value) and value >= 0, f"Invalid value for {gene}/{sample}: {path}")
                values[sample] = value
            result[gene] = values
    require(result, f"Empty gene matrix: {path}")
    return result


def normalize_per_million(matrix, samples):
    totals = {sample: sum(row[sample] for row in matrix.values()) for sample in samples}
    require(all(value > 0 for value in totals.values()), "Cannot normalize an empty RNA sample")
    return {gene: {sample: values[sample] * 1_000_000 / totals[sample] for sample in samples}
            for gene, values in matrix.items()}


def read_ribo(path, expected_samples):
    fields = ["sample_id", "condition", "gene_id", "p_site_count", "cds_length_nt", "p_sites_per_kb", "rpm"]
    result, lengths, conditions = defaultdict(dict), {}, {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == fields, "Invalid Ribo gene-count header")
        for row in reader:
            sample, gene = row["sample_id"], row["gene_id"]
            require(sample in expected_samples, f"Unexpected Ribo sample: {sample}")
            require(gene not in result[sample], f"Duplicate Ribo sample/gene: {sample}/{gene}")
            count, length = float(row["p_site_count"]), int(row["cds_length_nt"])
            require(math.isfinite(count) and count >= 0 and length > 0, f"Invalid Ribo metric: {sample}/{gene}")
            require(gene not in lengths or lengths[gene] == length, f"Inconsistent CDS length: {gene}")
            require(sample not in conditions or conditions[sample] == row["condition"], f"Inconsistent condition: {sample}")
            result[sample][gene], lengths[gene], conditions[sample] = count, length, row["condition"]
    require(set(result) == set(expected_samples), "Ribo samples disagree with configuration")
    genes = set(next(iter(result.values())))
    require(genes and all(set(values) == genes for values in result.values()), "Ribo gene sets differ across samples")
    return result, lengths, conditions


def ribo_cds_tpm(counts, lengths):
    result = {sample: {} for sample in counts}
    for sample, values in counts.items():
        rates = {gene: count / (lengths[gene] / 1000) for gene, count in values.items()}
        total = sum(rates.values())
        require(total > 0, f"Cannot normalize empty Ribo sample: {sample}")
        result[sample] = {gene: value * 1_000_000 / total for gene, value in rates.items()}
    return result


def mean(values):
    return sum(values) / len(values)


def descriptive_change(rna_by_condition, ribo_by_condition, cfg, pseudocount=None):
    """Return a condition contrast from already-normalized per-sample values."""
    pseudocount = cfg["pseudocount"] if pseudocount is None else pseudocount
    numerator, denominator = cfg["contrast"]
    values = {}
    for condition in (numerator, denominator):
        rna_mean, ribo_mean = mean(rna_by_condition[condition]), mean(ribo_by_condition[condition])
        if rna_mean < cfg["min_rna_tpm"] or ribo_mean < cfg["min_ribo_tpm"]:
            return None
        values[condition] = (rna_mean, ribo_mean)
    top, base = values[numerator], values[denominator]
    return math.log2((top[1] + pseudocount) / (top[0] + pseudocount)) - \
        math.log2((base[1] + pseudocount) / (base[0] + pseudocount))


def consistent_sign(values):
    """Require every sensitivity value to be available and to retain one sign."""
    return bool(values) and all(value is not None for value in values) and \
        (all(value > 0 for value in values) or all(value < 0 for value in values))


def sensitivity_rows(shared, rna_tpm, ribo_tpm, ribo_counts, rna_samples, ribo_samples,
                     rna_conditions, ribo_conditions, cfg):
    """Expose count support and one-library-removed TE sensitivity per gene."""
    conditions = cfg["contrast"]
    rna_groups = {condition: [sample for sample in rna_samples if rna_conditions[sample] == condition]
                  for condition in conditions}
    ribo_groups = {condition: [sample for sample in ribo_samples if ribo_conditions[sample] == condition]
                   for condition in conditions}
    # A singleton group has no valid leave-one-out result; retain it in the
    # descriptive output but never call it sensitivity-stable.
    removal_sets = [("rna", condition, sample) for condition, samples in rna_groups.items()
                    if len(samples) >= 2 for sample in samples]
    removal_sets += [("ribo", condition, sample) for condition, samples in ribo_groups.items()
                     if len(samples) >= 2 for sample in samples]
    rows = []
    for gene in shared:
        rna_values = {condition: [rna_tpm[gene][sample] for sample in samples]
                      for condition, samples in rna_groups.items()}
        ribo_values = {condition: [ribo_tpm[sample][gene] for sample in samples]
                       for condition, samples in ribo_groups.items()}
        baseline = descriptive_change(rna_values, ribo_values, cfg)
        pseudocount_values = [descriptive_change(rna_values, ribo_values, cfg, value)
                              for value in cfg["candidate_pseudocounts"]]
        leave_one_out = []
        for assay, condition, sample in removal_sets:
            rna_loo = {key: list(values) for key, values in rna_values.items()}
            ribo_loo = {key: list(values) for key, values in ribo_values.items()}
            if assay == "rna":
                rna_loo[condition] = [rna_tpm[gene][sid] for sid in rna_groups[condition] if sid != sample]
            else:
                ribo_loo[condition] = [ribo_tpm[sid][gene] for sid in ribo_groups[condition] if sid != sample]
            leave_one_out.append(descriptive_change(rna_loo, ribo_loo, cfg))
        counts = [ribo_counts[sample][gene] for sample in ribo_samples]
        count_supported = min(counts) >= cfg["candidate_min_ribo_psites_per_sample"]
        stable = consistent_sign(pseudocount_values) and consistent_sign(leave_one_out)
        rows.append(dict(gene_id=gene, eligible=baseline is not None,
                         min_ribo_psites_all_samples=min(counts),
                         ribo_samples_passing_count_gate=sum(count >= cfg["candidate_min_ribo_psites_per_sample"] for count in counts),
                         count_supported=count_supported, te_change=baseline,
                         pseudocount_sign_consistent=consistent_sign(pseudocount_values),
                         loo_sign_consistent=consistent_sign(leave_one_out),
                         candidate_stable=count_supported and stable,
                         te_change_pseudocount_min=min(pseudocount_values) if pseudocount_values and all(value is not None for value in pseudocount_values) else None,
                         te_change_pseudocount_max=max(pseudocount_values) if pseudocount_values and all(value is not None for value in pseudocount_values) else None,
                         te_change_loo_min=min(leave_one_out) if leave_one_out and all(value is not None for value in leave_one_out) else None,
                         te_change_loo_max=max(leave_one_out) if leave_one_out and all(value is not None for value in leave_one_out) else None))
    return rows


def differential_te_readiness(rna_samples, ribo_samples, rna_conditions, ribo_conditions,
                              sensitivity, ribo_qc, cfg):
    """Document why a model is, or is not, safe to add for this particular run."""
    conditions = cfg["contrast"]
    group_sizes = {f"rna:{condition}": sum(rna_conditions[sample] == condition for sample in rna_samples)
                   for condition in conditions}
    group_sizes.update({f"ribo:{condition}": sum(ribo_conditions[sample] == condition for sample in ribo_samples)
                        for condition in conditions})
    assigned = {row["sample_id"]: int(row["assigned_cds_psites"]) for row in ribo_qc}
    rrna = {row["sample_id"]: float(row["rrna_fraction"]) for row in ribo_qc}
    supported = sum(row["eligible"] and row["count_supported"] for row in sensitivity)
    reasons = []
    if min(group_sizes.values()) < 3:
        reasons.append("fewer_than_three_biological_replicates_per_assay_condition")
    if min(assigned.values()) < 500_000:
        reasons.append("one_or_more_ribo_libraries_have_fewer_than_500000_assigned_cds_psites")
    if max(rrna.values()) > 0.90:
        reasons.append("ribo_rrna_fraction_exceeds_90_percent")
    if supported < 1_000:
        reasons.append("fewer_than_1000_genes_pass_the_per_sample_ribo_count_gate")
    return dict(model="assay_by_condition_negative_binomial_interaction", status="not_run",
                recommended=not reasons, conditions=conditions, group_sizes=group_sizes,
                count_supported_genes=supported, candidate_count_gate=cfg["candidate_min_ribo_psites_per_sample"],
                ribo_assigned_cds_psites=assigned, ribo_rrna_fraction=rrna, blocking_reasons=reasons,
                interpretation=("A differential-TE model should use replicate-level raw counts with assay-aware normalization "
                                "and report an assay×condition interaction, not a t-test on condition-level TE."))


def integrate(task):
    output = Path(task["directory"])
    output.mkdir(parents=True, exist_ok=True)
    rna_samples = [row["sample_id"] for row in task["rna_samples"]]
    rna_conditions = {row["sample_id"]: row["condition"] for row in task["rna_samples"]}
    ribo_samples = [row["sample_id"] for row in task["ribo_samples"]]
    ribo_expected_conditions = {row["sample_id"]: row["condition"] for row in task["ribo_samples"]}
    rna_raw_tpm = read_matrix(task["rna_tpm"], rna_samples)
    rna_counts = read_matrix(task["rna_counts"], rna_samples)
    ribo_counts, cds_lengths, ribo_conditions = read_ribo(task["ribo_counts"], ribo_samples)
    require(ribo_conditions == ribo_expected_conditions, "Ribo result conditions disagree with samples.tsv")
    shared = sorted(set(rna_raw_tpm).intersection(cds_lengths))
    require(shared, "RNA and Ribo have no shared gene IDs")
    # Scale both assays over the same coding-gene universe. RNA-only rRNA and
    # other biotypes must not change the denominator used for descriptive TE.
    rna_shared = {gene: rna_raw_tpm[gene] for gene in shared}
    ribo_shared = {sample: {gene: values[gene] for gene in shared} for sample, values in ribo_counts.items()}
    shared_lengths = {gene: cds_lengths[gene] for gene in shared}
    rna_tpm = normalize_per_million(rna_shared, rna_samples)
    ribo_tpm = ribo_cds_tpm(ribo_shared, shared_lengths)
    conditions = task["config"]["contrast"]
    for condition in conditions:
        require(any(value == condition for value in rna_conditions.values()), f"Missing RNA condition: {condition}")
        require(any(value == condition for value in ribo_conditions.values()), f"Missing Ribo condition: {condition}")

    with (output / "rna_sample_tpm.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *rna_samples])
        writer.writerows([gene, *(f"{rna_tpm[gene][sample]:.10f}" for sample in rna_samples)] for gene in shared)
    with (output / "ribo_sample_cds_tpm.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *ribo_samples])
        writer.writerows([gene, *(f"{ribo_tpm[sample][gene]:.10f}" for sample in ribo_samples)] for gene in shared)

    condition_rows, by_gene = [], defaultdict(dict)
    cfg = task["config"]
    for condition in conditions:
        rs = [sample for sample in rna_samples if rna_conditions[sample] == condition]
        bs = [sample for sample in ribo_samples if ribo_conditions[sample] == condition]
        for gene in shared:
            rv = mean([rna_tpm[gene][sample] for sample in rs])
            bv = mean([ribo_tpm[sample][gene] for sample in bs])
            reasons = []
            if rv < cfg["min_rna_tpm"]:
                reasons.append("rna_below_min")
            if bv < cfg["min_ribo_tpm"]:
                reasons.append("ribo_below_min")
            eligible = not reasons
            te = math.log2((bv + cfg["pseudocount"]) / (rv + cfg["pseudocount"])) if eligible else None
            row = dict(gene_id=gene, condition=condition, rna_tpm_mean=rv, ribo_cds_tpm_mean=bv,
                       n_rna_samples=len(rs), n_ribo_samples=len(bs), eligible=eligible,
                       filter_reason=";".join(reasons), te_log2=te)
            condition_rows.append(row)
            by_gene[gene][condition] = row
    with (output / "condition_te.tsv").open("w", newline="") as handle:
        fields = ["gene_id", "condition", "rna_tpm_mean", "ribo_cds_tpm_mean", "n_rna_samples",
                  "n_ribo_samples", "eligible", "filter_reason", "te_log2"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in condition_rows:
            clean = dict(row, rna_tpm_mean=f"{row['rna_tpm_mean']:.10f}", ribo_cds_tpm_mean=f"{row['ribo_cds_tpm_mean']:.10f}",
                         eligible=str(row["eligible"]).lower(), te_log2="" if row["te_log2"] is None else f"{row['te_log2']:.10f}")
            writer.writerow(clean)
    numerator, denominator = conditions
    contrast_rows = []
    for gene in shared:
        top, base = by_gene[gene][numerator], by_gene[gene][denominator]
        eligible = top["eligible"] and base["eligible"]
        contrast_rows.append(dict(gene_id=gene, contrast=f"{numerator}_vs_{denominator}", eligible=eligible,
            filter_reason="" if eligible else "condition_not_eligible",
            rna_log2_change=math.log2((top["rna_tpm_mean"] + cfg["pseudocount"]) / (base["rna_tpm_mean"] + cfg["pseudocount"])) if eligible else None,
            ribo_log2_change=math.log2((top["ribo_cds_tpm_mean"] + cfg["pseudocount"]) / (base["ribo_cds_tpm_mean"] + cfg["pseudocount"])) if eligible else None,
            te_log2_change=top["te_log2"] - base["te_log2"] if eligible else None))
    with (output / "te_contrast.tsv").open("w", newline="") as handle:
        fields = ["gene_id", "contrast", "eligible", "filter_reason", "rna_log2_change", "ribo_log2_change", "te_log2_change"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in contrast_rows:
            writer.writerow({key: (str(value).lower() if isinstance(value, bool) else "" if value is None else
                                   f"{value:.10f}" if isinstance(value, float) else value) for key, value in row.items()})
    sensitivity = sensitivity_rows(shared, rna_tpm, ribo_tpm, ribo_shared, rna_samples, ribo_samples,
                                   rna_conditions, ribo_conditions, cfg)
    with (output / "te_candidate_qc.tsv").open("w", newline="") as handle:
        fields = ["gene_id", "eligible", "min_ribo_psites_all_samples", "ribo_samples_passing_count_gate",
                  "count_supported", "te_change", "pseudocount_sign_consistent", "loo_sign_consistent",
                  "candidate_stable", "te_change_pseudocount_min", "te_change_pseudocount_max",
                  "te_change_loo_min", "te_change_loo_max"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sensitivity:
            writer.writerow({key: (str(value).lower() if isinstance(value, bool) else "" if value is None else
                                   str(int(value)) if key in ("min_ribo_psites_all_samples", "ribo_samples_passing_count_gate") else
                                   f"{value:.10f}" if isinstance(value, float) else value) for key, value in row.items()})
    ribo_qc = []
    with Path(task["ribo_qc"]).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["sample_id", "condition", "input_reads", "rrna_fraction", "mapping_rate",
                                      "assigned_cds_psites", "frame0_fraction"], "Invalid Ribo QC header")
        ribo_qc = list(reader)
    require({row["sample_id"] for row in ribo_qc} == set(ribo_samples), "Ribo QC samples disagree with integration")
    readiness = differential_te_readiness(rna_samples, ribo_samples, rna_conditions, ribo_conditions,
                                          sensitivity, ribo_qc, cfg)
    (output / "differential_te_readiness.json").write_text(json.dumps(readiness, indent=2) + "\n")
    summary = dict(reference_id=task["reference_id"], synthetic=cfg["synthetic"], method="descriptive_condition_level",
                   config=cfg, conditions=conditions, rna_samples=rna_samples, ribo_samples=ribo_samples,
                   rna_input_genes=len(rna_raw_tpm), ribo_input_genes=len(cds_lengths),
                   rna_genes=len(rna_tpm), ribo_genes=len(shared_lengths), shared_genes=len(shared),
                   normalization_universe="shared_coding_genes",
                   eligible_condition_rows=sum(row["eligible"] for row in condition_rows),
                   eligible_contrast_genes=sum(row["eligible"] for row in contrast_rows),
                   candidate_count_supported_genes=sum(row["eligible"] and row["count_supported"] for row in sensitivity),
                   candidate_stable_genes=sum(row["eligible"] and row["candidate_stable"] for row in sensitivity),
                   no_pairing_assumed=True, differential_te=False, differential_te_readiness=readiness,
                   inputs_sha256={key: sha256(Path(task[key])) for key in ("rna_tpm", "rna_counts", "ribo_counts", "ribo_qc", "reference_manifest", "config_file")})
    artifacts = [output / name for name in ("rna_sample_tpm.tsv", "ribo_sample_cds_tpm.tsv", "condition_te.tsv", "te_contrast.tsv",
                                            "te_candidate_qc.tsv", "differential_te_readiness.json")]
    summary["artifacts_sha256"] = {str(path): sha256(path) for path in artifacts}
    (output / "analysis.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    integrate(json.loads(parser.parse_args().task))
