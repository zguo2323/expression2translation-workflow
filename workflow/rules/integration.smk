INTEGRATION_CONFIG = config.get("integration_config", "config/integration.yaml")
INTEGRATION_CODE = ["src/integration/run.py", "src/integration/config.py", "src/catalog/build.py", "src/report/build.py"]
INTEGRATION_SETTINGS = json.dumps(INTEGRATION, sort_keys=True)

rule integrate_expression_translation:
    input:
        rna=RNA_REPORT,
        ribo=RIBO_REPORT,
        config=INTEGRATION_CONFIG,
        manifest=RNA["reference_manifest"],
        code=INTEGRATION_CODE
    output:
        directory(INTEGRATION_ROOT)
    params:
        settings=INTEGRATION_SETTINGS,
        task=json.dumps(dict(directory=INTEGRATION_ROOT, config=INTEGRATION,
            config_file=INTEGRATION_CONFIG, reference_manifest=RNA["reference_manifest"],
            reference_id=RNA_REFERENCE["reference_id"],
            rna_tpm=f"{RNA_REPORT}/gene_tpm.tsv", rna_counts=f"{RNA_REPORT}/gene_estimated_counts.tsv",
            ribo_counts=f"{RIBO_REPORT}/gene_counts.tsv",
            rna_samples=[dict(sample_id=sid, condition=RNA_SAMPLES[sid]["condition"]) for sid in RNA["sample_ids"]],
            ribo_samples=[dict(sample_id=sid, condition=RIBO_SAMPLES[sid]["condition"]) for sid in RIBO["sample_ids"]]))
    conda: "../../envs/workflow.yaml"
    log: config["paths"]["logs"] + "/integration/integrate.log"
    shell: "python -m src.integration.run --task {params.task:q} > {log:q} 2>&1"

rule build_catalog:
    input:
        integration=INTEGRATION_ROOT,
        rna=RNA_REPORT,
        ribo=RIBO_REPORT,
        samples=config["samples"],
        manifest=RNA["reference_manifest"],
        annotation=RNA_REFERENCE["files"]["annotation"]["path"],
        config=INTEGRATION_CONFIG,
        code=INTEGRATION_CODE
    output:
        database=CATALOG,
        provenance=CATALOG_PROVENANCE
    params:
        settings=INTEGRATION_SETTINGS,
        task=json.dumps(dict(database=CATALOG, provenance=CATALOG_PROVENANCE,
            integration_analysis=f"{INTEGRATION_ROOT}/analysis.json", integration_config=INTEGRATION_CONFIG,
            reference_manifest=RNA["reference_manifest"], annotation=RNA_REFERENCE["files"]["annotation"]["path"],
            rna_counts=f"{RNA_REPORT}/gene_estimated_counts.tsv", rna_tpm=f"{RNA_REPORT}/gene_tpm.tsv",
            ribo_counts=f"{RIBO_REPORT}/gene_counts.tsv", ribo_tpm=f"{INTEGRATION_ROOT}/ribo_sample_cds_tpm.tsv",
            condition_te=f"{INTEGRATION_ROOT}/condition_te.tsv", te_contrast=f"{INTEGRATION_ROOT}/te_contrast.tsv",
            ribo_qc=f"{RIBO_REPORT}/qc_metrics.tsv",
            rna_qc=[dict(sample_id=sid, path=f"{RNA_ROOT}/quant/{sid}/aux_info/meta_info.json") for sid in RNA["sample_ids"]],
            samples=[RNA_SAMPLES[sid] for sid in RNA["sample_ids"]] + [RIBO_SAMPLES[sid] for sid in RIBO["sample_ids"]],
            artifacts=[
                dict(path=f"{RNA_REPORT}/analysis.json", kind="rna_analysis"),
                dict(path=f"{RNA_REPORT}/provenance.json", kind="rna_provenance"),
                dict(path=f"{RIBO_REPORT}/qc_metrics.tsv", kind="ribo_qc"),
                dict(path=f"{RIBO_REPORT}/provenance.json", kind="ribo_provenance"),
                dict(path=f"{INTEGRATION_ROOT}/analysis.json", kind="integration_analysis"),
                dict(path=f"{INTEGRATION_ROOT}/condition_te.tsv", kind="condition_te"),
                dict(path=f"{INTEGRATION_ROOT}/te_contrast.tsv", kind="te_contrast")]))
    conda: "../../envs/workflow.yaml"
    log: config["paths"]["logs"] + "/integration/catalog.log"
    shell: "python -m src.catalog.build --task {params.task:q} > {log:q} 2>&1"

rule render_report:
    input:
        database=CATALOG,
        provenance=CATALOG_PROVENANCE,
        integration=INTEGRATION_ROOT,
        config=INTEGRATION_CONFIG,
        code=INTEGRATION_CODE
    output:
        markdown=REPORT_ROOT + "/report.md",
        html=FINAL_REPORT,
        plot=REPORT_ROOT + "/te_changes.svg",
        manifest=REPORT_ROOT + "/report_manifest.json"
    params:
        settings=INTEGRATION_SETTINGS,
        task=json.dumps(dict(database=CATALOG, catalog_provenance=CATALOG_PROVENANCE,
            integration_analysis=f"{INTEGRATION_ROOT}/analysis.json", directory=REPORT_ROOT, config=INTEGRATION))
    conda: "../../envs/workflow.yaml"
    log: config["paths"]["logs"] + "/integration/report.log"
    shell: "python -m src.report.build --task {params.task:q} > {log:q} 2>&1"
