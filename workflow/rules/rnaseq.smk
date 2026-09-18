RNA_CONFIG = config.get("rna_config", "config/rnaseq.yaml")
RNA_CODE = ["src/rnaseq/run.py", "src/rnaseq/config.py", "src/validation/reference.py", "src/validation/sequences.py"]
RNA_REF_FILES = [entry["path"] for entry in RNA_REFERENCE["files"].values()]

rule salmon_index:
    input:
        manifest=RNA["reference_manifest"],
        reference=RNA_REF_FILES,
        code=RNA_CODE
    output:
        directory(RNA_ROOT + "/salmon_index")
    threads: config["resources"]["threads"]
    params:
        settings=json.dumps(dict(index_k=RNA["index_k"], synthetic=RNA["synthetic"])),
        task=lambda wc, threads: json.dumps(dict(action="index", manifest=RNA["reference_manifest"],
            directory=RNA_ROOT + "/salmon_index", k=RNA["index_k"], threads=threads, synthetic=RNA["synthetic"]))
    conda: "../../envs/salmon.yaml"
    log: config["paths"]["logs"] + "/rnaseq/index.log"
    shell: "python -m src.rnaseq.run --task {params.task:q} > {log:q} 2>&1"

rule salmon_quant:
    input:
        processed=QC_ROOT + "/fastp/{sample}",
        gate=QC_ROOT + "/validation/{sample}.trimmed.json",
        index=RNA_ROOT + "/salmon_index",
        manifest=RNA["reference_manifest"],
        code=RNA_CODE
    output:
        directory(RNA_ROOT + "/quant/{sample}")
    threads: config["resources"]["threads"]
    params:
        settings=json.dumps({k: RNA[k] for k in ("library_type", "library_evidence", "min_mapping_rate", "min_library_compatibility", "synthetic")}),
        task=lambda wc, threads: json.dumps(dict(action="quant", directory=f"{RNA_ROOT}/quant/{wc.sample}",
            index=RNA_ROOT + "/salmon_index/index", reads=trim_reads(wc.sample), threads=threads,
            library_type=RNA["library_type"], library_evidence=RNA["library_evidence"], manifest=RNA["reference_manifest"],
            min_mapping_rate=RNA["min_mapping_rate"], min_library_compatibility=RNA["min_library_compatibility"], synthetic=RNA["synthetic"]))
    conda: "../../envs/salmon.yaml"
    log: config["paths"]["logs"] + "/rnaseq/{sample}.salmon.log"
    shell: "python -m src.rnaseq.run --task {params.task:q} > {log:q} 2>&1"

rule rnaseq_gene_analysis:
    input:
        quant=[f"{RNA_ROOT}/quant/{sid}" for sid in RNA["sample_ids"]],
        tx2gene=RNA_REFERENCE["files"]["tx2gene"]["path"],
        script="src/rnaseq/gene_analysis.R",
        provenance="src/rnaseq/provenance.py",
        config=RNA_CONFIG
    output:
        directory(RNA_REPORT)
    params:
        python=sys.executable,
        task=json.dumps(dict(config=RNA, reference_id=RNA_REFERENCE["reference_id"],
            tx2gene=RNA_REFERENCE["files"]["tx2gene"]["path"], directory=RNA_REPORT,
            samples=[dict(sample_id=sid, condition=RNA_SAMPLES[sid]["condition"],
                          quant=f"{RNA_ROOT}/quant/{sid}/quant.sf") for sid in RNA["sample_ids"]]))
    conda: "../../envs/rnaseq-stats.yaml"
    log: config["paths"]["logs"] + "/rnaseq/gene_analysis.log"
    shell:
        "Rscript --vanilla {input.script:q} {params.task:q} > {log:q} 2>&1 && "
        "{params.python:q} -m src.rnaseq.provenance --task {params.task:q} >> {log:q} 2>&1"
