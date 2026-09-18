RIBO_CONFIG = config.get("ribo_config", "config/riboseq.yaml")
RIBO_CODE = ["src/riboseq/run.py", "src/riboseq/config.py", "src/riboseq/features.py",
             "src/validation/reference.py", "src/validation/sequences.py"]
RIBO_REF = RIBO_REFERENCE["files"]
RIBO_SETTINGS = {key: RIBO[key] for key in ("strand_mode", "strand_evidence", "offset_table", "offset_evidence",
    "min_mapq", "require_unique", "min_mapping_rate", "max_rrna_fraction", "min_frame0_fraction",
    "metagene_upstream", "metagene_downstream", "synthetic")}

rule ribo_rrna_index:
    input:
        reference=RIBO_REF["rrna"]["path"],
        manifest=RIBO["reference_manifest"],
        code=RIBO_CODE
    output:
        directory(RIBO_ROOT + "/reference/rrna_index")
    threads: config["resources"]["threads"]
    params:
        settings=json.dumps(dict(role="rrna", synthetic=RIBO["synthetic"])),
        task=lambda wc, threads: json.dumps(dict(action="index", role="rrna", reference=RIBO_REF["rrna"]["path"],
            directory=RIBO_ROOT + "/reference/rrna_index", threads=threads, synthetic=RIBO["synthetic"]))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/rrna_index.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"

rule ribo_transcriptome_index:
    input:
        reference=RIBO_REF["transcriptome"]["path"],
        manifest=RIBO["reference_manifest"],
        code=RIBO_CODE
    output:
        directory(RIBO_ROOT + "/reference/transcriptome_index")
    threads: config["resources"]["threads"]
    params:
        settings=json.dumps(dict(role="transcriptome", synthetic=RIBO["synthetic"])),
        task=lambda wc, threads: json.dumps(dict(action="index", role="transcriptome", reference=RIBO_REF["transcriptome"]["path"],
            directory=RIBO_ROOT + "/reference/transcriptome_index", threads=threads, synthetic=RIBO["synthetic"]))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/transcriptome_index.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"

rule ribo_transcript_features:
    input:
        annotation=RIBO_REF["annotation"]["path"],
        transcriptome=RIBO_REF["transcriptome"]["path"],
        manifest=RIBO["reference_manifest"],
        code=RIBO_CODE
    output:
        table=RIBO_ROOT + "/reference/transcript_features.tsv",
        report=RIBO_ROOT + "/reference/transcript_features.json"
    params:
        task=lambda wc, output, input: json.dumps(dict(action="features", annotation=str(input.annotation),
            transcriptome=str(input.transcriptome), output=str(output.table), report=str(output.report)))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/transcript_features.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"

rule ribo_decontaminate:
    input:
        processed=QC_ROOT + "/fastp/{sample}",
        gate=QC_ROOT + "/validation/{sample}.trimmed.json",
        index=RIBO_ROOT + "/reference/rrna_index",
        code=RIBO_CODE
    output:
        directory(RIBO_ROOT + "/decontaminated/{sample}")
    threads: config["resources"]["threads"]
    params:
        settings=json.dumps(dict(max_rrna_fraction=RIBO["max_rrna_fraction"], synthetic=RIBO["synthetic"])),
        task=lambda wc, threads: json.dumps(dict(action="decontaminate", sample_id=wc.sample,
            reads=trim_reads(wc.sample)[0], index=RIBO_ROOT + "/reference/rrna_index/index",
            directory=f"{RIBO_ROOT}/decontaminated/{wc.sample}", threads=threads,
            max_rrna_fraction=RIBO["max_rrna_fraction"], synthetic=RIBO["synthetic"]))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/{sample}.decontaminate.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"

rule ribo_align:
    input:
        reads=RIBO_ROOT + "/decontaminated/{sample}",
        index=RIBO_ROOT + "/reference/transcriptome_index",
        code=RIBO_CODE
    output:
        directory(RIBO_ROOT + "/alignment/{sample}")
    threads: config["resources"]["threads"]
    params:
        settings=json.dumps(dict(min_mapping_rate=RIBO["min_mapping_rate"], synthetic=RIBO["synthetic"])),
        task=lambda wc, threads: json.dumps(dict(action="align", sample_id=wc.sample,
            reads=f"{RIBO_ROOT}/decontaminated/{wc.sample}/clean.fastq.gz",
            index=RIBO_ROOT + "/reference/transcriptome_index/index", directory=f"{RIBO_ROOT}/alignment/{wc.sample}",
            threads=threads, min_mapping_rate=RIBO["min_mapping_rate"], synthetic=RIBO["synthetic"]))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/{sample}.align.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"

rule ribo_count_sites:
    input:
        alignment=RIBO_ROOT + "/alignment/{sample}",
        features=RIBO_ROOT + "/reference/transcript_features.tsv",
        offsets=RIBO["offset_table"],
        config=RIBO_CONFIG,
        code=RIBO_CODE
    output:
        directory(RIBO_ROOT + "/sites/{sample}")
    params:
        settings=json.dumps(RIBO_SETTINGS, sort_keys=True),
        task=lambda wc: json.dumps(dict(action="count", sample_id=wc.sample,
            bam=f"{RIBO_ROOT}/alignment/{wc.sample}/alignment.bam", features=RIBO_ROOT + "/reference/transcript_features.tsv",
            directory=f"{RIBO_ROOT}/sites/{wc.sample}", **RIBO_SETTINGS))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/{sample}.sites.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"

rule ribo_summary:
    input:
        sites=[f"{RIBO_ROOT}/sites/{sid}" for sid in RIBO["sample_ids"]],
        config=RIBO_CONFIG,
        manifest=RIBO["reference_manifest"]
    output:
        directory(RIBO_REPORT)
    params:
        settings=json.dumps(RIBO_SETTINGS, sort_keys=True),
        task=json.dumps(dict(action="summarize", directory=RIBO_REPORT, reference_id=RIBO_REFERENCE["reference_id"],
            reference_manifest=RIBO["reference_manifest"], config_file=RIBO_CONFIG,
            synthetic=RIBO["synthetic"], config=RIBO_SETTINGS,
            samples=[dict(sample_id=sid, condition=RIBO_SAMPLES[sid]["condition"],
                decontam_metrics=f"{RIBO_ROOT}/decontaminated/{sid}/metrics.json",
                alignment_metrics=f"{RIBO_ROOT}/alignment/{sid}/metrics.json",
                site_metrics=f"{RIBO_ROOT}/sites/{sid}/metrics.json",
                gene_counts=f"{RIBO_ROOT}/sites/{sid}/gene_counts.tsv") for sid in RIBO["sample_ids"]]))
    conda: "../../envs/riboseq.yaml"
    log: config["paths"]["logs"] + "/riboseq/summary.log"
    shell: "python -m src.riboseq.run --task {params.task:q} > {log:q} 2>&1"
