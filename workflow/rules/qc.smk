QC_CODE = ["src/qc/run.py", "src/qc/config.py", "src/validation/reads.py", "src/validation/validate.py"]
QC_CONFIG = config.get("qc_config", "config/qc.yaml")
QC_READ_KEYS = {sid: [k for k in ("fastq_1", "fastq_2") if row[k]] for sid, row in QC_SAMPLES.items()}
wildcard_constraints:
    sample="|".join(QC_SAMPLES),
    read="fastq_[12]"
RAW_ZIPS = [f"{QC_ROOT}/fastqc/raw/{sid}/{sid}_{key}_raw_fastqc.zip"
            for sid, keys in QC_READ_KEYS.items() for key in keys]
TRIM_ZIPS = [f"{QC_ROOT}/fastqc/trimmed/{sid}/{sid}_{key}_trimmed_fastqc.zip"
             for sid, keys in QC_READ_KEYS.items() for key in keys]

def reads_for(sid):
    return [QC_SAMPLES[sid][key] for key in QC_READ_KEYS[sid]]

def trim_reads(sid):
    return [f"{QC_ROOT}/fastp/{sid}/{key}.fastq.gz" for key in QC_READ_KEYS[sid]]

rule validate_fastq:
    input:
        reads=lambda wc: reads_for(wc.sample),
        metadata=VALIDATION_REPORT,
        code=QC_CODE
    output:
        QC_ROOT + "/validation/{sample}.json"
    params:
        python=sys.executable
    log:
        config["paths"]["logs"] + "/qc/{sample}.validation.log"
    shell:
        "{params.python:q} -m src.validation.reads --reads {input.reads:q} --output {output:q} > {log:q} 2>&1"

rule fastqc_raw:
    input:
        read=lambda wc: QC_SAMPLES[wc.sample][wc.read],
        gate=QC_ROOT + "/validation/{sample}.json",
        code=QC_CODE
    output:
        zip=QC_ROOT + "/fastqc/raw/{sample}/{sample}_{read}_raw_fastqc.zip",
        html=QC_ROOT + "/fastqc/raw/{sample}/{sample}_{read}_raw_fastqc.html",
        provenance=QC_ROOT + "/fastqc/raw/{sample}/{sample}_{read}_raw.provenance.json"
    params:
        synthetic=QC["synthetic"],
        task=lambda wc, input: json.dumps(dict(tool="fastqc", reads=[str(input.read)],
            directory=f"{QC_ROOT}/fastqc/raw/{wc.sample}", label=f"{wc.sample}_{wc.read}_raw", synthetic=QC["synthetic"]))
    conda:
        "../../envs/qc.yaml"
    log:
        config["paths"]["logs"] + "/qc/{sample}.{read}.raw.log"
    shell:
        "python -m src.qc.run --task {params.task:q} > {log:q} 2>&1"

if QC["mode"] == "trim":
    rule fastp:
        input:
            reads=lambda wc: reads_for(wc.sample),
            gate=QC_ROOT + "/validation/{sample}.json",
            code=QC_CODE,
            policy=QC_CONFIG
        output:
            directory(QC_ROOT + "/fastp/{sample}")
        threads: min(config["resources"]["threads"], 16)
        params:
            task=lambda wc, threads: json.dumps(dict(tool="fastp", assay=QC_SAMPLES[wc.sample]["assay"],
                reads=reads_for(wc.sample), outputs=trim_reads(wc.sample), directory=f"{QC_ROOT}/fastp/{wc.sample}",
                label=wc.sample, threads=threads, policy=QC["policies"][QC_SAMPLES[wc.sample]["assay"]], synthetic=QC["synthetic"]))
        conda:
            "../../envs/qc.yaml"
        log:
            config["paths"]["logs"] + "/qc/{sample}.fastp.log"
        shell:
            "python -m src.qc.run --task {params.task:q} > {log:q} 2>&1"

    rule validate_trimmed:
        input:
            processed=QC_ROOT + "/fastp/{sample}",
            code=QC_CODE
        output:
            QC_ROOT + "/validation/{sample}.trimmed.json"
        params:
            python=sys.executable,
            reads=lambda wc: trim_reads(wc.sample)
        log:
            config["paths"]["logs"] + "/qc/{sample}.trimmed.validation.log"
        shell:
            "{params.python:q} -m src.validation.reads --reads {params.reads:q} --output {output:q} > {log:q} 2>&1"

    rule fastqc_trimmed:
        input:
            gate=QC_ROOT + "/validation/{sample}.trimmed.json",
            processed=QC_ROOT + "/fastp/{sample}",
            code=QC_CODE
        output:
            zip=QC_ROOT + "/fastqc/trimmed/{sample}/{sample}_{read}_trimmed_fastqc.zip",
            html=QC_ROOT + "/fastqc/trimmed/{sample}/{sample}_{read}_trimmed_fastqc.html",
            provenance=QC_ROOT + "/fastqc/trimmed/{sample}/{sample}_{read}_trimmed.provenance.json"
        params:
            task=lambda wc: json.dumps(dict(tool="fastqc", reads=[f"{QC_ROOT}/fastp/{wc.sample}/{wc.read}.fastq.gz"],
                directory=f"{QC_ROOT}/fastqc/trimmed/{wc.sample}", label=f"{wc.sample}_{wc.read}_trimmed", synthetic=QC["synthetic"]))
        conda:
            "../../envs/qc.yaml"
        log:
            config["paths"]["logs"] + "/qc/{sample}.{read}.trimmed.log"
        shell:
            "python -m src.qc.run --task {params.task:q} > {log:q} 2>&1"

rule multiqc:
    input:
        reports=RAW_ZIPS + (TRIM_ZIPS + [f"{QC_ROOT}/fastp/{sid}" for sid in QC_SAMPLES] if QC["mode"] == "trim" else []),
        code=QC_CODE
    output:
        html=QC_REPORT,
        data=directory(QC_ROOT + "/multiqc/multiqc_report_data"),
        provenance=QC_ROOT + "/multiqc/e2t_qc.provenance.json"
    params:
        task=lambda wc, input: json.dumps(dict(tool="multiqc", inputs=list(input.reports), directory=f"{QC_ROOT}/multiqc",
            label="e2t_qc", synthetic=QC["synthetic"])),
        synthetic=QC["synthetic"]
    conda:
        "../../envs/qc.yaml"
    log:
        config["paths"]["logs"] + "/qc/multiqc.log"
    shell:
        "python -m src.qc.run --task {params.task:q} > {log:q} 2>&1"
