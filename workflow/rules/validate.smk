rule validate_metadata:
    input:
        samples=config["samples"],
        manifest=config["source_manifest"],
        sources=preflight["source_files"],
        code=[
            "src/__init__.py",
            "src/validation/__init__.py",
            "src/validation/validate.py",
            "workflow/Snakefile",
            "workflow/rules/validate.smk",
            "workflow/schemas/config.schema.json",
            "workflow/schemas/samples.schema.json",
            "envs/workflow.yaml",
            "envs/requirements.txt",
        ]
    output:
        VALIDATION_REPORT
    params:
        python=sys.executable,
        effective_config=json.dumps(config, sort_keys=True)
    threads: 1
    resources:
        mem_mb=config["resources"]["mem_mb"]
    log:
        str(Path(config["paths"]["logs"]) / "validate_metadata.log")
    shell:
        "{params.python:q} -m src.validation.validate --config-json {params.effective_config:q} "
        "--output {output:q} > {log:q} 2>&1"
