"""Independent synthetic sample designs shared by real-tool smoke scripts."""
import csv
from pathlib import Path


def write_samples(target, config, assays=("rnaseq", "riboseq"), conditions=("Young", "Middle"), replicates=("1", "2")):
    rows = []
    for assay in assays:
        for condition in conditions:
            for rep in replicates:
                sid = f"{condition.lower()}_{'rna' if assay == 'rnaseq' else 'ribo'}_{rep}"
                paired = assay == "rnaseq"
                rows.append(dict(sample_id=sid, assay=assay, condition=condition, replicate=rep,
                                 layout="PAIRED" if paired else "SINGLE",
                                 fastq_1=f"data/raw/synthetic_{sid}{'_1' if paired else ''}.fastq.gz",
                                 fastq_2=f"data/raw/synthetic_{sid}_2.fastq.gz" if paired else ""))
    config.update(source_validation="local", samples="config/samples.tsv", dataset={"organism": "Saccharomyces cerevisiae"})
    for key in ("source_manifest", "sources", "study_samples"):
        config.pop(key, None)
    config["design"] = dict(conditions=list(conditions), assays={"rnaseq": "PAIRED", "riboseq": "SINGLE"},
                            pairing="unconfirmed", integration_level="condition")
    with (Path(target) / config["samples"]).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows
