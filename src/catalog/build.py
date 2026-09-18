"""Build an integrity-checked SQLite catalog from RNA, Ribo and integration results."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3

from src.integration.run import read_matrix
from src.validation.sequences import open_text
from src.validation.validate import require, sha256


SCHEMA = """
PRAGMA foreign_keys = ON;
CREATE TABLE IF NOT EXISTS "references" (
  reference_id TEXT PRIMARY KEY, provider TEXT NOT NULL, release TEXT NOT NULL,
  manifest_sha256 TEXT NOT NULL CHECK(length(manifest_sha256)=64)
);
CREATE TABLE IF NOT EXISTS samples (
  sample_id TEXT PRIMARY KEY, run_accession TEXT NOT NULL UNIQUE, geo_accession TEXT NOT NULL,
  biosample TEXT NOT NULL, assay TEXT NOT NULL CHECK(assay IN ('rnaseq','riboseq')),
  condition TEXT NOT NULL, replicate INTEGER NOT NULL CHECK(replicate>0),
  layout TEXT NOT NULL CHECK(layout IN ('PAIRED','SINGLE'))
);
CREATE TABLE IF NOT EXISTS workflow_runs (
  run_id TEXT PRIMARY KEY, generated_at TEXT NOT NULL, reference_id TEXT NOT NULL,
  config_sha256 TEXT NOT NULL CHECK(length(config_sha256)=64), synthetic INTEGER NOT NULL CHECK(synthetic IN (0,1)),
  analysis_level TEXT NOT NULL CHECK(analysis_level='descriptive_condition'),
  FOREIGN KEY(reference_id) REFERENCES "references"(reference_id)
);
CREATE TABLE IF NOT EXISTS run_samples (
  run_id TEXT NOT NULL, sample_id TEXT NOT NULL, PRIMARY KEY(run_id,sample_id),
  FOREIGN KEY(run_id) REFERENCES workflow_runs(run_id) ON DELETE CASCADE,
  FOREIGN KEY(sample_id) REFERENCES samples(sample_id)
);
CREATE TABLE IF NOT EXISTS genes (
  reference_id TEXT NOT NULL, gene_id TEXT NOT NULL, contig TEXT NOT NULL,
  start INTEGER NOT NULL CHECK(start>=0), end INTEGER NOT NULL CHECK(end>start),
  strand TEXT NOT NULL CHECK(strand IN ('+','-')), biotype TEXT NOT NULL,
  PRIMARY KEY(reference_id,gene_id), FOREIGN KEY(reference_id) REFERENCES "references"(reference_id)
);
CREATE TABLE IF NOT EXISTS qc_metrics (
  run_id TEXT NOT NULL, sample_id TEXT NOT NULL, stage TEXT NOT NULL,
  total_reads INTEGER, retained_reads INTEGER, mapping_rate REAL, rrna_fraction REAL, frame0_fraction REAL,
  PRIMARY KEY(run_id,sample_id,stage),
  FOREIGN KEY(run_id,sample_id) REFERENCES run_samples(run_id,sample_id) ON DELETE CASCADE
);
CREATE TABLE IF NOT EXISTS rna_gene_metrics (
  run_id TEXT NOT NULL, sample_id TEXT NOT NULL, reference_id TEXT NOT NULL, gene_id TEXT NOT NULL,
  estimated_count REAL NOT NULL CHECK(estimated_count>=0), gene_tpm REAL NOT NULL CHECK(gene_tpm>=0),
  PRIMARY KEY(run_id,sample_id,gene_id),
  FOREIGN KEY(run_id,sample_id) REFERENCES run_samples(run_id,sample_id) ON DELETE CASCADE,
  FOREIGN KEY(reference_id,gene_id) REFERENCES genes(reference_id,gene_id)
);
CREATE TABLE IF NOT EXISTS ribo_gene_metrics (
  run_id TEXT NOT NULL, sample_id TEXT NOT NULL, reference_id TEXT NOT NULL, gene_id TEXT NOT NULL,
  p_site_count REAL NOT NULL CHECK(p_site_count>=0), cds_length_nt INTEGER NOT NULL CHECK(cds_length_nt>0),
  cds_tpm REAL NOT NULL CHECK(cds_tpm>=0),
  PRIMARY KEY(run_id,sample_id,gene_id),
  FOREIGN KEY(run_id,sample_id) REFERENCES run_samples(run_id,sample_id) ON DELETE CASCADE,
  FOREIGN KEY(reference_id,gene_id) REFERENCES genes(reference_id,gene_id)
);
CREATE TABLE IF NOT EXISTS integration_metrics (
  run_id TEXT NOT NULL, reference_id TEXT NOT NULL, condition TEXT NOT NULL, gene_id TEXT NOT NULL,
  rna_tpm_mean REAL NOT NULL CHECK(rna_tpm_mean>=0), ribo_cds_tpm_mean REAL NOT NULL CHECK(ribo_cds_tpm_mean>=0),
  n_rna_samples INTEGER NOT NULL CHECK(n_rna_samples>0), n_ribo_samples INTEGER NOT NULL CHECK(n_ribo_samples>0),
  eligible INTEGER NOT NULL CHECK(eligible IN (0,1)), filter_reason TEXT NOT NULL, te_log2 REAL,
  PRIMARY KEY(run_id,condition,gene_id), FOREIGN KEY(run_id) REFERENCES workflow_runs(run_id) ON DELETE CASCADE,
  FOREIGN KEY(reference_id,gene_id) REFERENCES genes(reference_id,gene_id)
);
CREATE TABLE IF NOT EXISTS integration_contrasts (
  run_id TEXT NOT NULL, reference_id TEXT NOT NULL, contrast TEXT NOT NULL, gene_id TEXT NOT NULL,
  eligible INTEGER NOT NULL CHECK(eligible IN (0,1)), filter_reason TEXT NOT NULL,
  rna_log2_change REAL, ribo_log2_change REAL, te_log2_change REAL,
  PRIMARY KEY(run_id,contrast,gene_id), FOREIGN KEY(run_id) REFERENCES workflow_runs(run_id) ON DELETE CASCADE,
  FOREIGN KEY(reference_id,gene_id) REFERENCES genes(reference_id,gene_id)
);
CREATE TABLE IF NOT EXISTS artifacts (
  run_id TEXT NOT NULL, path TEXT NOT NULL, kind TEXT NOT NULL,
  sha256 TEXT NOT NULL CHECK(length(sha256)=64), bytes INTEGER NOT NULL CHECK(bytes>=0),
  PRIMARY KEY(run_id,path), FOREIGN KEY(run_id) REFERENCES workflow_runs(run_id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_integration_te ON integration_metrics(run_id, condition, eligible, te_log2);
CREATE INDEX IF NOT EXISTS idx_contrast_te ON integration_contrasts(run_id, contrast, eligible, te_log2_change);
"""


def gene_features(annotation):
    genes = {}
    with open_text(annotation) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            require(len(fields) == 9, f"GTF line {number}: expected 9 columns")
            contig, _, feature, start, end, _, strand, _, attributes = fields
            if feature not in ("gene", "exon"):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', attributes))
            require("gene_id" in attrs and strand in ("+", "-"), f"Invalid gene annotation at line {number}")
            gene = attrs["gene_id"]
            biotype = attrs.get("gene_biotype", attrs.get("gene_type", "unknown"))
            start, end = int(start) - 1, int(end)
            if gene in genes:
                old = genes[gene]
                require((old["contig"], old["strand"], old["biotype"]) == (contig, strand, biotype),
                        f"Conflicting gene annotation: {gene}")
                old["start"], old["end"] = min(old["start"], start), max(old["end"], end)
            else:
                genes[gene] = dict(gene_id=gene, contig=contig, start=start, end=end, strand=strand, biotype=biotype)
    require(genes, "No genes found in annotation")
    return genes


def table_rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def nullable(value):
    return None if value == "" else float(value)


def build(task):
    analysis = json.loads(Path(task["integration_analysis"]).read_text())
    manifest = json.loads(Path(task["reference_manifest"]).read_text())
    require(analysis["reference_id"] == manifest["reference_id"], "Integration/reference identity mismatch")
    identity = json.dumps(dict(reference_id=analysis["reference_id"], config=analysis["config"],
                               inputs=analysis["inputs_sha256"]), sort_keys=True, separators=(",", ":"))
    run_id = "e2t-" + hashlib.sha256(identity.encode()).hexdigest()[:20]
    config_sha = sha256(Path(task["integration_config"]))
    selected = {row["sample_id"]: row for row in task["samples"]}
    expected_ids = analysis["rna_samples"] + analysis["ribo_samples"]
    require(set(selected) == set(expected_ids), "Catalog selected samples disagree with integration")
    genes = gene_features(task["annotation"])
    reference_id = manifest["reference_id"]
    rna_counts = read_matrix(task["rna_counts"], analysis["rna_samples"])
    rna_tpm = read_matrix(task["rna_tpm"], analysis["rna_samples"])
    ribo_rows = table_rows(task["ribo_counts"])
    ribo_tpm = read_matrix(task["ribo_tpm"], analysis["ribo_samples"])
    condition_rows = table_rows(task["condition_te"])
    contrast_rows = table_rows(task["te_contrast"])
    all_metric_genes = set(rna_counts) | {row["gene_id"] for row in ribo_rows}
    require(all_metric_genes <= set(genes), "Quantification contains genes absent from GTF")

    database = Path(task["database"])
    database.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA foreign_keys=ON")
    try:
        with connection:
            connection.executescript(SCHEMA)
            connection.execute("INSERT INTO \"references\" VALUES (?,?,?,?) ON CONFLICT(reference_id) DO UPDATE SET provider=excluded.provider, release=excluded.release, manifest_sha256=excluded.manifest_sha256",
                               (reference_id, manifest["provider"], manifest["release"], sha256(Path(task["reference_manifest"]))))
            for row in selected.values():
                connection.execute("INSERT INTO samples VALUES (?,?,?,?,?,?,?,?) ON CONFLICT(sample_id) DO UPDATE SET run_accession=excluded.run_accession, geo_accession=excluded.geo_accession, biosample=excluded.biosample, assay=excluded.assay, condition=excluded.condition, replicate=excluded.replicate, layout=excluded.layout",
                                   (row["sample_id"], row["run_accession"], row["geo_accession"], row["biosample"], row["assay"], row["condition"], int(row["replicate"]), row["layout"]))
            connection.execute("INSERT INTO workflow_runs VALUES (?,?,?,?,?,?) ON CONFLICT(run_id) DO UPDATE SET generated_at=excluded.generated_at",
                               (run_id, datetime.now(timezone.utc).isoformat(), reference_id, config_sha, int(analysis["synthetic"]), "descriptive_condition"))
            for sample in expected_ids:
                connection.execute("INSERT OR IGNORE INTO run_samples VALUES (?,?)", (run_id, sample))
            for gene in genes.values():
                connection.execute("INSERT INTO genes VALUES (?,?,?,?,?,?,?) ON CONFLICT(reference_id,gene_id) DO UPDATE SET contig=excluded.contig,start=excluded.start,end=excluded.end,strand=excluded.strand,biotype=excluded.biotype",
                                   (reference_id, gene["gene_id"], gene["contig"], gene["start"], gene["end"], gene["strand"], gene["biotype"]))
            for entry in task["rna_qc"]:
                meta = json.loads(Path(entry["path"]).read_text())
                connection.execute("INSERT OR REPLACE INTO qc_metrics VALUES (?,?,?,?,?,?,?,?)",
                                   (run_id, entry["sample_id"], "salmon", int(meta["num_processed"]), int(meta["num_mapped"]),
                                    float(meta["percent_mapped"]) / 100, None, None))
            for row in table_rows(task["ribo_qc"]):
                connection.execute("INSERT OR REPLACE INTO qc_metrics VALUES (?,?,?,?,?,?,?,?)",
                                   (run_id, row["sample_id"], "riboseq", int(row["input_reads"]), int(row["assigned_cds_psites"]),
                                    float(row["mapping_rate"]), float(row["rrna_fraction"]), float(row["frame0_fraction"])))
            for gene, values in rna_counts.items():
                for sample in analysis["rna_samples"]:
                    connection.execute("INSERT OR REPLACE INTO rna_gene_metrics VALUES (?,?,?,?,?,?)",
                                       (run_id, sample, reference_id, gene, values[sample], rna_tpm[gene][sample]))
            for row in ribo_rows:
                gene, sample = row["gene_id"], row["sample_id"]
                connection.execute("INSERT OR REPLACE INTO ribo_gene_metrics VALUES (?,?,?,?,?,?,?)",
                                   (run_id, sample, reference_id, gene, float(row["p_site_count"]), int(row["cds_length_nt"]), ribo_tpm[gene][sample]))
            for row in condition_rows:
                connection.execute("INSERT OR REPLACE INTO integration_metrics VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    (run_id, reference_id, row["condition"], row["gene_id"], float(row["rna_tpm_mean"]),
                     float(row["ribo_cds_tpm_mean"]), int(row["n_rna_samples"]), int(row["n_ribo_samples"]),
                     int(row["eligible"] == "true"), row["filter_reason"], nullable(row["te_log2"])))
            for row in contrast_rows:
                connection.execute("INSERT OR REPLACE INTO integration_contrasts VALUES (?,?,?,?,?,?,?,?,?)",
                    (run_id, reference_id, row["contrast"], row["gene_id"], int(row["eligible"] == "true"), row["filter_reason"],
                     nullable(row["rna_log2_change"]), nullable(row["ribo_log2_change"]), nullable(row["te_log2_change"])))
            for artifact in task["artifacts"]:
                path = Path(artifact["path"])
                require(path.is_file(), f"Missing catalog artifact: {path}")
                connection.execute("INSERT OR REPLACE INTO artifacts VALUES (?,?,?,?,?)",
                                   (run_id, str(path), artifact["kind"], sha256(path), path.stat().st_size))
        foreign = connection.execute("PRAGMA foreign_key_check").fetchall()
        quick = connection.execute("PRAGMA quick_check").fetchone()[0]
        require(not foreign and quick == "ok", f"SQLite integrity failed: {foreign or quick}")
        counts = {table: connection.execute(f"SELECT count(*) FROM {table} WHERE run_id=?", (run_id,)).fetchone()[0]
                  for table in ("run_samples", "qc_metrics", "rna_gene_metrics", "ribo_gene_metrics",
                                "integration_metrics", "integration_contrasts", "artifacts")}
    finally:
        connection.close()
    provenance = dict(run_id=run_id, generated_at=datetime.now(timezone.utc).isoformat(), reference_id=reference_id,
                      synthetic=analysis["synthetic"], config_sha256=config_sha, database=str(database), counts=counts,
                      sqlite_integrity="ok", inputs_sha256={key: sha256(Path(task[key])) for key in
                          ("integration_analysis", "reference_manifest", "annotation", "rna_counts", "rna_tpm",
                           "ribo_counts", "ribo_tpm", "condition_te", "te_contrast", "ribo_qc")},
                      rna_qc_sha256={entry["sample_id"]: sha256(Path(entry["path"])) for entry in task["rna_qc"]},
                      artifact_inputs_sha256={entry["path"]: sha256(Path(entry["path"])) for entry in task["artifacts"]})
    target = Path(task["provenance"])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    build(json.loads(parser.parse_args().task))
