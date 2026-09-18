import csv
from contextlib import chdir
import copy
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

import yaml

from src.catalog.build import build as build_catalog
from src.integration.run import integrate
from src.integration.config import load_integration
from src.report.build import build as build_report
from src.validation.validate import ValidationError


class IntegrationCatalogTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.rna_ids = ["young_rna", "middle_rna"]
        self.ribo_ids = ["young_ribo", "middle_ribo"]
        self.samples = [
            dict(sample_id="young_rna", run_accession="SRR1", geo_accession="GSM1", biosample="SAM1", assay="rnaseq", condition="Young", replicate="1", layout="PAIRED"),
            dict(sample_id="middle_rna", run_accession="SRR2", geo_accession="GSM2", biosample="SAM2", assay="rnaseq", condition="Middle", replicate="1", layout="PAIRED"),
            dict(sample_id="young_ribo", run_accession="SRR3", geo_accession="GSM3", biosample="SAM3", assay="riboseq", condition="Young", replicate="1", layout="SINGLE"),
            dict(sample_id="middle_ribo", run_accession="SRR4", geo_accession="GSM4", biosample="SAM4", assay="riboseq", condition="Middle", replicate="1", layout="SINGLE"),
        ]
        self.config = dict(synthetic=True, normalization="rna_gene_tpm__ribo_cds_tpm",
            aggregation="arithmetic_mean_by_condition", pseudocount=1.0, min_rna_tpm=1.0,
            min_ribo_tpm=1.0, contrast=["Middle", "Young"], differential_te=False, report_title="E2T test")
        (self.root / "integration.yaml").write_text(yaml.safe_dump(self.config, sort_keys=False))
        (self.root / "rna_tpm.tsv").write_text("gene_id\tyoung_rna\tmiddle_rna\ng1\t100\t100\ng2\t900\t900\ng3\t0\t0\n")
        (self.root / "rna_counts.tsv").write_text("gene_id\tyoung_rna\tmiddle_rna\ng1\t10\t10\ng2\t90\t90\ng3\t0\t0\n")
        with (self.root / "ribo_counts.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "condition", "gene_id", "p_site_count", "cds_length_nt", "p_sites_per_kb", "rpm"])
            for sample, condition, values in (("young_ribo", "Young", [10, 90, 0]), ("middle_ribo", "Middle", [20, 80, 0])):
                for gene, count in zip(("g1", "g2", "g3"), values):
                    writer.writerow([sample, condition, gene, count, 300, count / .3, count * 10000])
        (self.root / "manifest.json").write_text(json.dumps(dict(reference_id="ref1", provider="test", release="v1")))
        self.integration = self.root / "integration"
        integrate(dict(directory=str(self.integration), config=self.config, config_file=str(self.root / "integration.yaml"),
            reference_manifest=str(self.root / "manifest.json"), reference_id="ref1",
            rna_tpm=str(self.root / "rna_tpm.tsv"), rna_counts=str(self.root / "rna_counts.tsv"),
            ribo_counts=str(self.root / "ribo_counts.tsv"),
            rna_samples=[dict(sample_id="young_rna", condition="Young"), dict(sample_id="middle_rna", condition="Middle")],
            ribo_samples=[dict(sample_id="young_ribo", condition="Young"), dict(sample_id="middle_ribo", condition="Middle")]))

    def catalog_task(self):
        annotation = self.root / "annotation.gtf"
        annotation.write_text("".join(f'chr\ttest\tgene\t{index*400+1}\t{index*400+300}\t.\t+\t.\tgene_id "{gene}"; gene_biotype "protein_coding";\n'
                                      for index, gene in enumerate(("g1", "g2", "g3"))))
        for sample in self.rna_ids:
            (self.root / f"{sample}.meta.json").write_text(json.dumps(dict(num_processed=100, num_mapped=90, percent_mapped=90)))
        (self.root / "ribo_qc.tsv").write_text("sample_id\tcondition\tinput_reads\trrna_fraction\tmapping_rate\tassigned_cds_psites\tframe0_fraction\n"
            "young_ribo\tYoung\t110\t0.0909\t1.0\t100\t0.9\n"
            "middle_ribo\tMiddle\t110\t0.0909\t1.0\t100\t0.9\n")
        for name in ("rna_analysis.json", "rna_provenance.json", "ribo_provenance.json"):
            (self.root / name).write_text("{}\n")
        return dict(database=str(self.root / "catalog.db"), provenance=str(self.root / "catalog.json"),
            integration_analysis=str(self.integration / "analysis.json"), integration_config=str(self.root / "integration.yaml"),
            reference_manifest=str(self.root / "manifest.json"), annotation=str(annotation),
            rna_counts=str(self.root / "rna_counts.tsv"), rna_tpm=str(self.root / "rna_tpm.tsv"),
            ribo_counts=str(self.root / "ribo_counts.tsv"), ribo_tpm=str(self.integration / "ribo_sample_cds_tpm.tsv"),
            condition_te=str(self.integration / "condition_te.tsv"), te_contrast=str(self.integration / "te_contrast.tsv"),
            ribo_qc=str(self.root / "ribo_qc.tsv"), samples=self.samples,
            rna_qc=[dict(sample_id=sample, path=str(self.root / f"{sample}.meta.json")) for sample in self.rna_ids],
            artifacts=[dict(path=str(self.integration / "analysis.json"), kind="integration_analysis"),
                       dict(path=str(self.root / "rna_analysis.json"), kind="rna_analysis"),
                       dict(path=str(self.root / "rna_provenance.json"), kind="rna_provenance"),
                       dict(path=str(self.root / "ribo_provenance.json"), kind="ribo_provenance")])

    def test_normalization_filtering_and_contrast(self):
        with (self.integration / "condition_te.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        g1 = {(row["gene_id"], row["condition"]): row for row in rows}
        self.assertAlmostEqual(float(g1[("g1", "Young")]["rna_tpm_mean"]), 100000)
        self.assertAlmostEqual(float(g1[("g1", "Young")]["ribo_cds_tpm_mean"]), 100000)
        self.assertEqual(g1[("g3", "Young")]["eligible"], "false")
        with (self.integration / "te_contrast.tsv").open() as handle:
            contrasts = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertAlmostEqual(float(contrasts["g1"]["te_log2_change"]), 1.0, places=4)
        self.assertEqual(contrasts["g3"]["te_log2_change"], "")

    def test_integration_config_blocks_inference_and_design_mismatch(self):
        rna = dict(synthetic=True, reference_manifest="manifest.json")
        ribo = dict(synthetic=True, reference_manifest="manifest.json")
        rna_samples = {row["sample_id"]: row for row in self.samples if row["assay"] == "rnaseq"}
        ribo_samples = {row["sample_id"]: row for row in self.samples if row["assay"] == "riboseq"}
        reference = {"reference_id": "ref1"}
        with chdir(self.root):
            loaded = load_integration({"integration_config": "integration.yaml"}, rna, rna_samples, ribo, ribo_samples, reference, reference)
        self.assertEqual(loaded["contrast"], ["Middle", "Young"])
        changed = copy.deepcopy(self.config)
        changed["differential_te"] = True
        (self.root / "integration.yaml").write_text(yaml.safe_dump(changed, sort_keys=False))
        with chdir(self.root), self.assertRaisesRegex(ValidationError, "descriptive"):
            load_integration({"integration_config": "integration.yaml"}, rna, rna_samples, ribo, ribo_samples, reference, reference)
        changed["differential_te"] = False
        (self.root / "integration.yaml").write_text(yaml.safe_dump(changed, sort_keys=False))
        ribo_samples["middle_ribo"]["condition"] = "Other"
        with chdir(self.root), self.assertRaisesRegex(ValidationError, "same configured conditions"):
            load_integration({"integration_config": "integration.yaml"}, rna, rna_samples, ribo, ribo_samples, reference, reference)

    def test_matrix_sample_identity_is_strict(self):
        path = self.root / "bad_order.tsv"
        path.write_text("gene_id\tmiddle_rna\tyoung_rna\ng1\t1\t1\n")
        from src.integration.run import read_matrix
        with self.assertRaisesRegex(ValidationError, "samples/order"):
            read_matrix(path, self.rna_ids)

    def test_catalog_integrity_idempotency_and_reports(self):
        task = self.catalog_task()
        build_catalog(task)
        build_catalog(task)
        connection = sqlite3.connect(task["database"])
        try:
            self.assertEqual(connection.execute("PRAGMA foreign_key_check").fetchall(), [])
            self.assertEqual(connection.execute("SELECT count(*) FROM workflow_runs").fetchone()[0], 1)
            self.assertEqual(connection.execute("SELECT count(*) FROM run_samples").fetchone()[0], 4)
            self.assertEqual(connection.execute("SELECT count(*) FROM integration_metrics").fetchone()[0], 6)
            self.assertEqual(connection.execute("SELECT count(*) FROM integration_contrasts WHERE eligible=1").fetchone()[0], 2)
        finally:
            connection.close()
        report_dir = self.root / "report"
        build_report(dict(database=task["database"], catalog_provenance=task["provenance"],
            integration_analysis=str(self.integration / "analysis.json"), directory=str(report_dir), config=self.config))
        markdown = (report_dir / "report.md").read_text()
        html = (report_dir / "report.html").read_text()
        self.assertIn("descriptive", markdown)
        self.assertIn("SYNTHETIC TEST", html)
        self.assertTrue((report_dir / "report_manifest.json").is_file())
        self.assertTrue((report_dir / "te_changes.svg").read_text().startswith("<svg"))


if __name__ == "__main__":
    unittest.main()
