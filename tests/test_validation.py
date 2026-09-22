"""Failure tests isolate fixture copies; no real reads or network required."""

import copy
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

from src.validation.validate import ValidationError, validate_metadata


ROOT = Path(__file__).resolve().parents[1]


class MetadataValidationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        for directory in ("metadata", "config", "workflow", "src", "envs"):
            shutil.copytree(ROOT / directory, self.root / directory,
                            ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
        self.config = yaml.safe_load((self.root / "config/config.yaml").read_text())

    def edit_samples(self, change):
        path = self.root / self.config["samples"]
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        change(rows)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def fails(self, message):
        with self.assertRaisesRegex(ValidationError, message):
            validate_metadata(self.config, self.root)

    def test_valid_without_fastq_or_reference(self):
        report = validate_metadata(self.config, self.root)
        self.assertEqual(report["sample_count"], 8)
        self.assertEqual(report["expected_fastq_count"], 12)
        self.assertEqual(report["archive_bytes"], 10395980191)
        self.assertEqual(report["bases"], 22453938457)
        self.assertFalse(report["analysis_ready"])
        self.assertFalse((self.root / "results").exists())
        self.assertTrue(all(r["sra_assay_type"] == "OTHER" for r in report["samples"] if r["assay"] == "riboseq"))

    def local_mode(self):
        self.config["source_validation"] = "local"
        for key in ("source_manifest", "sources", "study_samples"):
            self.config.pop(key, None)
        shutil.rmtree(self.root / "metadata")
        def clean(rows):
            for row in rows:
                for key in ("run_accession", "geo_accession", "biosample"):
                    row[key] = ""
        self.edit_samples(clean)

    def test_four_run_subset_preserves_full_study_sources(self):
        source = self.root / self.config["study_samples"]
        before = source.read_bytes()
        ids = {"young_rna_2", "middle_rna_1", "young_ribo_1", "middle_ribo_1"}
        self.edit_samples(lambda rows: rows.__setitem__(slice(None), [r for r in rows if r["sample_id"] in ids]))
        report = validate_metadata(self.config, self.root)
        self.assertEqual(report["sample_count"], 4)
        self.assertEqual(report["archive_bytes"], 4384989091)
        self.assertEqual(before, source.read_bytes())

    def test_study_paths_and_conditions_independent_of_active_subset(self):
        self.config["paths"]["raw"] = "data/custom"
        self.config["design"]["conditions"] = ["Young"]
        def subset(rows):
            rows[:] = [r for r in rows if r["condition"] == "Young"]
            for row in rows:
                for key in ("fastq_1", "fastq_2"):
                    row[key] = row[key].replace("data/raw/", "data/custom/")
        self.edit_samples(subset)
        self.assertEqual(validate_metadata(self.config, self.root)["sample_count"], 4)

    def test_local_optional_columns_absent_and_unknown_sizes(self):
        self.local_mode()
        path = self.root / self.config["samples"]
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        fields = [k for k in rows[0] if k not in ("run_accession", "geo_accession", "biosample")]
        with path.open("w") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        report = validate_metadata(self.config, self.root)
        self.assertEqual(report["sample_count"], 8)
        self.assertIsNone(report["archive_bytes"])
        self.assertIsNone(report["bases"])
        self.assertIsNone(report["source_manifest"])
        self.assertEqual(report["source_files"], [])
        self.config["paths"]["results"] = "results/local"
        result = subprocess.run([sys.executable, "-m", "src.validation.validate", "--config-json", json.dumps(self.config),
                                 "--output", "results/local/validation/metadata.json"],
                                cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        snapshot = json.loads((self.root / "results/local/validation/metadata.json").read_text())
        self.assertIsNone(snapshot["provenance"]["source_manifest_sha256"])

    def test_local_variable_counts_replicate_labels_and_single_condition(self):
        self.local_mode()
        self.config["design"]["conditions"] = ["Young"]
        def subset(rows):
            rows[:] = [r for r in rows if r["condition"] == "Young"][:3]
            for index, row in enumerate(rows):
                row["replicate"] = f"donor_{index + 7}"
        self.edit_samples(subset)
        report = validate_metadata(self.config, self.root)
        self.assertEqual(report["sample_count"], 3)

    def test_local_invalid_accession_and_duplicate_biosample(self):
        self.local_mode()
        self.edit_samples(lambda rows: rows[0].update(run_accession="OTHER"))
        self.fails("run_accession")
        self.edit_samples(lambda rows: rows[0].update(run_accession=""))
        self.edit_samples(lambda rows: [r.update(biosample="SAMN1") for r in rows[:2]])
        self.fails("Duplicate BioSample")

    def test_study_mapping_still_checked_against_geo(self):
        active = self.config["samples"]
        self.config["samples"] = self.config["study_samples"]
        def swap(rows):
            rows[0]["replicate"], rows[1]["replicate"] = rows[1]["replicate"], rows[0]["replicate"]
        self.edit_samples(swap)
        self.config["samples"] = active
        self.fails("GEO title")

    @unittest.skipUnless(shutil.which("snakemake"), "Snakemake needed for local DAG test")
    def test_local_dag_without_source_files(self):
        self.local_mode()
        (self.root / "config/config.yaml").write_text(yaml.safe_dump(self.config))
        result = subprocess.run([shutil.which("snakemake"), "--snakefile", "workflow/Snakefile", "--cores", "1", "--dry-run"],
                                cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("validate_metadata", result.stdout + result.stderr)

    def test_duplicate_run(self):
        self.edit_samples(lambda rows: rows[1].update(run_accession=rows[0]["run_accession"]))
        self.fails("Duplicate run_accession")

    def test_duplicate_sample(self):
        self.edit_samples(lambda rows: rows[1].update(sample_id=rows[0]["sample_id"]))
        self.fails("Duplicate sample_id")

    def test_wrong_layout(self):
        self.edit_samples(lambda rows: rows[0].update(layout="SINGLE"))
        self.fails("assay/layout mismatch")

    def test_missing_replicate(self):
        self.edit_samples(lambda rows: rows[1].update(replicate="1"))
        self.fails("Duplicate biological identity")

    def test_wrong_accession(self):
        self.edit_samples(lambda rows: rows[0].update(run_accession="SRR99999999"))
        self.fails("study run_accession mismatch")

    def test_wrong_geo_mapping(self):
        def swap(rows):
            rows[0]["geo_accession"], rows[1]["geo_accession"] = rows[1]["geo_accession"], rows[0]["geo_accession"]
        self.edit_samples(swap)
        self.fails("study geo_accession mismatch")

    def test_wrong_biosample(self):
        self.edit_samples(lambda rows: rows[0].update(biosample="SAMN99999999"))
        self.fails("study biosample mismatch")

    def test_other_is_not_project_assay(self):
        self.edit_samples(lambda rows: rows[4].update(assay="OTHER"))
        self.fails("assay")

    def test_paired_missing_read2(self):
        self.edit_samples(lambda rows: rows[0].update(fastq_2=""))
        self.fails("PAIRED requires fastq_2")

    def test_single_with_read2(self):
        self.edit_samples(lambda rows: rows[4].update(fastq_2="data/raw/extra.fastq.gz"))
        self.fails("SINGLE requires empty fastq_2")

    def test_duplicate_fastq(self):
        self.edit_samples(lambda rows: rows[1].update(fastq_1=rows[0]["fastq_1"]))
        self.fails("duplicate FASTQ path")

    def test_missing_source(self):
        (self.root / self.config["sources"]["sra"]).unlink()
        self.fails("Missing source file")

    def test_changed_source(self):
        with (self.root / self.config["sources"]["sra"]).open("a") as handle:
            handle.write("\n")
        self.fails("checksum/size mismatch")

    def test_invalid_config(self):
        self.config["resources"]["threads"] = 0
        self.fails("resources.threads")

    def test_paired_te_cannot_be_enabled(self):
        self.config["integration"]["differential_te"] = True
        self.fails("differential_te")

    def test_source_crosscheck_catches_wrong_replicate(self):
        def swap(rows):
            rows[0]["replicate"], rows[1]["replicate"] = rows[1]["replicate"], rows[0]["replicate"]
        self.edit_samples(swap)
        self.fails("study replicate mismatch")

    def test_path_traversal(self):
        self.edit_samples(lambda rows: rows[0].update(fastq_1="../outside.fastq.gz"))
        self.fails("repository-relative")

    def test_output_cannot_overwrite_metadata(self):
        self.config["paths"]["results"] = "metadata"
        self.fails("protected directory")

    def test_config_override_snapshot(self):
        config = copy.deepcopy(self.config)
        config["resources"]["threads"] = 2
        result = subprocess.run([sys.executable, "-m", "src.validation.validate",
                                 "--config-json", json.dumps(config),
                                 "--output", "results/validation/metadata.json"],
                                cwd=self.root, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads((self.root / "results/validation/metadata.json").read_text())
        self.assertEqual(report["provenance"]["effective_config"], config)

    @unittest.skipUnless(shutil.which("snakemake"), "Snakemake needed for DAG integration test")
    def test_snakemake_dry_run_and_invalid_sheet(self):
        command = [shutil.which("snakemake"), "--snakefile", "workflow/Snakefile", "--cores", "1", "--dry-run"]
        good = subprocess.run(command, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(good.returncode, 0, good.stdout + good.stderr)
        self.assertIn("validate_metadata", good.stdout + good.stderr)
        self.assertFalse((self.root / "results/validation/metadata.json").exists())
        self.edit_samples(lambda rows: rows[1].update(run_accession=rows[0]["run_accession"]))
        bad = subprocess.run(command, cwd=self.root, capture_output=True, text=True)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("Duplicate run_accession", bad.stdout + bad.stderr)
        self.assertFalse((self.root / "results/validation/metadata.json").exists())


if __name__ == "__main__":
    unittest.main()
