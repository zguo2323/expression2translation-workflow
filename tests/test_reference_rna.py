import copy
import csv
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import yaml

from scripts.prepare_reference import prepare
from src.rnaseq.config import load_rna
from src.validation.reference import validate_rna_reference
from src.validation.sequences import reconstruct
from src.validation.validate import ValidationError, sha256


ROOT = Path(__file__).resolve().parents[1]


class ReferenceSequenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "genome.fa").write_text(">chr1\nAAAACCCCGGGGTTTT\n")
        rows = []
        for tid, gid, strand, biotype, exons in (
            ("tp", "gp", "+", "protein_coding", [(1, 4), (9, 12)]),
            ("tn", "gn", "-", "protein_coding", [(1, 4), (9, 12)]),
            ("tr", "gr", "+", "rRNA", [(5, 8)])):
            for start, end in exons:
                rows.append(f'chr1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; gene_biotype "{biotype}";\n')
        (self.root / "annotation.gtf").write_text("".join(rows))
        shutil.copytree(ROOT / "workflow/schemas", self.root / "workflow/schemas")
        (self.root / "scripts").mkdir()
        shutil.copy(ROOT / "scripts/prepare_reference.py", self.root / "scripts")
        (self.root / "src/validation").mkdir(parents=True)
        shutil.copy(ROOT / "src/validation/sequences.py", self.root / "src/validation")

    def parse(self):
        return reconstruct(self.root / "genome.fa", self.root / "annotation.gtf")

    def prepared(self):
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            prepare("genome.fa", "annotation.gtf", "derived", "test", "synthetic", "1", {"genome": "test", "annotation": "test"})
        finally:
            os.chdir(previous)
        return self.root / "derived/manifest.json"

    def test_positive_negative_splicing_and_rrna(self):
        _, tx = self.parse()
        self.assertEqual(tx["tp"]["sequence"], "AAAAGGGG")
        self.assertEqual(tx["tn"]["sequence"], "CCCCTTTT")
        self.assertEqual(tx["tr"]["sequence"], "CCCC")
        self.prepared()
        result = validate_rna_reference("derived/manifest.json", self.root)
        self.assertIn("tx2gene", result["files"])

    def test_unknown_contig(self):
        path = self.root / "annotation.gtf"
        path.write_text(path.read_text().replace("chr1", "unknown"))
        with self.assertRaisesRegex(ValidationError, "contig absent"):
            self.parse()

    def test_coordinate_overflow(self):
        path = self.root / "annotation.gtf"
        path.write_text(path.read_text().replace("\t9\t12\t", "\t9\t22\t"))
        with self.assertRaisesRegex(ValidationError, "out of range"):
            self.parse()

    def test_overlapping_exons(self):
        path = self.root / "annotation.gtf"
        path.write_text(path.read_text().replace("\t9\t12\t", "\t4\t12\t"))
        with self.assertRaisesRegex(ValidationError, "Overlapping"):
            self.parse()

    def test_conflicting_gene_mapping(self):
        path = self.root / "annotation.gtf"
        path.write_text(path.read_text().replace('gene_id "gp"', 'gene_id "other"', 1))
        with self.assertRaisesRegex(ValidationError, "Conflicting transcript"):
            self.parse()

    def test_sequence_disagreement_even_with_matching_checksum(self):
        manifest_path = self.prepared()
        manifest = json.loads(manifest_path.read_text())
        path = self.root / manifest["files"]["transcriptome"]["path"]
        path.write_text(path.read_text().replace("AAAAGGGG", "AAAAGGGT"))
        manifest["files"]["transcriptome"]["sha256"] = sha256(path)
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValidationError, "Transcriptome sequences"):
            validate_rna_reference("derived/manifest.json", self.root)

    def test_tx2gene_disagreement_even_with_matching_checksum(self):
        manifest_path = self.prepared()
        manifest = json.loads(manifest_path.read_text())
        path = self.root / manifest["files"]["tx2gene"]["path"]
        path.write_text(path.read_text().replace("tp\tgp", "tp\twrong_gene"))
        manifest["files"]["tx2gene"]["sha256"] = sha256(path)
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValidationError, "tx2gene"):
            validate_rna_reference("derived/manifest.json", self.root)


class RNAConfigTests(unittest.TestCase):
    def setUp(self):
        self.config = yaml.safe_load((ROOT / "config/config.yaml").read_text())
        self.rna = yaml.safe_load((ROOT / "config/rnaseq.yaml").read_text())
        with (ROOT / "config/samples.tsv").open() as handle:
            self.samples = list(csv.DictReader(handle, delimiter="\t"))

    def load(self):
        with patch("src.rnaseq.config.yaml.safe_load", return_value=self.rna):
            return load_rna(self.config, self.samples)

    def test_unresolved_library_blocks(self):
        with self.assertRaisesRegex(ValidationError, "library_type/evidence"):
            self.load()

    def test_unconfirmed_de_blocks(self):
        self.rna.update(library_type="A", library_evidence="explicit inference test")
        with self.assertRaisesRegex(ValidationError, "design not confirmed"):
            self.load()

    def test_ribo_selection_blocks(self):
        self.rna["sample_ids"] = ["young_ribo_1"]
        with self.assertRaisesRegex(ValidationError, "paired-end RNA"):
            self.load()

    def test_insufficient_replicates_blocks(self):
        self.rna.update(sample_ids=["young_rna_1", "middle_rna_1"], library_type="IU", library_evidence="test", design_confirmed=True, design_evidence="test")
        with self.assertRaisesRegex(ValidationError, ">=2"):
            self.load()
