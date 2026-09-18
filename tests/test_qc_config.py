import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import yaml

from src.qc.config import load_qc, validate_policy
from src.validation.reference import validate_reference
from src.validation.validate import ValidationError


class QCConfigTests(unittest.TestCase):
    def setUp(self):
        self.policy = yaml.safe_load(Path("config/qc.yaml").read_text())["policies"]["riboseq"]
        self.policy.update(confirmed=True, evidence="synthetic", adapter_mode="none", barcode_mode="none", min_length=20)

    def test_confirmed_none_valid(self):
        validate_policy(self.policy, "riboseq")

    def test_unresolved_barcode_and_adapter_block(self):
        for key in ("adapter_mode", "barcode_mode", "min_length"):
            policy = copy.deepcopy(self.policy)
            policy[key] = None
            with self.subTest(key=key), self.assertRaises(ValidationError):
                validate_policy(policy, "riboseq")

    def test_conflicting_trim_and_read2_block(self):
        for key, value in (("trim_front1", 4), ("adapter_r2", "ACGT")):
            policy = copy.deepcopy(self.policy)
            policy[key] = value
            with self.subTest(key=key), self.assertRaises(ValidationError):
                validate_policy(policy, "riboseq")

    def test_raw_needs_files_but_not_trim_policy(self):
        config = yaml.safe_load(Path("config/config.yaml").read_text())
        rows = [{"sample_id": "test", "assay": "riboseq", "fastq_1": "missing.fastq.gz", "fastq_2": ""}]
        with self.assertRaisesRegex(ValidationError, "Missing FASTQ"):
            load_qc(config, rows)
        with patch("pathlib.Path.is_file", return_value=True):
            qc, samples = load_qc(config, rows)
        self.assertEqual(qc["mode"], "raw")
        self.assertEqual(list(samples), ["test"])

    def test_reference_identity_then_changed_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shutil.copytree("workflow/schemas", root / "workflow/schemas")
            files = {}
            for role in ("genome", "annotation", "transcriptome", "rrna"):
                (root / role).write_text("synthetic test\n")
                files[role] = dict(path=role, sha256=hashlib.sha256((root / role).read_bytes()).hexdigest(),
                                   source="synthetic", retrieved_at="2026-09-16", derivation=None)
            manifest = dict(reference_id="synthetic", provider="test", release="test-v1", files=files)
            (root / "reference.json").write_text(json.dumps(manifest))
            self.assertEqual(validate_reference("reference.json", root), manifest)
            (root / "genome").write_text("changed")
            with self.assertRaisesRegex(ValidationError, "checksum"):
                validate_reference("reference.json", root)
