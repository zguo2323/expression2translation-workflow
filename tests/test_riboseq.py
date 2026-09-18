import copy
import csv
from contextlib import chdir
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import yaml

from src.riboseq.config import load_offsets, load_ribo
from src.riboseq.features import derive_cds_features, read_features, write_features
from src.riboseq.run import psite_coordinate
from src.validation.validate import ValidationError


class RiboFeatureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_positive_and_negative_cds_to_transcript_coordinates(self):
        gtf = self.root / "annotation.gtf"
        gtf.write_text(
            'chr\tx\texon\t100\t109\t.\t+\t.\tgene_id "gp"; transcript_id "tp";\n'
            'chr\tx\texon\t200\t209\t.\t+\t.\tgene_id "gp"; transcript_id "tp";\n'
            'chr\tx\tCDS\t104\t109\t.\t+\t0\tgene_id "gp"; transcript_id "tp";\n'
            'chr\tx\tCDS\t200\t205\t.\t+\t0\tgene_id "gp"; transcript_id "tp";\n'
            'chr\tx\texon\t300\t309\t.\t-\t.\tgene_id "gm"; transcript_id "tm";\n'
            'chr\tx\texon\t400\t409\t.\t-\t.\tgene_id "gm"; transcript_id "tm";\n'
            'chr\tx\tCDS\t302\t309\t.\t-\t0\tgene_id "gm"; transcript_id "tm";\n'
            'chr\tx\tCDS\t400\t403\t.\t-\t0\tgene_id "gm"; transcript_id "tm";\n')
        rows = derive_cds_features(gtf, {"tp": 20, "tm": 20})
        by_id = {row["transcript_id"]: row for row in rows}
        self.assertEqual((by_id["tp"]["cds_start"], by_id["tp"]["cds_end"]), (4, 16))
        self.assertEqual((by_id["tm"]["cds_start"], by_id["tm"]["cds_end"]), (6, 18))
        path = self.root / "features.tsv"
        write_features(path, rows)
        self.assertEqual(read_features(path)["tm"]["gene_id"], "gm")

    def test_non_contiguous_cds_is_rejected(self):
        gtf = self.root / "bad.gtf"
        gtf.write_text('chr\tx\texon\t1\t20\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
                       'chr\tx\tCDS\t2\t4\t.\t+\t0\tgene_id "g"; transcript_id "t";\n'
                       'chr\tx\tCDS\t7\t9\t.\t+\t0\tgene_id "g"; transcript_id "t";\n')
        with self.assertRaisesRegex(ValidationError, "not contiguous"):
            derive_cds_features(gtf, {"t": 20})

    def test_psite_coordinates_both_alignment_directions(self):
        self.assertEqual(psite_coordinate(100, 128, False, 28, 12), 112)
        self.assertEqual(psite_coordinate(100, 128, True, 28, 12), 115)
        with self.assertRaisesRegex(ValidationError, "outside"):
            psite_coordinate(100, 128, False, 28, 28)


class RiboConfigTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "offsets.tsv").write_text("read_length\tpsite_offset\n28\t12\n")
        self.ribo = yaml.safe_load(Path("config/riboseq.yaml").read_text())
        self.ribo.update(strand_mode="forward", strand_evidence="synthetic protocol", offset_table="offsets.tsv",
                         offset_evidence="synthetic planted P-sites")
        (self.root / "ribo.yaml").write_text(yaml.safe_dump(self.ribo, sort_keys=False))
        self.config = {"ribo_config": "ribo.yaml"}
        self.samples = [dict(sample_id=sid, assay="riboseq", layout="SINGLE") for sid in self.ribo["sample_ids"]]

    def load(self):
        with chdir(self.root), patch("src.riboseq.config.validate_rna_reference", return_value={"reference_id": "test", "files": {}}):
            return load_ribo(self.config, self.samples)

    def test_valid_ribo_config_and_offsets(self):
        ribo, selected, _, offsets = self.load()
        self.assertEqual(set(selected), set(ribo["sample_ids"]))
        self.assertEqual(offsets, {28: 12})

    def test_unresolved_method_fields_block(self):
        for key in ("strand_mode", "strand_evidence", "offset_table", "offset_evidence"):
            changed = copy.deepcopy(self.ribo)
            changed[key] = None
            (self.root / "ribo.yaml").write_text(yaml.safe_dump(changed, sort_keys=False))
            with self.subTest(key=key), self.assertRaises(ValidationError):
                self.load()

    def test_offsets_require_full_valid_unique_rows(self):
        for body in ("read_length\tpsite_offset\n28\t28\n", "read_length\tpsite_offset\n28\t12\n28\t13\n", "x\ty\n28\t12\n"):
            path = self.root / "bad.tsv"
            path.write_text(body)
            with self.subTest(body=body), self.assertRaises((ValidationError, ValueError)):
                load_offsets(path)


if __name__ == "__main__":
    unittest.main()
