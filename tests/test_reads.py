import gzip
from pathlib import Path
import tempfile
import unittest
from src.validation.reads import validate_reads
from src.validation.validate import ValidationError
from src.qc.config import validate_policy


class ReadsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fastq(self, name, text):
        path = self.root / name
        with gzip.open(path, "wt") as f:
            f.write(text)
        return path

    def test_pair_and_counts(self):
        a = self.fastq("a.gz", "@one/1\nACGT\n+\nIIII\n")
        b = self.fastq("b.gz", "@one/2\nTGCA\n+\nIIII\n")
        result = validate_reads([a, b])
        self.assertEqual(result["files"][0]["reads"], 1)
        self.assertEqual(result["files"][1]["bases"], 4)

    def test_mismatch_ids(self):
        a = self.fastq("a.gz", "@one\nAC\n+\nII\n")
        b = self.fastq("b.gz", "@two\nAC\n+\nII\n")
        with self.assertRaisesRegex(ValidationError, "IDs/order"):
            validate_reads([a, b])

    def test_mismatch_count(self):
        a = self.fastq("a.gz", "@one\nAC\n+\nII\n")
        b = self.fastq("b.gz", "@one\nAC\n+\nII\n@two\nAC\n+\nII\n")
        with self.assertRaisesRegex(ValidationError, "counts differ"):
            validate_reads([a, b])

    def test_invalid_record_and_empty(self):
        for content in ("", "@one\nAC\n+\nI\n", "@one\nAC\n+\n", "x\nAC\n+\nII\n"):
            with self.subTest(content=content), self.assertRaises(ValidationError):
                validate_reads([self.fastq("bad.gz", content)])

    def test_truncated_gzip(self):
        a = self.fastq("a.gz", "@one\nAC\n+\nII\n")
        a.write_bytes(a.read_bytes()[:-5])
        with self.assertRaisesRegex(ValidationError, "gzip"):
            validate_reads([a])

    def test_unconfirmed_policy_blocks(self):
        import yaml
        policy = yaml.safe_load(Path("config/qc.yaml").read_text())["policies"]["riboseq"]
        with self.assertRaisesRegex(ValidationError, "confirmed"):
            validate_policy(policy, "riboseq")
