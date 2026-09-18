"""Exercise the batch with fake tools, including real process interruption."""
import os
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        (self.root / "runs.txt").write_text("SRR1\nSRR2")  # deliberately no final newline
        (self.root / "samples.tsv").write_text("run_accession\nSRR1\nSRR2\n")
        self.command = [sys.executable, str(ROOT / "scripts/download_reads.py"), "--accessions", "runs.txt", "--samples", "samples.tsv"]
        self.env = {**os.environ, "PATH": str(self.bin) + os.pathsep + os.environ["PATH"]}

    def fake(self, body):
        tool = self.bin / "prefetch"
        tool.write_text(f"#!{sys.executable}\n" + body)
        tool.chmod(0o755)

    def run_batch(self):
        return subprocess.run(self.command, cwd=self.root, env=self.env, capture_output=True, text=True)

    def test_last_line_and_existing_nested_resume_path(self):
        nested = self.root / "data/sra/SRR1/SRR1"
        nested.mkdir(parents=True)
        (nested / "SRR1.sra.tmp").write_text("partial")
        self.fake("import sys\nfrom pathlib import Path\nr=sys.argv[1]\np=Path(sys.argv[-1])/r\np.mkdir(parents=True,exist_ok=True)\n(p/(r+'.sra')).write_text('complete')\n")
        result = self.run_batch()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((nested / "SRR1.sra").exists())
        self.assertTrue((self.root / "data/sra/SRR2/SRR2.sra").exists())
        self.assertEqual((nested / "SRR1.sra.tmp").read_text(), "partial")

    def test_failure_stops_batch(self):
        self.fake("import sys\nsys.exit(7)\n")
        result = self.run_batch()
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("Processing SRR2", result.stdout)

    def test_interrupt_stops_child_and_batch(self):
        self.fake("import time,os\nfrom pathlib import Path\nPath('child.pid').write_text(str(os.getpid()))\ntime.sleep(30)\n")
        child = subprocess.Popen(self.command, cwd=self.root, env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            deadline = time.monotonic() + 5
            while not (self.root / "child.pid").exists() and time.monotonic() < deadline:
                time.sleep(0.02)
            self.assertTrue((self.root / "child.pid").exists())
            child.send_signal(signal.SIGINT)
            out, err = child.communicate(timeout=8)
            self.assertEqual(child.returncode, 130, err)
            self.assertNotIn("Processing SRR2", out)
            pid = int((self.root / "child.pid").read_text())
            with self.assertRaises(ProcessLookupError):
                os.kill(pid, 0)
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()

    def test_existing_lock_not_removed(self):
        directory = self.root / "data/sra/SRR1"
        directory.mkdir(parents=True)
        lock = directory / "SRR1.sra.lock"
        lock.touch()
        self.fake("raise AssertionError('must not run')\n")
        result = self.run_batch()
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(lock.exists())

    def test_convert_skip_verified_and_refuse_changed_outputs(self):
        (self.root / "runs.txt").write_text("SRR1")
        (self.root / "samples.tsv").write_text("run_accession\tlayout\tfastq_1\tfastq_2\nSRR1\tSINGLE\tdata/raw/SRR1.fastq.gz\t\n")
        self.fake("import sys\nfrom pathlib import Path\np=Path(sys.argv[-1])/'SRR1'\np.mkdir(parents=True,exist_ok=True)\n(p/'SRR1.sra').write_text('complete')\n")
        validate = self.bin / "vdb-validate"
        validate.write_text(f"#!{sys.executable}\n")
        validate.chmod(0o755)
        dump = self.bin / "fasterq-dump"
        dump.write_text(f"#!{sys.executable}\nimport sys\nfrom pathlib import Path\np=Path(sys.argv[sys.argv.index('--outdir')+1])/'SRR1.fastq'\np.write_text('@one\\nACGT\\n+\\nIIII\\n')\n")
        dump.chmod(0o755)
        self.command += ["--convert"]
        first = self.run_batch()
        self.assertEqual(first.returncode, 0, first.stderr)
        target = self.root / "data/raw/SRR1.fastq.gz"
        original = target.read_bytes()
        event = json.loads((self.root / "logs/download/SRR1.status.json").read_text())
        self.assertEqual(event["status"], "converted")
        self.assertEqual(event["reads"]["files"][0]["path"], "data/raw/SRR1.fastq.gz")
        second = self.run_batch()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("Skipping verified", second.stdout)
        self.assertEqual(target.read_bytes(), original)
        target.write_bytes(b"corrupted")
        third = self.run_batch()
        self.assertNotEqual(third.returncode, 0)
        self.assertEqual(target.read_bytes(), b"corrupted")
