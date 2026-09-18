"""Small argument-list adapters; commands and versions are written alongside artifacts."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import tempfile

from src.qc.config import validate_policy


def run_task(task):
    output = Path(task["directory"])
    output.mkdir(parents=True, exist_ok=True)
    tool = task["tool"]
    version = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True)
    command = [tool]
    if tool == "fastqc":
        # Unique per-sample/read output name avoids sample collisions in MultiQC.
        with tempfile.TemporaryDirectory(prefix="input-", dir=output) as staging:
            link = Path(staging) / f"{task['label']}.fastq.gz"
            link.symlink_to(Path(task["reads"][0]).resolve())
            command += ["--threads", "1", "--outdir", str(output), "--noextract", str(link)]
            subprocess.run(command, check=True)
    elif tool == "fastp":
        policy = task["policy"]
        validate_policy(policy, task["assay"])
        command += ["--in1", task["reads"][0], "--out1", task["outputs"][0],
                    "--thread", str(task["threads"]), "--json", str(output / "fastp.json"),
                    "--html", str(output / "fastp.html"), "--report_title", task["label"],
                    "--disable_trim_poly_g", "--dont_eval_duplication",
                    "--length_required", str(policy["min_length"]),
                    "--qualified_quality_phred", str(policy["quality_phred"]),
                    "--unqualified_percent_limit", str(policy["unqualified_percent"]),
                    "--trim_front1", str(policy["trim_front1"])]
        if len(task["reads"]) == 2:
            command += ["--in2", task["reads"][1], "--out2", task["outputs"][1], "--trim_front2", str(policy["trim_front2"])]
        if policy["adapter_mode"] == "none":
            command += ["--disable_adapter_trimming"]
        else:
            command += ["--adapter_sequence", policy["adapter_r1"]]
            if len(task["reads"]) == 2:
                command += ["--adapter_sequence_r2", policy["adapter_r2"]]
        subprocess.run(command, check=True)
    elif tool == "multiqc":
        # Explicit input list prevents stale reports from other sample subsets/modes leaking in.
        command += [*task["inputs"], "--outdir", str(output), "--filename", "multiqc_report.html",
                    "--force", "--title", ("SYNTHETIC TEST — " if task["synthetic"] else "") + task["label"], "--dirs"]
        subprocess.run(command, check=True)
    else:
        raise ValueError(f"Unsupported tool: {tool}")
    (output / f"{task['label']}.provenance.json").write_text(json.dumps({
        "tool": tool, "version": (version.stdout + version.stderr).strip(), "command": command,
        "parameters": task, "timestamp": datetime.now(timezone.utc).isoformat(),
        "synthetic": task["synthetic"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    run_task(json.loads(parser.parse_args().task))
