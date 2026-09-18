"""Record SHA-256 for gene-analysis inputs and all generated artifacts."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from src.validation.validate import sha256


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    task = json.loads(parser.parse_args().task)
    out = Path(task["directory"])
    inputs = [s["quant"] for s in task["samples"]] + [task["tx2gene"], task["config"]["reference_manifest"], "src/rnaseq/gene_analysis.R"]
    provenance = dict(generated_at=datetime.now(timezone.utc).isoformat(), parameters=task,
        inputs_sha256={p: sha256(Path(p)) for p in inputs},
        artifacts_sha256={str(p): sha256(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "provenance.json"})
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
