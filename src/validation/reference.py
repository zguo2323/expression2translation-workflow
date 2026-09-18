"""Reference file identity contract; sequence/annotation compatibility is checked later."""
import json
import csv
from pathlib import Path
from src.validation.validate import check_schema, relative_path, require, sha256


def validate_reference(manifest_path, root=Path(".")):
    root = Path(root).resolve()
    manifest = json.loads(relative_path(root, str(manifest_path)).read_text())
    check_schema(manifest, root / "workflow/schemas/reference.schema.json", "reference")
    for role, entry in manifest["files"].items():
        path = relative_path(root, entry["path"])
        require(path.is_file(), f"Missing reference {role}: {path}")
        require(sha256(path) == entry["sha256"], f"Reference checksum mismatch: {role}")
    return manifest


def validate_rna_reference(manifest_path, root=Path(".")):
    from src.validation.sequences import reconstruct, fasta
    root = Path(root).resolve()
    manifest = validate_reference(manifest_path, root)
    files = manifest["files"]
    require("tx2gene" in files, "RNA reference requires tx2gene")
    genome, transcripts = reconstruct(root / files["genome"]["path"], root / files["annotation"]["path"])
    observed = fasta(root / files["transcriptome"]["path"])
    require(observed == {tid: t["sequence"] for tid, t in transcripts.items()}, "Transcriptome sequences/IDs disagree with genome/GTF exons")
    rrna = {tid: t["sequence"] for tid, t in transcripts.items() if t["biotype"].lower() == "rrna"}
    require(bool(rrna) and fasta(root / files["rrna"]["path"]) == rrna, "rRNA sequences disagree with genome/GTF")
    with (root / files["tx2gene"]["path"]).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["transcript_id", "gene_id"], "Invalid tx2gene header")
        rows = list(reader)
    mapping = {r["transcript_id"]: r["gene_id"] for r in rows}
    require(len(rows) == len(mapping) and mapping == {tid: t["gene_id"] for tid, t in transcripts.items()}, "tx2gene does not exactly match GTF")
    return manifest
