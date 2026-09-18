"""Render self-contained Markdown and HTML reports from the SQLite catalog."""
import argparse
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import sqlite3

from src.validation.validate import require, sha256


def md_table(headers, rows):
    clean = lambda value: str(value).replace("|", "\\|")
    return "\n".join(["| " + " | ".join(map(clean, headers)) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |",
                      *("| " + " | ".join(map(clean, row)) + " |" for row in rows)])


def html_table(headers, rows):
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def fmt(value, digits=4):
    return "NA" if value is None else f"{value:.{digits}f}" if isinstance(value, float) else value


def svg_bars(path, rows):
    width, height, left, right, top, bottom = 760, 380, 70, 20, 35, 90
    values = [row[1] for row in rows]
    bound = max(max(abs(value) for value in values), 1e-9)
    center = top + (height - top - bottom) / 2
    scale = (height - top - bottom) / 2 / bound
    slot = (width - left - right) / max(len(rows), 1)
    blocks = []
    for index, (gene, value) in enumerate(rows):
        x = left + index * slot + slot * 0.15
        y = center - max(value, 0) * scale
        h = abs(value) * scale
        color = "#d95f02" if value >= 0 else "#1b9e77"
        blocks.append(f'<rect x="{x:.2f}" y="{y if value >= 0 else center:.2f}" width="{slot*0.7:.2f}" height="{h:.2f}" fill="{color}"/>')
        blocks.append(f'<text x="{x+slot*0.35:.2f}" y="{height-bottom+12}" text-anchor="end" transform="rotate(-55 {x+slot*0.35:.2f},{height-bottom+12})" font-size="10">{html.escape(gene)}</text>')
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="{width/2}" y="22" text-anchor="middle" font-size="16">Largest descriptive TE changes</text>',
           f'<line x1="{left}" y1="{center:.2f}" x2="{width-right}" y2="{center:.2f}" stroke="black"/>', *blocks,
           f'<text x="{left-5}" y="{top+4}" text-anchor="end" font-size="10">+{bound:.2f}</text>',
           f'<text x="{left-5}" y="{height-bottom+4}" text-anchor="end" font-size="10">-{bound:.2f}</text>', '</svg>']
    Path(path).write_text("\n".join(svg) + "\n")


def build(task):
    database, output = Path(task["database"]), Path(task["directory"])
    output.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    try:
        runs = connection.execute("SELECT * FROM workflow_runs ORDER BY generated_at DESC").fetchall()
        require(runs, "Catalog contains no workflow run")
        run = runs[0]
        run_id = run["run_id"]
        samples = connection.execute("SELECT s.sample_id,s.assay,s.condition,s.replicate,s.run_accession FROM samples s JOIN run_samples rs USING(sample_id) WHERE rs.run_id=? ORDER BY s.assay,s.condition,s.replicate", (run_id,)).fetchall()
        qc = connection.execute("SELECT sample_id,stage,total_reads,retained_reads,mapping_rate,rrna_fraction,frame0_fraction FROM qc_metrics WHERE run_id=? ORDER BY sample_id,stage", (run_id,)).fetchall()
        summary = connection.execute("SELECT count(*) AS rows,sum(eligible) AS eligible,count(DISTINCT gene_id) AS genes FROM integration_metrics WHERE run_id=?", (run_id,)).fetchone()
        contrast_name = connection.execute("SELECT contrast FROM integration_contrasts WHERE run_id=? LIMIT 1", (run_id,)).fetchone()[0]
        top = connection.execute("SELECT gene_id,te_log2_change,rna_log2_change,ribo_log2_change FROM integration_contrasts WHERE run_id=? AND eligible=1 ORDER BY abs(te_log2_change) DESC,gene_id LIMIT 15", (run_id,)).fetchall()
        artifacts = connection.execute("SELECT kind,path,sha256,bytes FROM artifacts WHERE run_id=? ORDER BY kind,path", (run_id,)).fetchall()
        reference = connection.execute("SELECT * FROM \"references\" WHERE reference_id=?", (run["reference_id"],)).fetchone()
    finally:
        connection.close()
    config = task["config"]
    qc_rows = [[row["sample_id"], row["stage"], row["total_reads"], row["retained_reads"], fmt(row["mapping_rate"]),
                fmt(row["rrna_fraction"]), fmt(row["frame0_fraction"])] for row in qc]
    sample_rows = [[row["sample_id"], row["run_accession"], row["assay"], row["condition"], row["replicate"]] for row in samples]
    top_rows = [[row["gene_id"], fmt(row["te_log2_change"]), fmt(row["rna_log2_change"]), fmt(row["ribo_log2_change"])] for row in top]
    artifact_rows = [[row["kind"], row["path"], row["bytes"], row["sha256"]] for row in artifacts]
    group_counts = {}
    for row in samples:
        group_counts[(row["assay"], row["condition"])] = group_counts.get((row["assay"], row["condition"]), 0) + 1
    limited_replication = any(value < 2 for value in group_counts.values())
    replication_limit = ("- At least one assay/condition has only one library, so biological variance cannot be estimated."
                         if limited_replication else
                         "- This run has at least two libraries per assay/condition, but the current TE method remains descriptive and does not fit a differential-TE model.")
    title = config["report_title"] + (" — SYNTHETIC TEST" if run["synthetic"] else "")
    generated = datetime.now(timezone.utc).isoformat()
    sections = [f"# {title}", "", f"Generated: `{generated}`  ", f"Run: `{run_id}`  ",
                f"Reference: `{reference['reference_id']}` ({reference['provider']} release {reference['release']})", "",
                "## Scope and interpretation", "",
                "This report integrates RNA and Ribo libraries at the condition level. RNA and Ribo BioSamples are not treated as paired. All TE values and changes are descriptive; no differential-TE p-values or FDR are produced.", "",
                f"Normalization: RNA gene TPM is rescaled to one million per sample; Ribo P-site counts are divided by CDS length and scaled to one million per sample. Condition values are arithmetic means. TE = `log2((Ribo CDS TPM + {config['pseudocount']}) / (RNA TPM + {config['pseudocount']}))`.", "",
                f"Low-expression gates: RNA TPM ≥ {config['min_rna_tpm']}; Ribo CDS TPM ≥ {config['min_ribo_tpm']} in each condition.", "",
                "## Selected libraries", "", md_table(["Sample", "Run", "Assay", "Condition", "Replicate"], sample_rows), "",
                "## QC summary", "", md_table(["Sample", "Stage", "Input", "Retained/assigned", "Mapping", "rRNA", "Frame 0"], qc_rows), "",
                "## Integration summary", "", f"- Shared genes: {summary['genes']}", f"- Condition × gene rows: {summary['rows']}",
                f"- Rows passing both expression gates: {summary['eligible']}", f"- Contrast: `{contrast_name}` (descriptive numerator minus denominator)", "",
                "## Largest absolute descriptive TE changes", "", md_table(["Gene", "ΔTE log2", "RNA log2 change", "Ribo log2 change"], top_rows), "",
                "## Reproducibility artifacts", "", md_table(["Kind", "Path", "Bytes", "SHA-256"], artifact_rows), "",
                "## Limitations", "", "- Cross-assay biological pairing is unconfirmed; replicate numbers are not pair IDs.",
                replication_limit,
                "- TE depends on the stated normalization, CDS definition, low-expression thresholds and pseudocount.",
                "- Synthetic output validates software behavior only and has no biological interpretation.", ""]
    markdown = "\n".join(sections)
    (output / "report.md").write_text(markdown)
    plot_rows = [(row["gene_id"], row["te_log2_change"]) for row in reversed(top[:12])]
    svg_bars(output / "te_changes.svg", plot_rows or [("none", 0.0)])
    style = "body{font-family:system-ui,sans-serif;max-width:1120px;margin:2rem auto;padding:0 1rem;color:#172033}h1,h2{color:#19324d}table{border-collapse:collapse;width:100%;font-size:.86rem;margin:1rem 0 2rem}th,td{border:1px solid #ccd5df;padding:.45rem;text-align:left;vertical-align:top}th{background:#edf3f8}code{background:#eef2f5;padding:.1rem .25rem}.notice{background:#fff4d6;border-left:4px solid #d88b00;padding:1rem}img{max-width:100%}"
    html_doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>{style}</style></head><body>
<h1>{html.escape(title)}</h1><p>Generated: <code>{html.escape(generated)}</code><br>Run: <code>{html.escape(run_id)}</code><br>Reference: <code>{html.escape(reference['reference_id'])}</code> ({html.escape(reference['provider'])} release {html.escape(reference['release'])})</p>
<h2>Scope and interpretation</h2><div class="notice">Condition-level descriptive integration. RNA/Ribo BioSamples are not paired. No differential-TE p-values or FDR are produced.</div><p>RNA gene TPM and CDS-length-normalized Ribo P-sites are each scaled to one million per sample, then averaged by condition. TE uses pseudocount {config['pseudocount']}. Gates: RNA ≥ {config['min_rna_tpm']}, Ribo ≥ {config['min_ribo_tpm']}.</p>
<h2>Selected libraries</h2>{html_table(['Sample','Run','Assay','Condition','Replicate'],sample_rows)}
<h2>QC summary</h2>{html_table(['Sample','Stage','Input','Retained/assigned','Mapping','rRNA','Frame 0'],qc_rows)}
<h2>Integration summary</h2><ul><li>Shared genes: {summary['genes']}</li><li>Condition × gene rows: {summary['rows']}</li><li>Rows passing both gates: {summary['eligible']}</li><li>Contrast: <code>{html.escape(contrast_name)}</code></li></ul>
<h2>Largest absolute descriptive TE changes</h2><img src="te_changes.svg" alt="Descriptive TE change bar chart">{html_table(['Gene','ΔTE log2','RNA log2 change','Ribo log2 change'],top_rows)}
<h2>Reproducibility artifacts</h2>{html_table(['Kind','Path','Bytes','SHA-256'],artifact_rows)}
<h2>Limitations</h2><ul><li>Cross-assay pairing is unconfirmed.</li><li>{html.escape(replication_limit[2:])}</li><li>TE depends on normalization, CDS definition, filtering and pseudocount.</li><li>Synthetic results have no biological interpretation.</li></ul></body></html>"""
    (output / "report.html").write_text(html_doc)
    manifest = dict(run_id=run_id, generated_at=generated, synthetic=bool(run["synthetic"]),
                    inputs_sha256={"database": sha256(database), "catalog_provenance": sha256(Path(task["catalog_provenance"])),
                                   "integration_analysis": sha256(Path(task["integration_analysis"]))},
                    artifacts_sha256={name: sha256(output / name) for name in ("report.md", "report.html", "te_changes.svg")})
    (output / "report_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True)
    build(json.loads(parser.parse_args().task))
