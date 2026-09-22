# 描述性 TE、SQLite 与报告

阶段 4 将已经通过各自门禁的 RNA 和 Ribo 结果按统一 gene ID 联结，生成 condition-level 描述性 translation efficiency（TE）、SQLite catalog、跨分支 provenance 和 Markdown/HTML 报告。

## 统计契约

RNA 与 Ribo 的 BioSample ID 不同，工作流不把相同 replicate 编号解释为配对关系。处理顺序为：

1. 先取 RNA gene 与 Ribo CDS 的共同编码基因集合。RNA-only 的 rRNA 和其他 biotype 不进入 TE 归一化分母。
2. RNA 使用 tximport 汇总的 gene TPM，并在共同编码基因集合内将每个样本重新缩放至总和 1,000,000。
3. Ribo 使用共同编码基因的 CDS 内 P-site counts，先除以 CDS kb，再在每个样本内缩放至总和 1,000,000；输出名为 `ribo_cds_tpm`，表示 TPM-like density。
4. 每个 assay 分别按 condition 对样本取算术均值。
5. 每个 condition/gene 必须同时达到 `min_rna_tpm` 和 `min_ribo_tpm`；未达到的行保留过滤原因，TE 留空。
6. 合格行计算 `log2((ribo_cds_tpm_mean + pseudocount) / (rna_tpm_mean + pseudocount))`。
7. `te_contrast.tsv` 报告配置 numerator−denominator 的描述性 RNA、Ribo 和 TE log2 变化，不计算 p 值或 FDR。

为防止低 P-site 计数因伪计数和归一化而占据排序前列，工作流还会生成候选证据表。默认要求每个 Ribo library 至少 10 个原始 CDS P-sites，并检查 TE 方向在 `0.1`、`1`、`10` 三个 pseudocount 下，以及逐一移除每个 RNA 或 Ribo library 后是否保持一致。只有同时通过这些检查的行才标记为 `candidate_stable`；它仍然是描述性候选，不是差异 TE 显著性结果。

`differential_te_readiness.json` 单独审查是否可运行基于 replicate-level raw counts 的 assay×condition 负二项交互模型。审查会记录每组重复数、计数支持的基因数、各 Ribo library 的 CDS P-site 深度和 rRNA fraction。未达到门槛时 `status: not_run`，而非产生不可靠的 p 值。

所有单位、聚合方式、阈值、pseudocount 和 contrast 固定在 `config/integration.yaml` 并写入 `analysis.json`。MVP 明确禁止 `differential_te: true`。

这种计算适合展示同一数据集中的相对变化，不等同于带 assay×condition interaction 的 differential-TE 统计模型。本次默认完整 8-run，每个 condition × assay 两个 libraries；重复数量仍有限。4-run 仅为可选描述性示例。即使 8-run 到位，也不为本模型增加 differential TE p 值/FDR。

RNA/integration 目前只支持一个两组 contrast，两组均需具备所选 RNA 与 Ribo；活跃表可以登记更多条件，但必须显式选择该 contrast 子集。输入契约见 [data-input.md](data-input.md)。

## 运行

完整 8-run 配置在确认 RNA/Ribo 方法参数后运行：

```bash
export PATH="$PWD/.conda/riboseq/bin:$PWD/.conda/qc/bin:$PWD/.conda/salmon/bin:$PWD/.conda/rnaseq-stats/bin:$PATH"
.venv/bin/snakemake --cores 4 --config stage=integration --dry-run
.venv/bin/snakemake --cores 4 --config stage=integration
```

4-run 最小真实演示使用独立配置：

```bash
.venv/bin/snakemake --cores 4 \
  --config stage=integration \
           rna_config=config/rnaseq-minimal.yaml \
           ribo_config=config/riboseq-minimal.yaml \
  --dry-run
```

首次运行自有数据前仍需在最小 RNA/Ribo 配置和 `config/qc.yaml` 填入由该数据方法与 reads 支持的 adapter、library type、strand、显式 footprint lengths 和 P-site offsets。仓库默认值属于 GSE203147，不能复制给其他数据。RNA 最小配置固定 `run_deseq2: false`。

## 输出

| 路径 | 内容 |
| --- | --- |
| `results/integration/rna_sample_tpm.tsv` | RNA 每样本重缩放 gene TPM |
| `results/integration/ribo_sample_cds_tpm.tsv` | Ribo 每样本 CDS-length-normalized TPM-like density |
| `results/integration/condition_te.tsv` | condition/gene 均值、样本数、过滤原因与描述性 TE |
| `results/integration/te_contrast.tsv` | numerator−denominator 的描述性 RNA/Ribo/TE 变化 |
| `results/integration/te_candidate_qc.tsv` | 每基因的最小原始 P-site、pseudocount 与 leave-one-out 方向稳定性 |
| `results/integration/differential_te_readiness.json` | 差异 TE 交互模型的运行就绪审查与阻断原因 |
| `results/integration/analysis.json` | 有效方法、样本、gene 数、输入/产物 SHA-256 |
| `results/catalog.db` | 结构化 catalog，不包含 FASTQ/BAM |
| `results/provenance/catalog.json` | 确定性 run ID、配置与所有输入摘要、表行数、SQLite integrity |
| `results/report/report.md` | 可审阅的文本报告 |
| `results/report/report.html` | 自包含样式的 HTML 报告，引用本地 SVG |
| `results/report/te_changes.svg` | 最大绝对描述性 TE 变化图 |
| `results/report/report_manifest.json` | 报告输入与产物 SHA-256 |

## SQLite catalog

核心表：

- `references`、`genes`：固定 reference 和 GTF gene 坐标。
- `samples`、`workflow_runs`、`run_samples`：样本事实和确定性运行身份。
- `qc_metrics`：Salmon mapping 与 Ribo contamination/mapping/periodicity。
- `rna_gene_metrics`、`ribo_gene_metrics`：每样本定量。
- `integration_metrics`、`integration_contrasts`：condition TE 和描述性 contrast。
- `artifacts`：重要结果路径、bytes、SHA-256 和类型。

schema v2 中来源字段允许 SQL NULL，replicate 为文本。旧 schema 不自动覆盖或迁移：dry-run/写入前拒绝旧库，保留旧文件并设置新的 `paths.results` 目录重建。

数据库启用 foreign keys、CHECK/UNIQUE 约束和事务；相同输入/参数得到相同 run ID，重复写入使用主键更新，不复制结果。每次生成后执行 `foreign_key_check` 和 `quick_check`。

示例查询：

```sql
SELECT gene_id, te_log2_change, rna_log2_change, ribo_log2_change
FROM integration_contrasts
WHERE run_id = :run_id AND eligible = 1
ORDER BY abs(te_log2_change) DESC
LIMIT 20;
```

## 合成端到端验收

```bash
.venv/bin/python scripts/smoke_integration.py --directory results/new-synthetic-e2e-smoke
```

脚本独立生成 local 样本表、统一 reference、4 个 RNA PE 和 4 个 Ribo SE libraries，不使用真实 GEO/SRA 身份，运行 QC、Salmon、tximport、rRNA removal、Bowtie2/SAMtools、P-site、TE、SQLite 和报告。5 个基因预置 Middle Ribo 增强，最终必须占据 ΔTE 前五；同时检查 SQLite 外键、无修改重跑和 pseudocount 只使阶段 4 失效。synthetic 结果不具有生物学含义。
