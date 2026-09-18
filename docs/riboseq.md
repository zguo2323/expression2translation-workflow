# Ribo-seq 分支与 P-site 契约

当前实现：已确认的 Ribo single-end 预处理 → rRNA 去除 → 转录本比对 → P-site/CDS gene counts → read-length、frame periodicity、起始位点 metagene、汇总与 provenance。
它产出 assay 内可复核的 Ribo 定量，不把 RNA/Ribo 重复编号当成生物学配对，也不提前计算 TE。

## 方法与坐标定义

Ribo 与 RNA 使用同一个 Ensembl release 115 / R64-1-1 manifest。Bowtie2 先对 manifest 中的 rRNA FASTA 比对并保留未比对 reads，再对同一 manifest 的 transcriptome 比对；SAMtools 生成排序 BAM 和索引。

CDS 由 manifest 对应 GTF 动态转换为转录本 0-based、右开坐标。转换显式处理负基因组链和多 exon transcript，并要求每个 coding gene 只有一个 transcript；当前真实注释检查得到 6,600 个 coding transcripts/genes。这个限制适合本项目的酵母 reference，若以后换成多 isoform 注释，需要先定义 canonical transcript 或 gene-level ambiguity 策略。

P-site offset 表格式：

```text
read_length	psite_offset
28	12
29	12
```

`psite_offset` 表示从 read 5' 端到 P-site 的 0-based 距离。正向 alignment 使用 `reference_start + offset`，反向 alignment 使用 `reference_end - 1 - offset`。通过 MAPQ、唯一性、链方向和 CIGAR 过滤后的每一种 read length 都必须有 offset；缺失长度会停止计数。带 indel 或 clipping 的 alignment 不参与 P-site 计数，并在过滤指标中单列。

`strand_mode` 可为 `forward`、`reverse` 或有明确证据的 `unstranded`。转录本参考全部以 5'→3' 方向保存，因此 forward 表示 read 按 transcript orientation 比对。

## 真实数据启动门槛

在运行 `stage=riboseq` 前需要完成：

1. 4 个 Ribo FASTQ 完整下载并通过输入校验。
2. `config/qc.yaml` 设为 `mode: trim`，Ribo adapter、barcode/前端裁剪和长度阈值有明确证据。
3. `config/riboseq.yaml` 中的 `strand_mode`、`strand_evidence`、`offset_table`、`offset_evidence` 已填写。
4. offset 来自本数据的 read-length、起始位点 metagene和 3-nt periodicity 评估；不能把 synthetic 的 28 nt/12 nt 参数复制到真实配置。
5. 根据真实数据复核 rRNA fraction、mapping rate 和 frame-0 fraction 门槛。当前默认值是工程门禁，不是通用生物学标准。

默认配置保留这些方法字段为 `null`，因此误把阶段切换到 `riboseq` 会在 DAG 构建时失败，不会用猜测参数继续运行。

## 环境与运行

本机实测版本为 Bowtie2 2.5.4、SAMtools 1.21、pysam 0.22.1。隔离环境已创建于 `.conda/riboseq`；新机器可运行：

```bash
conda env create --prefix .conda/riboseq --file envs/riboseq.yaml
```

本机 Conda 23.9 低于 Snakemake 自动环境部署要求，所以实际测试用显式 PATH：

```bash
export PATH="$PWD/.conda/riboseq/bin:$PWD/.conda/qc/bin:$PATH"
.venv/bin/snakemake --cores 4 --config stage=riboseq --dry-run
.venv/bin/snakemake --cores 4 --config stage=riboseq
```

具备新版 Conda/mamba 时可从 workflow 环境使用 `--use-conda`，让每条 rule 按 `envs/qc.yaml` 和 `envs/riboseq.yaml` 自动部署。

## 输出

`results/riboseq/` 包含：

| 路径 | 内容 |
| --- | --- |
| `reference/` | rRNA/transcriptome Bowtie2 index、CDS transcript feature table及校验摘要 |
| `decontaminated/{sample}/` | 保留 FASTQ、输入/污染/保留 reads、rRNA fraction、版本和 SHA-256 |
| `alignment/{sample}/` | 排序 BAM/BAI、flagstat、mapping rate、版本和 SHA-256 |
| `sites/{sample}/gene_counts.tsv` | 每 gene 的 P-site count、CDS length、每 kb density 和 RPM |
| `sites/{sample}/read_lengths.tsv` | 各 read length 的 offset、比对数、可计数数和 CDS P-sites |
| `sites/{sample}/periodicity.*` | frame 0/1/2 counts、fraction 和 SVG |
| `sites/{sample}/start_metagene.*` | CDS start 周围的 P-site 表和 SVG |
| `summary/qc_metrics.tsv` | 4 个样本的 rRNA、mapping、assigned sites 和 frame-0 指标 |
| `summary/gene_counts.tsv` | condition/sample/gene 长表，供阶段 4 联合分析使用 |
| `summary/provenance.json` | 有效配置、reference/config、样本输入和汇总产物 SHA-256 |

gene counts 只统计 CDS 内 P-sites。RPM 和 `p_sites_per_kb` 是描述性输出；后续 TE 不直接用未经 assay 内归一化的原始 RNA/Ribo counts 相除。

## 合成验收

```bash
.venv/bin/python scripts/smoke_riboseq.py --directory results/new-synthetic-ribo-smoke
```

脚本生成 50 个 coding transcripts、1 个 rRNA transcript 和 4 个 synthetic Ribo libraries。P-sites 按已知 28 nt read length、12 nt offset植入 frame 0，并混入已知 rRNA reads；它检查真实 QC/Bowtie2/SAMtools 路径、4 个样本汇总、周期性 SVG、无修改重跑以及 `min_mapq` 变化后的重新调度。
这些数据只验证工程和坐标契约，不代表 GSE203147 的真实 read-length、offset、污染率或生物学结果。
