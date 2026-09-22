# Expression2Translation Workflow (E2T)

E2T 是一个由 Snakemake 驱动、可复现的 *Saccharomyces cerevisiae* bulk RNA-seq + Ribo-seq 联合分析 MVP。它从 raw reads 出发，生成 RNA/Ribo 质量控制、gene-level 定量、condition-level 描述性 translation efficiency（TE）、SQLite catalog、provenance 以及 Markdown/HTML 报告。

项目针对 GEO [GSE203147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203147) / BioProject PRJNA838843 / SRA study SRP375616 的 Young 与 Middle 样本。仓库只保存代码、小型来源元数据和配置；FASTQ、SRA、reference、BAM、日志和分析结果均由 Git 忽略。

## 当前状态

完整 MVP 工程链路已经实现，并通过真实分析工具参与的合成端到端测试：

```text
metadata validation
  → FASTQ validation → FastQC / fastp / MultiQC
  ├─ RNA: Salmon → tximport → optional DESeq2
  └─ Ribo: rRNA removal → Bowtie2 / SAMtools → P-site and CDS counts
  → condition-level descriptive TE
  → SQLite catalog + provenance + Markdown/HTML report
```

默认 `stage: metadata` 只校验元数据，不会下载 reads 或启动分析。本研究实例的 adapter、RNA library type、Ribo strand、显式 footprint 长度和 P-site offset 已由真实 reads 诊断写入配置；新数据集仍必须用自己的方法与 reads 证据替换这些项目参数，不能照搬 GSE203147 的值。

同类酵母项目可直接输入自有 FASTQ，无需 accession；样本数不固定。详见 [新用户输入指南](docs/data-input.md)和 [local 配置示例](config/local.example.yaml)。

## 数据设计

- 本次默认活跃集合为 8 个 libraries：4 个 paired-end RNA-seq 和 4 个 single-end Ribo-seq。
- Young/Middle 每个 condition、每种 assay 各有 2 个重复。
- SRA 将 Ribo assay 标为 `OTHER`，因此 [samples.tsv](config/samples.tsv) 人工记录 `assay=riboseq`。
- RNA 与 Ribo 的 BioSample ID 不同；E2T 不把相同 replicate 编号视为严格配对。
- Reference 固定为 Ensembl release 115 / R64-1-1，所有分支使用同一套 gene/transcript ID。
- 当前联合分析只报告 condition-level 描述性 TE，不报告 differential TE 的 p 值或 FDR。

## 安装

在仓库根目录创建 workflow 环境：

```bash
conda env create -f envs/workflow.yaml
conda activate e2t-workflow
```

也可以只为 Snakemake 和校验代码创建 Python 3.11 虚拟环境：

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r envs/requirements.txt
```

QC、Salmon、R/DESeq2 和 Ribo 工具分别固定在 `envs/qc.yaml`、`envs/salmon.yaml`、`envs/rnaseq-stats.yaml` 和 `envs/riboseq.yaml`。支持的 Conda/mamba 环境可由 Snakemake 使用 `--use-conda` 自动创建。

## 第一次运行

先执行不需要 FASTQ 或 reference 的 metadata 校验：

```bash
python -m src.validation.validate
snakemake --snakefile workflow/Snakefile --cores 1 --dry-run
snakemake --snakefile workflow/Snakefile --cores 1
```

成功后会生成 `results/validation/metadata.json` 和 `logs/validate_metadata.log`。这些文件只说明样本表与来源元数据一致，不表示 reads 或分析参数已经验收。

运行测试：

```bash
python -m unittest discover -s tests -v
```

## 实验室项目 registry（第一阶段）

安装本仓库后可使用 `e2t` 命令创建并登记多个独立项目。registry 记录项目状态、历史和最终产物索引；原始 FASTQ、BAM 与分析产物仍保存在项目目录中，不写入 SQLite。

```bash
python -m pip install -e .
# 在无网络且环境缺少 wheel 的开发机上：python setup.py develop
e2t init aging-yeast-2026 --root /srv/lab-e2t --title "Yeast aging" --owner alice
e2t set-status aging-yeast-2026 --root /srv/lab-e2t --stage qc --status running
e2t add-artifact aging-yeast-2026 --root /srv/lab-e2t \
  --path /srv/lab-e2t/projects/aging-yeast-2026/results/report/report.html --kind html_report
e2t status aging-yeast-2026 --root /srv/lab-e2t --json
```

目前 CLI 只管理项目目录和 registry，不会代替 Snakemake 调度分析任务。项目级 `validate/run/report`、容器和集群 profile 属于后续阶段。目录约定、数据边界和命令说明见 [实验室 CLI 指南](docs/lab-cli.md)。

## 获取 reads

下载脚本默认从活跃 samplesheet 提取有效 SRR/ERR/DRR，也可用 `--accessions` 指定子集，调用 SRA Toolkit 的 `prefetch`、`vdb-validate` 和 `fasterq-dump`，然后压缩并验证 FASTQ。先确认 SRA Toolkit 已在 `PATH`：

```bash
prefetch --version
fasterq-dump --version
```

查看计划，不发起网络请求：

```bash
python scripts/download_reads.py --dry-run
```

首次准备且文件尚不存在时，下载并转换活跃 runs（当前默认 8 个）：

```bash
python scripts/download_reads.py --convert --threads 4
```

本次 8-run 已由用户确认下载完成，已有文件应直接校验。以下 4-run 仅作为首次准备其他实例时的可选子集示例：

```bash
python scripts/download_reads.py \
  --accessions config/minimal-demo-accessions.txt \
  --convert --threads 4
```

脚本支持续传，按 run 写入状态、日志、reads 数量和 SHA-256。FASTQ 保存到 `data/raw/`，SRA archive 保存到 `data/sra/`；两者都不会被 Git 跟踪。详细说明见 [数据下载文档](docs/data-download.md)。

## 运行分析

真实数据分析前必须完成以下配置：

1. 在 [QC 配置](config/qc.yaml) 中确认 adapter、barcode/前端裁剪、最短长度和证据来源，并将 `mode` 设为 `trim`。
2. 在 RNA 配置中确认 Salmon library type；启用 DESeq2 时还需确认实验设计。
3. 在 Ribo 配置中确认 strand mode，并提供从真实 read-length、periodicity 和 metagene 证据得到的 P-site offset 表。
4. 下载并按 [RNA/reference 文档](docs/rnaseq.md) 构建、验证固定 reference。

每次实际运行前先 dry-run。例如完整 8-run 联合分析：

```bash
# 本机 Conda 23.9 不满足 Snakemake 自动部署的最低版本；使用已创建的隔离环境。
export PATH="$PWD/.conda/riboseq/bin:$PWD/.conda/qc/bin:$PWD/.conda/salmon/bin:$PWD/.conda/rnaseq-stats/bin:$PATH"
.venv/bin/snakemake --snakefile workflow/Snakefile \
  --cores 4 --config stage=integration --dry-run

.venv/bin/snakemake --snakefile workflow/Snakefile \
  --cores 4 --config stage=integration
```

在 Conda >=24.7.1 的机器上，可改用 `--use-conda` 让 Snakemake 按 `envs/*.yaml` 创建规则环境。

4-run 演示集合使用独立 RNA/Ribo 配置，只生成描述性结果：

```bash
snakemake --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  --config stage=integration \
           rna_config=config/rnaseq-minimal.yaml \
           ribo_config=config/riboseq-minimal.yaml \
  --dry-run
```

可用阶段为 `metadata`、`qc`、`rnaseq`、`riboseq` 和 `integration`。各分支的参数、输出和统计边界见 [QC](docs/qc.md)、[RNA-seq](docs/rnaseq.md)、[Ribo-seq](docs/riboseq.md) 与 [联合分析](docs/integration.md) 文档。

## 主要输出

| 路径 | 内容 |
| --- | --- |
| `results/qc/` | FASTQ 校验、FastQC、fastp、MultiQC 和工具 provenance |
| `results/rnaseq/` | Salmon 定量、gene counts/TPM、可选 DESeq2 与样本 QC |
| `results/riboseq/` | 去污染、比对、P-site/CDS counts、periodicity 与 metagene |
| `results/integration/` | condition-level RNA/Ribo abundance、描述性 TE 与 contrast |
| `results/catalog.db` | samples、reference、QC、定量、TE 和 artifacts catalog |
| `results/report/` | Markdown/HTML 报告、SVG 图和文件摘要 |
| `results/provenance/` | 配置、输入、软件和产物的可追溯记录 |

SQLite 只保存结构化结果和文件索引，不保存 FASTQ 或 BAM。

旧 SQLite catalog 需保留并在新的 `paths.results` 目录重建为 schema v2，详见 [输入迁移说明](docs/data-input.md)。

## 合成端到端验收

下面的命令生成小型合成 reference 和 reads，并执行完整工具链；不会使用真实 GEO 数据：

```bash
python scripts/smoke_integration.py --use-conda \
  --directory results/synthetic-e2e-smoke
```

GitHub Actions 也会运行 unit/DAG tests 及 QC、RNA、Ribo 和完整联合分析 smoke tests。合成结果只验证工程、坐标和统计契约，不具有生物学含义。

## 仓库结构

| 路径 | 内容 |
| --- | --- |
| `workflow/` | Snakemake 入口、rules 和 JSON schemas |
| `src/` | metadata/QC/RNA/Ribo/integration/catalog/report 实现 |
| `config/` | 样本表、阶段配置和最小演示配置 |
| `metadata/` | GEO/SRA 来源快照、manifest 和 reference manifest |
| `envs/` | 可复现的软件环境定义 |
| `scripts/` | 下载、reference 构建和 synthetic smoke 入口 |
| `tests/` | 单元测试、异常输入测试和 DAG 测试 |
| `docs/` | 架构、数据下载及各分析分支说明 |

总体设计和数据契约见 [架构说明](docs/architecture.md)，元数据来源与校验摘要见 [metadata README](metadata/README.md)。
