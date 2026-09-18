# E2T 当前架构与实施契约

本文定义 E2T MVP 的当前实施范围、数据契约和统计边界。

## 范围

酵母 bulk RNA-seq + Ribo-seq 联合分析 MVP，演示数据为 GSE203147 的 Young/Middle 8 个 runs。
不扩展 DNA assembly、单细胞、蛋白组、Web/Agent 或复杂新 ORF 分析。
配置管理分组，Python/Snakemake 不把 aging 条件写死；GEO/SRA 证据读取针对当前数据来源格式。

## 分阶段数据流

```text
当前实现：
config + samples.tsv + source manifest + GEO/SRA metadata
  → 只读 preflight（含 dry-run）
  → validate_metadata → results/validation/metadata.json

阶段 2–3 已实现并经 synthetic 真实工具验收：
统一 genome/GTF → transcriptome/rRNA/tx2gene → Salmon index/quant
  → tximport → DESeq2、PCA/相关性、provenance
统一 reference → Ribo QC → rRNA removal → transcriptome mapping
  → P-site/CDS gene counts、periodicity、metagene、provenance

当前完整 MVP 数据流：
真实 FASTQ + versioned reference + 已确认方法参数
  → analysis input validation
  ├─ RNA: FastQC/fastp → Salmon → tximport gene estimates → DESeq2
  └─ Ribo: QC/adapter/barcode → rRNA removal → mapping → P-site/QC/counts
  → condition-level gene integration → SQLite + provenance + Markdown/HTML report
```

当前默认 DAG 只有 `validate_metadata` 和 `all`。新增 `stage=qc` 分支可运行输入门禁、FastQC、可配置 fastp 和 MultiQC，详见 [QC 使用说明](qc.md)。下载使用独立显式命令，不在默认 DAG 自动启动。
验收通过 metadata 不意味着 raw reads 或统计设计已通过验证。

## 数据契约

`config/samples.tsv` 是项目的样本事实来源：

| 字段 | 约束 |
| --- | --- |
| `sample_id` | 唯一的项目内 library ID |
| `run_accession` | 唯一 SRR，与选中 accession 清单完全一致 |
| `geo_accession`、`biosample` | 当前数据一 library 对应一个 GSM/BioSample，与 CSV/SOFT 交叉校验 |
| `assay` | `rnaseq` / `riboseq`，后者由 GEO 证据人工标注 |
| `condition`、`replicate` | 在配置设计中且每个组合唯一；重复编号不表示跨 assay 配对 |
| `layout` | RNA PAIRED、Ribo SINGLE |
| `fastq_1/2` | `paths.raw` 下的相对压缩 FASTQ 路径；PAIRED 两个不同路径，SINGLE 第二列为空 |

样本表所有 FASTQ 路径必须唯一，禁止绝对路径和 `..` 路径穿越。
source manifest 校验所有列出的本地来源文件，包括完整 GEO 元数据中未选择的 Aged 记录。
源数据不做过滤改写；只有 samples.tsv 和 accession 清单决定本次分析范围。

## 配置与质量门禁

`stage=metadata`：允许 reference/预处理参数为 null，禁止启用 differential TE；仅访问元数据。
未来的分析输入门禁必须检查实际 FASTQ、格式/完整性、下载 checksum、参考成套一致性、gene ID 映射及每个处理阶段所需参数。
read QC 才能验证 strandedness、adapter/barcode 处理与 Ribo read-length、periodicity、metagene/offset。
不会用固定 offset 或作者的 assembly 名字代替实证检查和 release 锁定。

校验器对错误输入返回非零退出码；Snakemake 解析阶段同样调用此逻辑。
因此已有报告也不能掩盖新出现的坏样本表或来源文件改动。
报告以临时文件原子替换；源码、schema、样本、来源文件或有效配置变更会触发重跑。
所有命令从仓库根目录执行，当前 smoke test 使用运行 Snakemake 的同一 Python 环境。

## 联合分析与统计边界

RNA/Ribo 的 BioSample 不同，严格配对尚未确认。首先按 condition 汇总，并按统一 reference 的 gene ID 联结。
不将相同 replicate 编号当作 `pair_id`；SQLite integration 表设计使用 condition/contrast 与方法版本，不能强制非空 pair_id。
描述性 TE 需明确 assay 内归一化、聚合顺序、计数区域、低表达过滤和 pseudocount。
不直接把不同文库深度的原始 RNA/Ribo counts 比值作为可比较的 TE。
RNA Salmon estimated counts 与 Ribo site counts 的定量含义须分别记录；模型与输入准备在分析分支中再验证。
assay 内 DESeq2 和 differential TE 是不同问题，本阶段不产出显著性结论。

## 后续产物契约

- `results/qc/`：FastQC/fastp/MultiQC、read-length、去污染、mapping 等指标。
- `results/rnaseq/`：transcript/gene estimates、归一化 counts、PCA，设计允许时 DE。
- `results/riboseq/`：gene/site counts、density、periodicity、metagene 与 offset 证据。
- `results/integration/`：condition-level gene table、描述性 TE 和方法说明。
- `results/catalog.db`：samples、workflow runs、references、genes、QC、RNA/Ribo metrics、integration、artifacts。
- `results/provenance/`、`results/report/`：有效参数、工具/reference 版本、SHA-256、日志/DAG、Markdown/HTML 报告。

SQLite 只保存结构化结果和产物索引，不存储大 FASTQ/BAM；要求外键、唯一性、幂等写入和完整性检查。
这些目录/数据库是后续交付目标，第一阶段不创建虚假的分析结果。

阶段 2–4 更新：reference 已固定 Ensembl release 115 / R64-1-1，RNA、Ribo、描述性 TE、SQLite 和报告已通过完整 synthetic 端到端验收。真实运行仍受预处理、文库方向和 Ribo offset 证据门禁约束；4-run 最小演示不产生显著性结论。

## 第一阶段验收

1. 来源文件完整、8 行人工样本表与 GEO/SRA 一致。
2. 无 FASTQ/reference 时 metadata preflight 与 dry-run 成功。
3. 实际运行只生成校验报告和日志，记录 effective config 和来源/代码/环境摘要。
4. 错误 SRR、layout、重复、缺失来源、checksum、跨 assay 设计等输入被阻止。
5. 无修改重跑不重复生成报告；大数据路径被 Git 忽略。
