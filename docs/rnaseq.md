# RNA 分支与固定参考

当前实现：已确认的 RNA 预处理 → FASTQ 门禁 → Salmon decoy-aware index/quant → tximport → 可选 DESeq2、PCA/相关性、结果摘要和 SHA-256 provenance。
这是 assay 内 Middle vs Young 的 RNA 分析，不涉及跨 assay pairing 或 differential TE。

## 固定参考与方法依据

选择 Ensembl **release 115 / R64-1-1**，从同一 release 的 genome FASTA 和 GTF 派生资源。固定 manifest 位于 `metadata/reference-ensembl-115.json`，原始及派生大文件位于被忽略的 `resources/reference/ensembl-115-r64/`。
本机已实际获取并校验：17 个 contigs、7,127 个转录本/基因、24 个 rRNA 转录本。

选择此成套 GTF 资源是为了直接获得一致的 exon/transcript/gene ID 关系；不混用 SGD 另一 release 的注释。
[Ensembl 的组装说明](https://jun2026-fungi.ensembl.org/Saccharomyces_cerevisiae/Info/Annotation/)将 R64-1-1 标为 SGD S288c reference。
GEO SOFT 的实验菌株为 BY4743，原文记录 Assembly: SacCer3；本项目使用标准 S288c reference，不声称与作者具体 annotation release 或株系变异完全一致。

参考的生成规则：GTF 1-based inclusive 坐标转为 Python 切片，按基因组位置拼接 exons，负链整体反向互补；保留准确 transcript/gene ID，不自动删除版本后缀。
transcriptome 包含所有有 exon 注释的 transcript biotypes；rRNA 子集由 GTF biotype 提取。
门禁检查 contig/坐标、exon 重叠、重复 ID、跨基因冲突，以及 FASTA 序列、rRNA 和 tx2gene 与重建结果严格一致。
这证明文件彼此一致，不代表注释本身或所有 transcript isoform 都已由实验验证。

首次获取来源（已有文件时勿直接覆盖）：

```bash
mkdir -p resources/reference/ensembl-115-r64/source
curl -fL 'https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz' -o resources/reference/ensembl-115-r64/source/genome.fa.gz
curl -fL 'https://ftp.ensembl.org/pub/release-115/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz' -o resources/reference/ensembl-115-r64/source/annotation.gtf.gz
```

随后用 `scripts/prepare_reference.py` 生成派生文件；它拒绝覆盖已有输出目录。完整命令：

```bash
.venv/bin/python scripts/prepare_reference.py \
  --genome resources/reference/ensembl-115-r64/source/genome.fa.gz \
  --annotation resources/reference/ensembl-115-r64/source/annotation.gtf.gz \
  --output resources/reference/ensembl-115-r64/derived \
  --reference-id ensembl-115-r64 --provider Ensembl --release 115 \
  --genome-source https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
  --annotation-source https://ftp.ensembl.org/pub/release-115/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz
```

可独立核对仓库中已固定的 SHA-256 和序列一致性：

```bash
.venv/bin/python -c 'from src.validation.reference import validate_rna_reference; validate_rna_reference("metadata/reference-ensembl-115.json"); print("reference passed")'
```

manifest 时间记录为本次下载后的核验登记时间；上游未提供本项目使用的 SHA-256，这些摘要是本地计算的固定快照。`preparation.json` 记录生成代码摘要。

## 真实数据启动条件

`config/rnaseq.yaml` 明确选择 4 个 RNA libraries，RNA 阶段使用这些 ID 驱动 QC，不要求 Ribo FASTQ 已下载。
运行前需要：

1. 4 个 RNA libraries 都有完整两端 FASTQ。
2. `config/qc.yaml` 为 `mode: trim`，RNA policy 经确认且有方法证据；Ribo policy 可继续待定。
3. `library_type`/`library_evidence` 明确。允许 IU/ISF/ISR 或显式选择 A 自动推断；GEO 的 TruSeq Stranded mRNA 线索不自动等于已验证 ISR。
4. 启用 DESeq2 时，`design_confirmed: true` 且记录 assay 内生物学重复/分组证据；至少两组、每组两个 libraries。不同 assay 的重复不参与这个模型。

Salmon 使用完整 genome 作为 decoy，显式固定 k=31 和 keepDuplicates，保存实际版本/命令、mapping 和文库兼容性指标。
当前 mapping rate 门槛为 0.5，compatible fragment ratio 为 0.8，属于可配置工程门禁，需结合真实 reads 复核；不是通用生物学合格标准。
`A` 的推断结果保留在 Salmon metadata/lib_format_counts，仍需人工检查；门禁不把自动推断当作实验设计确认。
参数依据：[Salmon 文档](https://salmon.readthedocs.io/en/latest/salmon.html)。本项目固定实测 **Salmon 1.10.3**，没有随未固定的最新包变更索引实现。

## 运行环境与命令

已建立 `.conda/qc`、`.conda/salmon`、`.conda/rnaseq-stats`。新机器可按对应 `envs/*.yaml` 创建。
本机 Conda 23.9 低于 Snakemake 9.24 自动环境部署要求的 24.7.1，因此本次实测使用显式环境路径；未修改用户 base。

真实参数确认后，在仓库根目录运行：

```bash
export PATH="$PWD/.conda/qc/bin:$PWD/.conda/salmon/bin:$PWD/.conda/rnaseq-stats/bin:$PATH"
.venv/bin/snakemake --cores 4 --config stage=rnaseq --dry-run
.venv/bin/snakemake --cores 4 --config stage=rnaseq
```

具备满足版本要求的 Conda 的机器可用 `--use-conda` 自动部署规则环境，CI 配置采用此方式。本地此次未实测该部署方式。
R 使用 `--vanilla`，不加载个人启动配置，以减少全局 R 库对项目环境的影响。

## 结果及统计边界

`results/rnaseq/` 包含 index、每样本 quant.sf/Salmon metadata/provenance，以及 `gene_analysis/`：

| 文件 | 含义 |
| --- | --- |
| gene_estimated_counts.tsv | tximport 汇总的估计 fragment counts，允许小数，不是 TPM |
| gene_tpm.tsv、gene_average_length.tsv | gene abundance 与样本特异平均 transcript length |
| tximport.rds | 原始 tximport 对象 |
| gene_normalized_counts.tsv | DESeq2 归一化结果，仅含通过预过滤的基因 |
| differential_expression.tsv | baseMean、log2FoldChange、SE、stat、pvalue、padj；NA 保留 |
| pca.tsv、pca.pdf、sample_correlation.tsv | 使用 VST、blind=FALSE 的样本 QC |
| analysis.json、sessionInfo.txt | contrast、过滤、拟合类型、真实包版本、synthetic 标识 |
| provenance.json | quant/reference/tx2gene/脚本和输出产物的 SHA-256 |

统计路径为 `countsFromAbundance="no"` → `DESeqDataSetFromTximport`，由 DESeq2 使用长度信息；不拿 TPM 当 counts，也不手动重复归一化长度。
contrast 显式为 Middle / Young，正 log2FC 代表 Middle 更高；不强加 pair/batch 项，不做无证据的配对。
预过滤是 counts >= min_count 的样本数 >= min_samples。结果默认未经 LFC shrinkage，DESeq2 的 independent filtering 和 Cook's 处理产生的 NA 不填成显著值。
方法依据：[tximport 官方说明](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)。
`run_deseq2: false` 只汇总 gene counts/TPM/length，不产生 DE/PCA；当前样本选择仍需包含配置 contrast 两组。

## 合成数据与实测限制

```bash
.venv/bin/python scripts/smoke_rnaseq.py --directory results/new-synthetic-rna-smoke
```

脚本生成含正/负链及 intron 的 81 个 synthetic transcripts，4 个 RNA PE libraries。
其中 8 个基因增加、8 个减少、1 个 rRNA gene 不生成 reads；检查已知变化方向、零计数过滤、SHA-256 产物、无修改重跑及 library_type 变化后的重新调度。
synthetic fit_type 使用 mean，真实配置保留 parametric；实际 fitType 写入摘要。
合成数据中的显著基因数不能作为真实 GSE203147 分析结果，也不是灵敏度/假阳性率的全面基准。
真实 reads 方向、adapter、参考适配率、DE 结果仍待数据；GitHub Actions 仅已配置，未声称远端运行成功。
