# 输入自己的酵母 RNA/Ribo 数据

当前研究实例默认分析 Young/Middle 全部 8 runs（4 个 RNA PE、4 个 Ribo SE），原始映射保存在 `metadata/study_samples.tsv`。`config/samples.tsv` 只表示当前活跃 libraries；数量不限定为 4 或 8。工作流仍限于酵母 RNA paired-end、Ribo single-end 与既有 reference 契约。

## 一行一个生物学 library

必需列为 `sample_id, assay, condition, replicate, layout, fastq_1, fastq_2`，文件使用 TAB 分隔，不能用逗号。参考 [无 accession 的六样本示例](../config/samples-local.example.tsv)。

- `sample_id` 唯一，使用小写字母开头的小写字母、数字和下划线。
- `condition` 使用字母开头的字母、数字、点、下划线或连字符，须在 `design.conditions` 中。
- `replicate` 是 condition × assay 内独立生物学重复的标识，可以是 `1`、`7`、`donor_A`；不要求连续编号。该三元组必须唯一。
- 技术重复、拆分 lanes 应先在工作流外按既定方法整理为同一 library，记录处理来源；不能修改 replicate 编号冒充独立生物学重复。校验器能拒绝重复身份，不能从编号证明生物学独立性。
- RNA 固定 `rnaseq / PAIRED`，必须填写两端；Ribo 固定 `riboseq / SINGLE`，末列保留为空。路径必须位于 `paths.raw` 下，以仓库为根的相对路径，不允许 `..`、绝对路径、逃逸 symlink 或重复路径，后缀 `.fastq.gz` / `.fq.gz`。
- `run_accession`、`geo_accession`、`biosample` 可省略整列或留空。填写时检查格式；run 支持 SRR/ERR/DRR，非空 run/GSM 唯一。同一 assay 内 BioSample 不能重复计数；跨 assay 可复用，但不因此启用配对模型。
- MVP 暂不支持 `pair_id` 列或配对模型；重复编号相同不代表 RNA/Ribo 配对。

不再要求完整 condition × assay × replicate 笛卡尔积。旧 `design.replicates` 可保留供配置兼容，但不参与约束，建议移除。

## 两种来源模式

`source_validation: local` 只检查通用输入契约，不需要 GEO/SRA、下载记录或 source manifest。外部获取的 FASTQ 与项目下载器的 FASTQ 使用同一 reads 校验。无来源统计时 `archive_bytes`、`bases` 及总量为 JSON `null`；可选来源字段进入 SQLite 时使用 SQL NULL。

`source_validation: geo_sra` 是当前研究的显式证据适配器，需要 `study_samples`、`source_manifest`、`sources` 和研究编号。完整 study 表仍与原始 GEO/SRA 映射、标题及 checksum 交叉验证；活跃表可以是其子集，路径可按实际存放位置调整，但样本身份不能变。来源清单无需按子集删改，也不决定当前下载范围。其他研究若不符合这一 GEO 标题/CSV 格式，应使用 local 模式并自行保留来源记录。

## 新用户入口

1. 将自己的 FASTQ 放到 `data/raw/`，复制并编辑 [local.example.yaml](../config/local.example.yaml) 和六样本示例，填实际路径与条件。示例不附带 reads。保留 `source_validation: local`；示例显式清空研究编号与来源路径，防止 Snakemake 合并默认配置时保留当前研究身份。
2. 默认 RNA/Ribo/QC 的 `sample_ids: []` 自动选取活跃表中的相应样本。非空列表是显式子集；未知、重复或错误 assay/layout 会被拒绝。分析阶段 QC 的非空选择必须与分析样本完全一致，否则清空让工作流统一选择。
3. 按 [reference 文档](rnaseq.md)准备并验证固定 Ensembl 115 / R64-1-1 资源；不要覆盖已有 reference。
4. 先做 metadata dry-run，再做 raw QC dry-run。下面仅构建 DAG，不启动分析：

```bash
.venv/bin/python -m src.validation.validate --config config/local.example.yaml
.venv/bin/snakemake --snakefile workflow/Snakefile --cores 1 \
  --configfile config/local.example.yaml --dry-run
.venv/bin/snakemake --snakefile workflow/Snakefile --cores 2 \
  --configfile config/local.example.yaml --config stage=qc --dry-run
```

5. raw QC 的 FASTQ 全量校验记录 gzip/四行格式、PE IDs/顺序/数量、SHA-256，不要求 `.completed.json`。已有外部文件不要交给 `--convert` 重新转换。可单独运行统一校验并保留报告：

```bash
.venv/bin/python -m src.validation.reads \
  --reads data/raw/my_library_1.fastq.gz data/raw/my_library_2.fastq.gz \
  --output results/input-check/my_library.json
```

6. 方法和真实 reads 确认后，在 QC/RNA/Ribo 配置中填写 adapter、barcode handling、library type、strand、offset 及证据；切换 trim 后才允许定量。各分支实际运行前先通过对应 dry-run。不能复制 synthetic 的参数。

## 统计与子集

metadata/QC 支持一个或多个配置条件的合法子集。当前 RNA（包括描述性导入）需要所选样本恰好覆盖一个两组 contrast；integration 要求两组各同时有 RNA 和 Ribo。更多条件可在样本表中登记，但须显式选两组进入 RNA/integration，不支持任意多条件统计模型。

`run_deseq2: true` 需要两组各至少两个独立生物学重复、`design_confirmed: true` 和非空 `design_evidence`，还要通过输入与方法门禁；关闭时允许 counts/TPM 描述性输出，不生成 DE/PCA。当前 8-run 已依据论文与 GEO 分组确认 assay 内每组两个生物学重复；n=2 功效仍有限，且不增加跨 assay 配对项。TE 始终为 condition-level 描述，`differential_te: false`；不扩展 Ribo 单独差异表达。

四样本示例 [samples-minimal.tsv](../config/samples-minimal.tsv) 可直接替代活跃表，无需改变 study/source 文件：

```bash
.venv/bin/snakemake --snakefile workflow/Snakefile --cores 1 \
  --config samples=config/samples-minimal.tsv --dry-run
.venv/bin/python scripts/download_reads.py --samples config/samples-minimal.tsv --dry-run
```

进入分析时同时选择 `rnaseq-minimal.yaml` / `riboseq-minimal.yaml`，前者关闭 DESeq2，二者仍保留方法门禁。此子集只是可选示例，当前真实实例默认完整 8-run。

## 旧配置与结果

旧配置新增 `source_validation: geo_sra` 和 `study_samples: metadata/study_samples.tsv` 即可迁移当前研究；自有数据选 local 模式。SQLite schema v2 允许 NULL 来源、文本 replicate。旧 catalog 在 dry-run 和写入前被拒绝，保留原库，将 `paths.results` 改为新的相对目录重建；不自动迁移、删除或覆盖旧库。

Synthetic smoke 自行生成无 GEO/SRA 身份的固定测试设计，不依赖用户的活跃样本表；其结果只验证工程行为。
