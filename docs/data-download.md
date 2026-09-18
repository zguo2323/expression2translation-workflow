# 第一阶段之后：实际数据下载交接

仅下载 `metadata/source/SRR_Acc_List.txt` 中的 8 个 runs。仓库不分发 raw reads。
`metadata/source/GSE203147_RAW.tar` 是作者已处理计数表，不能代替本项目的 raw FASTQ。

## 存放位置

从仓库根目录看：SRA archive 可放 `data/sra/`，压缩后的 FASTQ 放 `data/raw/`。
这些目录已被 Git 忽略。按当前 `config/samples.tsv`，最终应有 12 个文件：

| sample_id | run | FASTQ（位于 `data/raw/`） |
| --- | --- | --- |
| young_rna_1 | SRR19240769 | `SRR19240769_1.fastq.gz`、`SRR19240769_2.fastq.gz` |
| young_rna_2 | SRR19240768 | `SRR19240768_1.fastq.gz`、`SRR19240768_2.fastq.gz` |
| middle_rna_1 | SRR19240772 | `SRR19240772_1.fastq.gz`、`SRR19240772_2.fastq.gz` |
| middle_rna_2 | SRR19240771 | `SRR19240771_1.fastq.gz`、`SRR19240771_2.fastq.gz` |
| young_ribo_1 | SRR19240765 | `SRR19240765.fastq.gz` |
| young_ribo_2 | SRR19240764 | `SRR19240764.fastq.gz` |
| middle_ribo_1 | SRR19240770 | `SRR19240770.fastq.gz` |
| middle_ribo_2 | SRR19240766 | `SRR19240766.fastq.gz` |

下载端可能为 single-end 使用 `_1` 后缀。保留 SRR 对应关系后，可统一文件名到本表，或修改 samples.tsv 为实际路径再校验。
`.fastq.gz` 必须是真正 gzip 压缩的 FASTQ，不能仅把未压缩文件改扩展名。
不要根据文件名把 Ribo 推断为 RNA，也不要将 RNA 两端拼接成单个输入。

## 需要保留的信息

- 每个下载文件的来源 URL/获取工具与版本、获取日期、文件大小。
- 下载端提供的 checksum 及其算法（如果有）；另记录本地 SHA-256。
- 如果从 SRA 转换，保留转换/拆分/压缩命令与日志；下载 MD5 不等于重压缩文件的 MD5。

来源库：[SRA study SRP375616](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP375616)、[GEO GSE203147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203147)。
当前使用 SRA Toolkit 3.4.1；下面的下载脚本使用本机 `prefetch --help` 核实过的参数。

## 可靠的命令行入口

不要同时对同一 run 启动第二个下载。
从仓库根目录运行，先检查将使用的路径（不启动网络请求）：

```bash
.venv/bin/python scripts/download_reads.py --dry-run
```

只续传/下载一个 run：

```bash
.venv/bin/python scripts/download_reads.py --run SRR19240765
```

按 accession 清单依次下载（某个失败就停止）：

```bash
.venv/bin/python scripts/download_reads.py
```

若要下载后立即做 SRA 校验、FASTQ 转换、gzip 压缩及 FASTQ 完整性验证：

```bash
.venv/bin/python scripts/download_reads.py --convert --threads 4
```

脚本捕获 Ctrl+C/SIGTERM，停止自己的子进程并退出整个循环；不会跳到下一个 SRR。
已有单层或双层 accession 目录会被识别；新下载使用 `data/sra/SRR.../`。
旧命令在本机产生的是 `data/sra/SRR.../SRR.../`，新脚本使用实际包含 `.sra` 的目录调用转换工具，避免将外层目录误当作 archive。
遇到已有 SRA lock 会停止并要求检查，不强制删除锁或已有下载。

日志追加写入 `logs/download/`，每个 run 有状态 JSON。完整转换产物另有 `.completed.json`，包含路径、数量、大小、SHA-256；重跑校验摘要一致后跳过。
已有 FASTQ 没有完成记录时拒绝覆盖；下载完成不等于转换完成。
转换在 `data/tmp/` 的独立目录完成；成功验证才移动压缩产物到 samples.tsv 路径。
失败的临时输出保留便于诊断。RNA 中额外的 unpaired FASTQ 保留在对应 staging 目录，状态记录包含位置，不作为 PE 输入。
压缩会短暂同时保留未压缩和压缩文件，磁盘空间仍需关注。

测试：模拟网络工具验证了错误退出、中断传播、末行无换行、双层目录识别及转换行为；本脚本尚未用完整真实 SRA 完成端到端下载/转换。

## 网络受限时的 4-run 最小真实演示

完整研究设计仍是 8 runs；为了先打通真实数据路径，可以只取每个 condition、每个 assay 各一个 library：

| condition | assay | sample | run | archive bytes |
| --- | --- | --- | --- | ---: |
| Young | RNA | young_rna_2 | SRR19240768 | 1,922,068,149 |
| Middle | RNA | middle_rna_1 | SRR19240772 | 951,523,787 |
| Young | Ribo | young_ribo_1 | SRR19240765 | 1,109,220,705 |
| Middle | Ribo | middle_ribo_1 | SRR19240770 | 402,176,450 |

选中 archive 合计 4,384,989,091 bytes（约 4.08 GiB），是完整 8-run archive 的 42.2%。
清单位于 `config/minimal-demo-accessions.txt`。先看计划：

```bash
.venv/bin/python scripts/download_reads.py \
  --accessions config/minimal-demo-accessions.txt \
  --convert --threads 4 --dry-run
```

确认后执行相同命令并去掉 `--dry-run`。脚本只处理该清单中的 4 条；已有且完成摘要一致的 run 会在复核后跳过。

分析时分别使用 `config/rnaseq-minimal.yaml` 和 `config/riboseq-minimal.yaml`。两份配置仍保留真实 library type、strand 和 offset 门禁，必须在真实 QC/方法核实后填写。RNA 最小配置固定 `run_deseq2: false`，因为每个 condition 只有一个 library；它只能产出定量和描述性比较，不能估计生物学变异、报告可靠 p 值/FDR 或声称差异表达。Ribo 与最终 TE 同样只作描述性展示。

最小演示通过后，可以续传其余 4 条并切回 `config/rnaseq.yaml`、`config/riboseq.yaml`，不需要推倒已有结果。

## 空间与后续入口

archive 合计约 10.40 GB（9.68 GiB），不是 FASTQ 解压大小。交接文档的 60–80 GB 预留量是初步估计，需结合下载方式和可用磁盘核实。
下载完成后，下一阶段先做真实 reads 完整性/配对校验和方法核实，再锁定 reference、adapter/barcode、strandedness、offset。
当前的 metadata 校验不会检查 FASTQ 内容，所以下载后再次出现 `status=passed` 也不代表 reads 已验收。
