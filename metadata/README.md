# 数据来源与选样

`source/` 保留用户提供的原始文件，第一阶段未重新下载或修改这些文件。
`source_manifest.json` 记录本地文件大小与 SHA-256；原始获取日期未知，填 `null`。
URL 用于定位来源，不表示已经验证过本地文件的下载历史。checksum 是本地快照校验，不能证明上游真实性。

| 文件 | 用途 |
| --- | --- |
| `SRR_Acc_List.txt` | 本项目选中的 8 个 SRR |
| `SraRunTable.csv` | 8 个 runs 的 layout、BioSample、GSM、年龄组、大小等原始证据 |
| `GSE203147_family.soft.gz` | GEO 样本标题、实验设计和处理方法；含未选中的 Aged 样本 |
| `GSE203147_series_matrix.txt.gz` | GEO series 元数据 |
| `GSE203147_RAW.tar` | 11 个作者处理后的 `.txt.gz` 计数表，约 430 KiB；**不是 raw reads/FASTQ** |

归档很小，可作为来源证据入 Git；当前校验只验证其 checksum，不把作者计数表作为工作流产出。
后续 FASTQ、SRA、BAM、reference/index 与 results 均使用被忽略的目录。

## 人工选样契约

仅选择 Young/Middle，每组各 2 个 RNA 和 2 个 Ribo library，原始映射见 `metadata/study_samples.tsv`；`config/samples.tsv` 默认保留全部 8 runs，也允许合法活跃子集。
Ribo 的 SRA `Assay Type=OTHER` 保留在来源文件；项目 `assay=riboseq` 是人工根据 GEO 标题和描述确认的。
校验器交叉检查 GEO title 的 assay/condition/replicate、GSM 与 BioSample/SRX 关系，以及 CSV 中的 SRR、GSM、物种和研究编号。
SOFT 中 Aged 记录保留作为来源，不属于当前 study 的 Young/Middle 8-run 设计。
通用自有 FASTQ 使用 `source_validation: local`，可无 accession；只有显式 `geo_sra` 模式要求这些来源证据。详见 [输入指南](../docs/data-input.md)。

8 个 archive 合计 **10,395,980,191 bytes = 10.40 GB ≈ 9.68 GiB**；合计 **22,453,938,457 bases ≈ 22.45 Gb**。
archive 大小不等于解压后的 FASTQ 大小。

## 方法线索，待第二阶段核实

SOFT 记录了 SacCer3、TruSeq Stranded mRNA kit，以及 adapter、样本 barcode sorting、random barcode trimming 和 STAR 处理。
这些线索不足以锁定本项目 reference release、Salmon library type、adapter 序列或 Ribo barcode/offset 参数。
RNA/Ribo 使用不同 BioSample；相同 replicate 编号不构成跨 assay 配对证据。

## 更新来源

来源校验失败时先查明文件是否被意外改动。只有有意更新来源时才更新 manifest 的大小和 SHA-256，并在版本历史中记录原因。
不要为了让校验通过而无条件重新生成 checksum。
