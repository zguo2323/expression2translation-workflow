# QC 使用与边界

## 当前功能

- `stage=metadata` 仍为默认值，未下载 FASTQ 时可运行。
- `stage=qc` 读取 `config/qc.yaml`，仅运行其中的 `sample_ids`；空列表表示全部 8 个。
- `mode=raw`：逐样本完整 FASTQ/gzip 校验 → raw FastQC → MultiQC。
- `mode=trim`：上述步骤 + fastp → trimmed FASTQ 校验 → trimmed FastQC → MultiQC。

paired-end 校验完整扫描两端，检查记录格式、非空、序列/质量长度、read ID/顺序及数量一致，最后计算压缩文件 SHA-256。
当前 FASTQ 解析契约是常见的每条 4 行 ASCII FASTQ，序列只接受 A/C/G/T/N；异常文件会失败，不尝试自动修复。
这些检查不代表 reads 已通过生物学 QC，也不判断平台的质量编码或 adapter/strandedness。

## 第一个 run 下载完成后

将 `config/qc.yaml` 中的样本选择改为已经完整转换的样本，例如：

```yaml
mode: raw
sample_ids: [young_ribo_1]
synthetic: false
```

使用已创建的项目 QC Conda 环境和现有 `.venv` 中的 Snakemake：

```bash
conda run --no-capture-output --prefix .conda/qc .venv/bin/snakemake --cores 2 --config stage=qc --dry-run
conda run --no-capture-output --prefix .conda/qc .venv/bin/snakemake --cores 2 --config stage=qc
```

其他机器先 `conda env create --prefix .conda/qc --file envs/qc.yaml`。
具备 Conda >=24.7.1 时，也可激活 workflow 环境后使用 `snakemake --cores 2 --config stage=qc --use-conda`，让规则声明的 `envs/qc.yaml` 自动部署。本机 Conda 23.9 使用前述显式环境方式。
本次本地实测使用前一种共享 QC 环境方式。

预处理参数的执行事实来源是 `config/qc.yaml` 的 `policies`；原始 config 中的 `preprocessing` 待定项保留为后续方法调查清单，不会自动转成 fastp 参数。
只有需要运行的 assay 必须设置 `confirmed: true` 并提供 `evidence`，明确 adapter 序列或已确认无需去 adapter、barcode 策略、前端裁剪长度、最短长度和质量过滤阈值。
当前只支持无 barcode 处理或固定前端裁剪，未知 barcode sorting/UMI 方案必须先核实和扩展实现。
禁用 fastp 的自动 poly-G 裁剪，避免工具依据 read header 隐式改变行为；不启用去重、合并或纠错。
strandedness 和 P-site offset 不影响 raw QC，后续定量步骤另设门禁。

## 输出

`results/qc/raw/` 或 `results/qc/trim/` 下包含：

- `validation/`：完整输入/预处理后 FASTQ 校验报告及 SHA-256。
- `fastqc/raw/`、`fastqc/trimmed/`：按样本/read 独立的 HTML/ZIP。
- `fastp/`：预处理 FASTQ、JSON、HTML（trim 模式）。
- `multiqc/multiqc_report.html`、`multiqc/multiqc_report_data/`：本次显式选择的样本与模块结果。
- 每个工具旁的 provenance JSON：版本、实际命令、配置参数、时间和 synthetic 标识。

MultiQC 只扫描当前 DAG 声明的结果，避免把其他样本子集或旧模式输出混入报告。
FASTQ 被全部过滤为空时门禁失败，需检查策略；不会静默推进下游。

## 合成数据验收

```bash
conda run --no-capture-output --prefix .conda/qc .venv/bin/python scripts/smoke_qc.py --directory results/new-synthetic-smoke
conda run --no-capture-output --prefix .conda/qc .venv/bin/python scripts/smoke_qc.py --mode raw --directory results/new-raw-smoke
```

脚本复制工作流到全新隔离目录，只生成 1 个 RNA PE 和 1 个 Ribo SE 的 synthetic FASTQ，各 200 个 records/pairs。
元数据仍用于测试现有样本接口；合成序列不能冒充对应 GEO 样本。MultiQC 标题和 provenance 标明 SYNTHETIC TEST。
检查 dry-run、真实工具输出、已知前缀/adapter 去除、reads/base 数量、MultiQC 解析模块及无修改重跑。
已有 smoke 目录不会覆盖，请指定新目录。
GitHub Actions 配置已提供；本地通过不代表 GitHub 上已经执行 CI。

## Reference 契约

`workflow/schemas/reference.schema.json` 和 `src/validation/reference.py` 定义 release/provider/reference ID、四种资源、来源/获取时间、衍生命令和 SHA-256。
当前仅检查文件身份与完整性，不检验 genome/GFF/transcriptome 的生物学兼容性，且 QC 不需要 reference。
阶段 2 已新增 `validate_rna_reference` 并锁定正式 reference：验证 contig、exon 坐标、负链拼接、transcript/rRNA 序列及 tx2gene 一致性，详见 [RNA 说明](rnaseq.md)。

工具行为参考：[fastp](https://github.com/OpenGene/fastp)、[MultiQC CLI](https://docs.seqera.io/multiqc/getting_started/running_multiqc)。
