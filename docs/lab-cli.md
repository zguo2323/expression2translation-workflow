# 实验室项目 CLI（第一阶段）

`e2t` 是既有 Snakemake 分析引擎外的一层项目管理 CLI。它用于保存多个课题的项目状态和最终结果索引；不读取、复制或上传 FASTQ，也不改变现有 `results/catalog.db` 的单项目分析职责。

## 安装

```bash
python -m pip install -e .
# 无网络且环境缺少 wheel 时：python setup.py develop
e2t --help
```

默认 storage root 为当前目录下的 `e2t-lab-data`；生产环境建议显式使用共享 Linux 文件系统路径，或设置 `E2T_HOME`。

## 项目结构

```text
<root>/
  registry.sqlite
  projects/<project-id>/
    project.json
    inputs/ raw/ work/ results/ logs/ provenance/
```

`raw/` 用于只读原始输入，`work/` 用于可重建中间文件，`results/`、`logs/` 与 `provenance/` 应作为项目交付和审计记录保留。共享 reference 和容器目录尚未由第一阶段 CLI 创建或管理。

## 命令

创建项目。project ID 仅允许小写字母、数字和连字符，且不能以连字符开始：

```bash
e2t init aging-yeast-2026 --root /srv/lab-e2t --title "Yeast aging" --owner alice
```

查看全部项目或一个项目的完整状态历史与 artifacts：

```bash
e2t list --root /srv/lab-e2t
e2t status aging-yeast-2026 --root /srv/lab-e2t --json
```

登记分析进度。`stage` 可取 `metadata`、`qc`、`rnaseq`、`riboseq` 或 `integration`；`status` 可取 `created`、`running`、`succeeded`、`failed` 或 `archived`：

```bash
e2t set-status aging-yeast-2026 --root /srv/lab-e2t \
  --stage integration --status succeeded --message "Report reviewed"
```

登记项目内已存在的最终文件。CLI 会计算 SHA-256，并拒绝项目目录以外的文件：

```bash
e2t add-artifact aging-yeast-2026 --root /srv/lab-e2t \
  --path /srv/lab-e2t/projects/aging-yeast-2026/results/report/report.html \
  --kind html_report
e2t check --root /srv/lab-e2t
```

## 数据与并发边界

registry 采用 SQLite，适合十余人规模下的短事务元数据操作。它保存项目 ID、所有者、状态历史、文件相对路径、大小和 SHA-256；不保存 FASTQ、BAM 或完整 gene-level 表。第一阶段没有任务队列、权限系统、远程对象存储上传、Docker daemon 或 Slurm 调度。后续版本应以共享 reference 的只读目录、容器 digest provenance 和 Slurm profile 补齐多用户计算环境。
