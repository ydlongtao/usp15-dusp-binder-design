# ESM3 阶段 B：OVO/Nextflow 评估适配器

阶段 B 建立一个独立的 OVO-facing descriptor/refolding 步骤，不修改 OVO Python 环境，也不替代 AF2、ProteinQC 或 USP4/USP11 反筛。

## 文件

- `esm3_eval.nf`：Nextflow DSL2 工作流，单 FASTA 输入，容器化运行。
- `scripts/run_esm3_eval.py`：生成 JSON、CSV 和 ESM3 structure-track PDB。
- `results/esm3_phase_b/`：本地运行结果目录，不提交大模型输出。

## 运行

```bash
nextflow run esm3_eval.nf \
  --input_fasta results/candidates.fasta \
  --output_dir results/esm3_phase_b \
  --adapter_script scripts/run_esm3_eval.py \
  --esm3_site /DATABANK/users/hflt/ovo/esm3_phase_a/site \
  --esm3_cache /DATABANK/users/hflt/ovo/esm3_phase_a/cache
```

在服务器上，`--esm3_site` 和 `--esm3_cache` 必须指向本地权重/依赖目录；不得把权重或令牌写入 Git。

## 输出语义

- `esm3_eval.json`：模型、设备、版本、种子、耗时、sequence/structure track 状态和序列摘要。
- `esm3_metrics.csv`：OVO 结果处理器可读取的扁平描述符表。
- `esm3_structure.pdb`：ESM3 structure-track 输出，仅作为正交结构诊断。

ESM3 结果不得解释为 iPAE、binder RMSD、pLDDT、ddG、KD 或结合自由能。阶段 B 的正式候选仍必须通过原有 AF2 正向门控和 USP4/USP11 反筛。

## 远端单候选验证（2026-08-20）

- 运行时：Nextflow 26.04.6 + `ovo-esm:latest` 基础镜像 + 独立 ESM 3.2.1 site 目录。
- 设备：Tesla V100-SXM2 32 GB，CUDA 12.4 PyTorch 2.5.1。
- 输入：一个 76 aa USP15 binder FASTA。
- sequence-track：通过；structure-track：通过。
- 输出：`esm3_eval.json`、`esm3_metrics.csv` 和 `esm3_structure.pdb`。
- Nextflow 发布目录：`/DATABANK/users/hflt/ovo/esm3_phase_a/results/phase_b_nf/candidates/`。

这证明阶段 B 的容器化适配器和 Nextflow 发布路径可运行；尚未对全部候选批量运行，也没有将 ESM3 结果用于硬性筛选。

## 阶段 B.2：批量评估

使用 `scripts/split_fasta_records.py` 将候选 FASTA 拆分为单候选输入，再以 `maxForks 1` 运行 `esm3_eval.nf`。每个候选独立输出 JSON、CSV 和结构 PDB，随后使用 `scripts/summarize_esm3_batch.py` 聚合结果。批量结果仍然只作为正交 descriptor，不自动晋级候选。
