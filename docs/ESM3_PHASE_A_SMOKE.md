# ESM3 阶段 A：独立 Smoke Test

## 输入

默认输入为 `results/candidates.fasta` 中的第一条 USP15 binder 序列。也可以通过 `--fasta` 指定单条 FASTA 文件。

## 检查内容

1. Python、PyTorch、ESM 包和 Hugging Face 运行时版本。
2. CUDA/CPU 设备和显存可见性。
3. `esm3-sm-open-v1` 模型加载。
4. sequence-track 生成。
5. structure-track 生成（若设备和运行时允许）。
6. 输出文件、耗时和可复现失败原因。

## 执行方式

```bash
python scripts/run_esm3_phase_a_smoke.py \
  --fasta results/candidates.fasta \
  --output-dir results/esm3_phase_a \
  --model esm3-sm-open-v1
```

脚本默认只处理第一条序列、单个随机种子，不启动 OVO 或 Nextflow。模型权重写入本地缓存，不提交 Git。

## 结果解释

- `passed`：模型和 sequence-track 成功；若 structure-track 也成功，则满足阶段 A 完整通过条件。
- `blocked`：依赖、许可证、模型下载、显存或运行时阻断；保留诊断信息，不伪造通过。
- 本阶段不产生 binder 结合、亲和力、选择性或动力学结论。

## 本地执行记录（2026-08-18）

- 模型：`esm3-sm-open-v1`，本地权重已成功加载。
- 输入：`results/candidates.fasta` 第一条序列，76 aa。
- 运行时：Python 3.14.5、PyTorch 2.13.0、ESM 3.2.1，macOS arm64。
- sequence-track：CPU、1 个生成 step，成功，耗时约 43.7 s。
- Apple MPS：当前 ESM3 实现返回 `Unsupported device type: mps`，因此标记为运行时限制，不宣称 MPS 通过。
- structure-track：尚未在本地 CPU 上运行；应在 Minerva H100 或兼容 CUDA 节点上继续验证。

本次结果只证明 ESM3 可以在本地环境完成模型加载和最小序列生成，不代表任何 USP15 结合、亲和力或选择性结论。
