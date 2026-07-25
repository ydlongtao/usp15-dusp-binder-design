# USP15 DUSP Binder Design

这是一个面向人 USP15 N 端 DUSP 结构域的计算型微型蛋白 binder 设计项目。工作流使用 OVO 提供的 RFdiffusion RFD1、LigandMPNN 和 AlphaFold2/ColabDesign 模块，通过 Nextflow + Docker 在 GPU 服务器上直接运行。

本项目仅覆盖计算设计、排序和导出，不包含湿实验，也不宣称候选一定抑制 USP15 酶活或阻断特定天然互作。

## 设计对象

- 目标蛋白：human USP15，UniProt `Q9Y4E8`
- 主结构：[3T9L](https://www.rcsb.org/structure/3T9L) chain A residues 6–134
- 界面参考：[6DJ9](https://www.rcsb.org/structure/6DJ9)
- 热点：`A50,A52,A53,A55,A57,A61`
- 设计长度：45–60 aa 和 61–80 aa
- RFdiffusion 权重：`Complex_base` 和 `Complex_beta`
- 反筛对象：USP4 [5CTR](https://www.rcsb.org/structure/5CTR) 与 USP11 [4MEL](https://www.rcsb.org/structure/4MEL)

四个设计池及全部阈值记录在 [`config/campaign.json`](config/campaign.json) 和 [`config/pools.tsv`](config/pools.tsv)。

## 参数优化与失败恢复

R1 的四个 100-backbone pilot 已完成，但首轮 Complex_base LigandMPNN + AF2 smoke 未通过既定 AF2 门控，因此不得直接扩量或放宽阈值。

后续保持 USP15 DUSP 靶点不变的分阶段优化方案，包括现有 Complex_beta 候选复核、生成热点子集、紧凑性过滤、scaffold-guided 小矩阵和扩量停止规则，记录在：

- [`docs/USP15_PARAMETER_OPTIMIZATION_PLAN.md`](docs/USP15_PARAMETER_OPTIMIZATION_PLAN.md)

R2 阶段 A 已完成：18/18 条 LigandMPNN + AF2 设计技术成功，但
0/18 通过既定 AF2 门控，且没有启动扩量。完整执行结果与阶段 B
停止条件已写回优化计划；不得覆盖 R1 参数和失败记录。

阶段 A 的独立链映射审计也已通过：18/18 个预测均保持 binder A /
target B，target CA RMSD 为 0.733–0.788 Å，独立复算的 binder RMSD
与 OVO 指标一致。因此阶段 A 的失败不是链号或 target-template 配置错误。

阶段 A 获得授权后，可使用独立矩阵执行已有 Complex_beta 候选的 18 条序列诊断：

```bash
cp config/r2_phase_a.tsv "$USP15_CAMPAIGN_DIR/config/"
cp scripts/run_r2_phase_a.sh scripts/summarize_r2_phase_a.py \
  "$USP15_CAMPAIGN_DIR/scripts/"

"$USP15_CAMPAIGN_DIR/scripts/run_r2_phase_a.sh"
```

脚本严格串行运行 6 个条件，输出
`r2/phase_a/reports/af2_metrics.csv` 和
`r2/phase_a/reports/phase_a_summary.json`，并明确禁止自动启动扩量。

R2 阶段 B 使用 [`config/r2_phase_b.tsv`](config/r2_phase_b.tsv) 的六条件、
每条件 50-backbone 矩阵。标准条件继续走 OVO Nextflow；scaffold-guided
条件使用 RFdiffusion 官方无-contig调用，再把原始 TRB 复制并规范化为 OVO
可读取的元数据。原始 TRB 和两次失败 smoke 都保留，不修改坐标。

R2 阶段 B/C 已全部完成。300 个 backbone 中有 20 个进入序列设计：
`B1=0, B2=3, B3=3, B4=4, S1=5, S2=5`。40 个 LigandMPNN 条件均
技术成功，共生成 120 条无 Cys 序列；120/120 条
`af2_model_1_multimer_tt_3rec` 预测完成，但 0/120 同时通过 iPAE、
target-aligned binder RMSD 和 binder pLDDT 门控。本轮状态为
`af2_gate_failed`，没有启动 1000-backbone 扩量，也没有放宽阈值。

```bash
"$USP15_CAMPAIGN_DIR/scripts/install_rfdiffusion_scaffold_checkpoint.sh"
"$USP15_CAMPAIGN_DIR/scripts/prepare_r2_phase_b_resources.sh"
"$USP15_CAMPAIGN_DIR/scripts/run_r2_phase_b_queue.sh"
```

六组 backbone 全部完成后，Phase C 自动从每组过滤结果中取最多 5 个
backbone，在 temperature `0.05` 和 `0.10` 下各生成 3 条无 Cys 序列，
并仅运行 `af2_model_1_multimer_tt_3rec`：

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_r2_phase_c.sh"
```

R3 已获授权转向 6DJ9 已知结合支架救援，完整决策树记录在
[`docs/USP15_R3_UBV_RESCUE_PLAN.md`](docs/USP15_R3_UBV_RESCUE_PLAN.md)。
原始 UbV、全 backbone LigandMPNN 和稳定 1UBQ 界面嫁接的最小 smoke
均未通过不变的 AF2 门控，因此这些结果不计为候选。当前分支使用 RFD1
对 3T9L–ubiquitin 复合物做低噪声 partial diffusion：

- `partial_T = 5, 10, 15`，每组 10 个 backbone；
- target B6–134 固定，binder 固定为 76 aa；
- 仍要求 `N_contact_hotspots >= 8` 和
  `N_hotspots_on_interface >= 4`；
- 每组最多取 3 个 backbone，每个生成 3 条 temperature `0.1`、
  omit-C 序列；
- 仍只用 `af2_model_1_multimer_tt_3rec` 做本阶段正向门控。

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_r3_partial_diffusion.sh"
```

至少三条正向通过者仍需完成 USP4/USP11 反筛和去冗余，才能成为最终
计算候选。

R3 已完成并判定未收敛：共 51 条目标模板 AF2 预测技术成功，0 条通过
三项联合门控。30/30 个 partial-diffusion backbone 虽然通过热点硬过滤，
但序列预测没有保持目标结合姿态；实验 6DJ9 阳性对照本身也被当前
`tt_3rec` 协议判为阴性。详细证据与 R4 所需决策见：

- [`docs/USP15_R3_RESULTS.md`](docs/USP15_R3_RESULTS.md)

在获得结构模板诊断或验证协议变更授权前，不继续放大同一分布，也不把
任何近门控结果标记为候选。

## 当前计算协议

本版本不使用 PyRosetta：

1. RFdiffusion 生成 RFD1 backbone。
2. OVO backbone metrics 执行硬过滤：
   - `N_contact_hotspots >= 8`
   - `N_hotspots_on_interface >= 4`
3. LigandMPNN 使用 ProteinMPNN weights，每个 R1 backbone 生成 3 条序列：
   - R1 temperature `0.1`
   - R2 诊断 temperature `0.05` 和 `0.10` 分别运行
   - omit `C`
   - 不设置氨基酸 bias
4. AF2 `model_1_multimer` target-template、3 recycles 执行 smoke 和主正筛。
5. 最终候选增加 `model_1_ptm` target-template 与 USP4/USP11 反筛。
6. 用界面 ΔSASA、碰撞、ProteinQC、ESM-IF、ProteinSol 和序列性质替代 Rosetta 指标。

OpenMM 相互作用能和埋藏极性原子检查仅用于排序，不能解释为 Rosetta ddG 或 Rosetta buried-unsatisfied hydrogen bonds。

## 运行环境

已测试环境：

- OVO 1.0.2
- Nextflow 25.10.4
- Docker + NVIDIA Container Toolkit
- 单张 NVIDIA V100 32 GB
- OVO Docker images：
  - `ovo-rfdiffusion`
  - `ovo-python-structure`
  - `ovo-ligandmpnn`
  - `ovo-colabdesign`

设置以下环境变量；不要将真实服务器路径写入仓库：

```bash
export USP15_CAMPAIGN_DIR=/path/to/USP15_DUSP_R1
export OVO_HOME_DIR=/path/to/initialized/ovo
export OVO_ENV_DIR=/path/to/conda/envs/ovo
```

准备 campaign 目录：

```bash
mkdir -p "$USP15_CAMPAIGN_DIR"/{config,inputs,scripts,reports}
cp config/* "$USP15_CAMPAIGN_DIR/config/"
cp scripts/* "$USP15_CAMPAIGN_DIR/scripts/"
chmod 750 "$USP15_CAMPAIGN_DIR"/scripts/*
```

## 准备目标结构

先从 RCSB 下载 `3T9L.pdb`，然后提取和验证 chain A 6–134：

```bash
python scripts/prepare_usp15_target.py \
  --source /path/to/3T9L.pdb \
  --output "$USP15_CAMPAIGN_DIR/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
  --report "$USP15_CAMPAIGN_DIR/reports/USP15_DUSP_3T9L_A6-134.validation.json"
```

脚本要求残基 6–134 连续、无插入码、主链原子完整且六个热点全部存在。

## 执行顺序

### 1. 四组 15-step preview

```bash
for pool in \
  USP15_R1_short_base \
  USP15_R1_short_beta \
  USP15_R1_long_base \
  USP15_R1_long_beta
do
  "$USP15_CAMPAIGN_DIR/scripts/run_preview.sh" "$pool"
done
```

如需保留新的随机 preview，传入第二个 attempt 参数：

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_preview.sh" USP15_R1_short_beta 2
```

### 2. 四组 100-backbone pilot

先启动第一组：

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_backbone_pilot.sh" USP15_R1_short_base
```

其余三组串行执行：

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_pilot_queue.sh"
```

服务器只有一张 GPU 时，不要并行启动多个 RFdiffusion 或 AF2 pool。

### 3. LigandMPNN + AF2 smoke

四组 pilot 均完成且至少有一个 backbone 通过硬过滤后运行：

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_sequence_af2_smoke.sh" USP15_R1_short_base
```

也可使用门控脚本等待 pilot 队列后自动选择第一个有合格 backbone 的 pool：

```bash
"$USP15_CAMPAIGN_DIR/scripts/run_smoke_after_pilots.sh"
```

### 4. 扩量门控

只有以下条件全部满足后才能扩展到每池 1000 个 backbone：

- 四组 preview 成功完成；
- 四组各 100-backbone pilot 成功完成；
- backbone 硬过滤实际产生合格设计；
- LigandMPNN 生成恰好 3 条无 Cys 序列；
- AF2 smoke 通过：
  - iPAE ≤ 10
  - target-aligned binder RMSD ≤ 2 Å
  - binder pLDDT ≥ 80

阈值失败时停止并报告，不自动放宽。

## 主要脚本

| 脚本 | 用途 |
|---|---|
| `prepare_usp15_target.py` | 清理并验证 3T9L target |
| `run_preview.sh` | 运行单个 15-step RFdiffusion preview |
| `validate_backbones.py` | 独立复核链长与热点接触 |
| `run_backbone_pilot.sh` | 生成 100 个完整 backbone 并运行 OVO 硬过滤 |
| `summarize_backbone_metrics.py` | 汇总 pilot 合格率 |
| `run_pilot_queue.sh` | 将剩余 pilot 串行排队 |
| `run_sequence_af2_smoke.sh` | 运行 LigandMPNN 3-sequence + AF2 smoke |
| `run_smoke_after_pilots.sh` | pilot 完成后的自动门控 |
| `audit_r2_phase_a_alignment.py` | 独立复算 Phase A 链映射与 RMSD |
| `prepare_r2_phase_b_resources.sh` | 无 PyRosetta 准备 scaffold/target SS 与邻接矩阵 |
| `install_rfdiffusion_scaffold_checkpoint.sh` | 下载并校验官方 scaffold checkpoint |
| `run_r2_phase_b_backbones.sh` | 运行一个 R2 backbone 条件及固定过滤 |
| `run_r2_phase_b_queue.sh` | 串行运行六个 50-backbone 条件 |
| `normalize_scaffold_trb.py` | 保留原始 TRB 并生成 OVO 兼容副本 |
| `build_r2_phase_c_matrix.py` | 从每组 top-5 构建序列/AF2 矩阵 |
| `run_r2_phase_c.sh` | 串行运行 R2 LigandMPNN 与 AF2 |
| `summarize_r2_phase_c.py` | 汇总固定 AF2 门控与失败原因 |
| `prepare_r3_ubv_controls.py` | 构建 6DJ9 UbV/3T9L 对照与诊断输入 |
| `prepare_r3_ubiquitin_scaffold.py` | 将稳定 1UBQ scaffold 放置到已知界面 |
| `run_r3_partial_diffusion.sh` | 串行执行 RFD1 partial diffusion、硬过滤、LigandMPNN 与 AF2 |
| `summarize_r3_partial_diffusion.py` | 汇总 R3 不变 AF2 门控 |

## 输出与可恢复性

每个 Nextflow 阶段独立保存：

- `trace.txt`
- `report.html`
- `nextflow.log`
- `nextflow.stdout.log`
- published PDB 和 CSV

脚本检测成功 trace 后会跳过已完成阶段，并使用 Nextflow `-resume` 恢复可缓存任务。

## 安全与科学边界

- 不提交 OVO token、SSH 地址、用户名、私有服务器路径或数据库。
- 登录凭据只在交互式会话中使用。
- `raw/`、`results/`、Nextflow `work/` 和模型文件不进入 Git。
- 计算结果只能作为后续实验验证的候选来源。
- 6DJ9 证明该 DUSP 表面可被 UbV 结合，但不保证设计会调控某个特定天然互作。
