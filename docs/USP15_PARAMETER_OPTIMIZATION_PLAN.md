# USP15 DUSP 参数优化计划（R2）

记录日期：2026-07-25

状态：**阶段 A 已于 2026-07-25 完成且未通过 AF2 门控；阶段 B 尚未启动。**

本文件记录 `USP15_DUSP_R1` 未通过 AlphaFold2 smoke gate 后的参数优化方案，供后续恢复计算时直接复用。R2 保持 USP15 DUSP 靶点、输入结构和最终验收阈值不变，不通过放宽门控来制造“合格”结果。

## 1. 不变条件

- 靶点：human USP15 N-terminal DUSP domain。
- 输入结构：3T9L chain A residues 6–134。
- 完整评价热点：`A50,A52,A53,A55,A57,A61`。
- RFdiffusion：RFD1。
- 单 GPU 上的 GPU-heavy 阶段必须串行。
- 不使用 PyRosetta。
- 不自动放宽以下 backbone 门控：
  - `N_contact_hotspots >= 8`
  - `N_hotspots_on_interface >= 4`
- 不自动放宽以下 AF2 门控：
  - iPAE `<= 10`
  - target-aligned binder RMSD `<= 2 Å`
  - binder pLDDT `>= 80`
- LigandMPNN 继续 omit `C`，不设置氨基酸 bias。
- 未出现完整 AF2 通过者前，不启动 1000-backbone 扩量。

R2 中用于 RFdiffusion 条件化的热点可以是完整六热点的子集，但所有设计必须继续用完整六热点计算最终接口指标。这是“生成参数优化”，不是更换靶点或降低评价标准。

## 2. R1 结果和失败诊断

四个 100-backbone pilot 均完成：

| Pool | 通过 backbone 门控 | 通过率 |
|---|---:|---:|
| `USP15_R1_short_base` | 52/100 | 52% |
| `USP15_R1_short_beta` | 19/100 | 19% |
| `USP15_R1_long_base` | 69/100 | 69% |
| `USP15_R1_long_beta` | 29/100 | 29% |
| 合计 | 169/400 | 42.25% |

LigandMPNN 和 `af2_model_1_multimer_tt_3rec` 技术流程运行成功，每个测试 backbone 均生成了 3 条无 Cys 序列，但两个 Complex_base backbone 的 6 条序列均未通过 AF2。复核 smoke 的最佳结果仍为：

- iPAE：25.74
- target-aligned binder RMSD：24.16 Å
- binder pLDDT：62.12

对应要求分别是 `<=10`、`<=2 Å`、`>=80`。因此失败不是“接近阈值”，不能用增加同分布样本或放宽阈值解决。

两个已测试 Complex_base backbone 都高度 α 螺旋化，其中一个接近连续单螺旋。当前主要假设是：

1. backbone 接触指标合格，但缺乏紧凑、可独立折叠的微型蛋白拓扑；
2. 任意选择第一个合格 backbone 不是可靠的 smoke 代表；
3. Complex_beta 中更紧凑、含转角或混合拓扑的候选尚未被充分测试；
4. 需要把“生成热点”和“评价热点”解耦，以减少过度约束，同时保留完整六热点评价。

## 3. 阶段 A：复用现有 beta backbone

在生成新 backbone 前，先对已有 Complex_beta 候选做低成本诊断。优先选择：

| 候选 | 已有 backbone 指标 |
|---|---|
| `USP15_R1_short_beta/rfdiffusion_55_standardized` | 6/6 热点、43 个热点接触、Rg 14.70 Å |
| `USP15_R1_short_beta/rfdiffusion_75_standardized` | 5/6 热点、38 个热点接触、Rg 9.46 Å |
| `USP15_R1_long_beta/rfdiffusion_76_standardized` | 5/6 热点、38 个热点接触、Rg 15.45 Å |

每个 backbone 运行两组 LigandMPNN：

- temperature `0.05`，3 条序列；
- temperature `0.10`，3 条序列；
- omit `C`；
- 不设置氨基酸 bias。

共生成 18 条序列，并全部运行：

- `af2_model_1_multimer_tt_3rec`
- target-template
- 3 recycles

阶段 A 的判定：

- 任意一条序列完整通过三项 AF2 门控：保留成功的模型、长度和 temperature 条件，进入针对性扩量准备。
- 18 条全部失败：进入阶段 B。
- 不用“接近通过”替代正式通过。

temperature `0.05` 是 R2 的诊断变量，不覆盖 R1 的 temperature `0.1` 记录。实施前应在 R2 配置和 Agent 指令中显式登记。

### 3.1 阶段 A 执行记录

阶段 A 已按上述矩阵完成：

- 6/6 个 LigandMPNN 条件技术成功；
- 18/18 条预期序列已生成；
- 所有 binder 序列均无 Cys；
- 18/18 条序列的 `af2_model_1_multimer_tt_3rec` 技术成功；
- 0/18 条序列完整通过三项 AF2 门控；
- 未启动任何扩量任务。

最佳单项指标分别为：

| 指标 | 最佳值 | 对应要求 |
|---|---:|---:|
| iPAE | 25.06 | `<=10` |
| target-aligned binder RMSD | 24.97 Å | `<=2 Å` |
| binder pLDDT | 57.49 | `>=80` |

这些最佳值来自不同或部分不同设计，不能组合成一个候选。由于所有结果与正式门控仍有很大距离，阶段 A 的结论为 `af2_gate_failed`。下一步只能在复核链映射和 target-template 对齐后，另行授权阶段 B；不得从阶段 A 自动扩量。

## 4. 阶段 B：新 RFdiffusion 小矩阵

不直接生成数千个 backbone。先运行 6 个条件，每个 50 个 backbone，共 300 个。

### 4.1 生成热点

完整六热点始终用于最终评价。RFdiffusion 生成阶段测试以下两个四热点子集：

- `H4-hydrophobic = A50,A53,A55,A61`
  - Trp50、Tyr53、Met55、Tyr61，强调疏水/芳香核心；
- `H4-balanced = A50,A52,A55,A61`
  - 在疏水锚点中加入 Lys52，测试极性定位作用。

### 4.2 参数矩阵

| 条件 | Binder 长度 | RFdiffusion 模型 | 生成热点 | Backbone 数 |
|---|---:|---|---|---:|
| B1 | 45–55 | Complex_beta | H4-hydrophobic | 50 |
| B2 | 45–55 | Complex_beta | H4-balanced | 50 |
| B3 | 56–70 | Complex_beta | H4-hydrophobic | 50 |
| B4 | 56–70 | Complex_beta | H4-balanced | 50 |
| S1 | 50–65 | scaffold-guided Complex_base | H4-hydrophobic | 50 |
| S2 | 60–75 | scaffold-guided Complex_base | H4-balanced | 50 |

共同参数：

- `diffuser.T=50`
- RFD1
- target 为 3T9L A6–134
- 不使用 auxiliary potentials
- 不使用 PyRosetta

S1/S2 的目的不是增加普通 Complex_base 样本，而是用 fold/scaffold conditioning 约束紧凑的三螺旋束或其他可折叠拓扑，避免再次得到连续单螺旋。实施前必须先确认当前 RFD1 镜像中 scaffold-guided PPI 所需文件和参数可用；验证失败时不得静默改用其他模型。

## 5. R2 backbone 过滤

保留 R1 两项硬门控，并增加紧凑性和拓扑过滤。

### 5.1 必须保留的门控

- `N_contact_hotspots >= 8`
- `N_hotspots_on_interface >= 4`
- 上述热点统计必须基于完整六热点，而不是生成子集。

### 5.2 新增过滤

- 45–55 aa：radius of gyration `<=15.5 Å`
- 56–70 aa：radius of gyration `<=18 Å`
- 至少两个由 loop 分隔的二级结构片段
- 最长连续螺旋不超过 binder 长度的 65%
- loop 比例建议为 5–40%
- `N_contact_interface / binder_length >= 1.0`
- 无链断裂、严重原子碰撞或明显游离末端
- 至少 4/6 完整热点参与界面

Rg 和二级结构阈值属于 R2 新增硬过滤，不替代原热点门控。首次 50-backbone pilot 后应报告它们各自的淘汰数量，避免无法判断哪个条件过严。

## 6. 阶段 C：序列设计和 AF2 选择

每个 R2 条件执行：

1. 生成 50 个 backbone。
2. 运行原热点门控和新增紧凑性过滤。
3. 按热点覆盖、接触密度、紧凑性和拓扑多样性选择前 5 个 backbone。
4. 每个 backbone 在 temperature `0.05` 和 `0.10` 下各生成 3 条序列。
5. 每个参数条件最多运行 30 条 AF2 smoke。
6. 保存所有失败原因和完整指标，不只保存通过者。

主判定仍为：

- iPAE `<=10`
- target-aligned binder RMSD `<=2 Å`
- binder pLDDT `>=80`

只有至少一条序列同时通过三项门控，该 RFdiffusion 条件才有资格扩量。扩量优先选择出现正式通过者且拓扑不单一的条件，而不是仅按 backbone 接触通过率排序。

## 7. 扩量和停止规则

### 7.1 允许扩量

某一 R2 条件至少出现一条完整 AF2 通过者后：

- 固定对应的模型、长度、生成热点和序列 temperature；
- 为该条件准备 1000-backbone 运行；
- 单 V100 上保持 GPU 阶段串行；
- 扩量前再次记录可复现的完整命令和随机种子策略。

不默认把所有六个 R2 条件都扩展到 1000。

### 7.2 必须停止

以下任一情况发生时停止，不自动放宽：

- 阶段 A 的 18 条序列和阶段 B/C 的 180 条序列均无完整 AF2 通过者；
- scaffold-guided 所需模型或许可不可用；
- 需要改变最终热点评价、AF2 阈值或靶点结构才能继续；
- 出现重复 OOM、链映射错误或 target-template 无法正确应用。

停止后可讨论新的 RFD1 拓扑条件、partial diffusion 或 binder diversification，但必须建立新 Round 并保留 R1/R2 的失败记录。

## 8. 实施前检查清单

- [x] 确认 R1 没有活动的扩量任务。
- [ ] 手工复核 AF2 输入与输出的 binder A / target B 链映射。
- [ ] 确认 target-template 对齐正确，排除指标计算错误。
- [ ] 为 R2 建立独立 Round 和输出目录，不覆盖 R1。
- [x] 在 `config/` 中新增阶段 A 配置，不修改 R1 结果。
- [x] 更新 `Agent.md`，明确阶段 A 的 temperature 诊断矩阵。
- [x] 验证 Complex_beta 阶段 A smoke。
- [x] 记录 LigandMPNN 每 backbone 的实际序列数、temperature 和 Cys 数量。
- [x] AF2 smoke 完整通过前不创建或启动 1000-backbone 队列。

## 9. 参考资料

- [OVO RFdiffusion binder design workflow](https://ovo.dichlab.org/docs/rfdiffusion/binder_design.html)
- [RosettaCommons RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [LigandMPNN](https://github.com/dauparas/LigandMPNN)
