# USP15 R10 结构图、MD 参数与 SPR/MST 验证方案

## 状态与解释边界

本文件为 R10 十个候选补充三类可交付物：

1. 十个 USP15 DUSP–binder 最佳正向 AF2 seed 的统一结构图；
2. 对应的可下载 PDB，作为后续分子动力学初始构象；
3. 一套预先定义、尚未执行的显式溶剂 MD 与 SPR/MST 起始方案。

结构图是 `model_2_ptm` interface-template 条件下的预测复合物，不是 MD
轨迹快照，也不是晶体、冷冻电镜或 NMR 实验结构。本文没有产生新的 MD
轨迹或结合自由能结果；后续只有在轨迹完成、质控和重复间一致性检查通过后，
才能把相应数据写成 MD 结果。

## 结构文件与链编号

- binder：预测文件 chain A，76 aa；
- USP15 DUSP：预测文件 chain B，129 aa，对应 3T9L chain A residues
  6–134；
- 3T9L 源热点 `A50/A52/A53/A55/A57/A61` 在预测文件中对应
  `B45/B47/B48/B50/B52/B56`；
- 每个候选使用 R10 manifest 指定的 best positive seed；
- 图片位于 [`figures/USP15_R10`](figures/USP15_R10)，PDB 位于
  [`structures/USP15_R10`](structures/USP15_R10)。

图片左侧显示整体复合物，右侧显示 DUSP 热点界面近景。USP15 DUSP、
binder 和六个热点分别以蓝灰、橙色和洋红色表示。十张图使用同一渲染参数，
便于横向比较，但不能通过视觉外观推断亲和力。

## 候选初始结构与理化信息

分子量、理论 pI 和 pH 7.4 净电荷由 Biopython `ProteinAnalysis` 根据
未加标签的 76-aa 序列计算，属于序列层面的近似值。实际构建中的标签、
末端修饰和缓冲液会改变这些数值。

| Rank | Best seed | 分子量 (Da) | 理论 pI | pH 7.4 预测净电荷 |
|---:|---:|---:|---:|---:|
| 1 | 2 | 8653.81 | 5.193 | -5.623 |
| 2 | 2 | 8563.77 | 4.904 | -4.733 |
| 3 | 2 | 8592.80 | 4.730 | -6.688 |
| 4 | 2 | 8441.53 | 4.910 | -4.693 |
| 5 | 2 | 8473.76 | 4.799 | -4.731 |
| 6 | 0 | 8638.97 | 5.888 | -1.705 |
| 7 | 2 | 8321.37 | 4.591 | -7.691 |
| 8 | 2 | 8571.62 | 4.630 | -7.729 |
| 9 | 2 | 8510.52 | 4.728 | -5.707 |
| 10 | 2 | 8434.64 | 5.205 | -3.319 |

## 建议的显式溶剂 MD 协议

机器可读的完整参数位于
[`config/usp15_r10_openmm_md.json`](../config/usp15_r10_openmm_md.json)。
建议使用 OpenMM 8.5 或更高版本，蛋白采用
`amber19/protein.ff19SB.xml`，水和离子采用 `amber19/opc.xml`。OpenMM
官方文档列出了 Amber19、ff19SB、OPC、PME 和约束等对应实现：
[OpenMM 8.5 Force Fields and Simulation Setup](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)。

### 体系准备

| 参数 | 固定值 |
|---|---|
| 输入 | 每个候选的 best positive seed PDB |
| 模拟构建 | 无标签 USP15 DUSP + 无标签 binder |
| 质子化 | pH 7.4；His 状态按局部氢键环境逐体系记录 |
| 末端 | 默认两性离子；只有实验样品确实封端时才改为封端 |
| 力场 | AMBER ff19SB |
| 水模型 | OPC |
| 周期盒 | dodecahedron，溶质到盒边至少 1.2 nm |
| 盐 | 0.15 M NaCl，并中和体系净电荷 |
| 长程静电 | PME |
| 非键截断 | 1.0 nm |
| 约束 | 含氢键 `HBonds`，刚性水 |

所有缺失的标准重原子和氢原子必须重建。不得为了让轨迹稳定而手工添加
binder–target 距离约束；位置约束只用于预平衡并在生产期完全移除。

### 最小化、平衡与生产

| 阶段 | 参数 |
|---|---|
| 能量最小化 | tolerance 10 kJ mol⁻¹ nm⁻¹；最多 20,000 iterations |
| NVT | 0.5 ns；蛋白重原子位置约束 1000 kJ mol⁻¹ nm⁻² |
| NPT-1 | 0.5 ns；位置约束 100 kJ mol⁻¹ nm⁻² |
| NPT-2 | 1.0 ns；位置约束 10 kJ mol⁻¹ nm⁻² |
| 生产系综 | NPT，300 K，1 bar |
| 积分器 | `LangevinMiddleIntegrator` |
| 时间步长 | 2 fs，不使用氢质量重分配 |
| 摩擦系数 | 1 ps⁻¹ |
| 压强耦合 | `MonteCarloBarostat`，每 25 steps 尝试 |
| 初筛采样 | 每个复合物 3 个独立重复 × 100 ns，共 3.0 μs |
| 扩展采样 | 入选候选每个重复延长至总计 500 ns |
| 轨迹输出 | 坐标和状态每 10 ps；checkpoint 每 1 ns |
| 随机性 | 每个体系、每个重复使用唯一且被记录的速度 seed |

100 ns × 3 是实验排序前的计算初筛，不足以证明热力学收敛。若 Rank 1–3
进入实验验证，应为其补做游离 binder 以及 apo USP15 DUSP 对照，并把复合物
重复延长到至少 500 ns；不同重复必须分别检查，不能先拼接再判断稳定性。

### 轨迹分析与预先定义的计算分流条件

生产轨迹前 20 ns 作为 burn-in。所有指标先按重复计算，再报告重复间中位数
与 95% bootstrap CI：

- 以 USP15 DUSP chain B Cα 对齐后的 binder Cα RMSD；
- binder 与 DUSP 的 backbone RMSF 和二级结构保持率；
- binder–target 质心距离、4.5 Å 界面接触数；
- 六个热点的逐残基接触占有率；
- 界面氢键占有率：供体–受体距离 ≤3.5 Å、角度 ≥135°；
- 盐桥占有率：相反电荷重原子距离 ≤4.0 Å；
- 界面 buried SASA。

为了避免分析完成后再挑阈值，建议预先使用以下仅用于实验顺序排序的条件：

- 任何重复均不出现持续性复合物分离；
- post-burn-in binder RMSD 中位数 ≤2.5 Å；
- 至少 70% 帧的 binder RMSD ≤3.0 Å；
- 至少 4/6 个热点的接触占有率 ≥50%；
- buried SASA 中位数 ≥600 Å²；
- 三个重复给出一致的是否保留界面的结论。

这些条件未经实验校准，不能称为结合阳性标准。MM/GBSA 如使用，只能作为
等间隔 post-burn-in 帧上的相对诊断，不能直接换算为 KD，也不能单独决定
候选淘汰。

## SPR 起始方案

十个 binder 只有约 8.3–8.7 kDa，建议把 USP15 DUSP 固定在传感器表面，
把 binder 作为 analyte。优先使用可定向捕获的生物素化或 His-tag USP15
DUSP；若使用胺偶联，先做低密度固定并验证 target 活性。binder 本身无
Cys，不建议为了首轮 SPR 临时引入 Cys 并直接把修饰体当作原序列。

| 项目 | 建议起始条件 |
|---|---|
| 运行缓冲液 | 10 mM HEPES pH 7.4、150 mM NaCl、0.05% Tween-20 |
| 还原剂 | 仅在 USP15 DUSP 稳定性需要时加入 0.5–1 mM TCEP，并在全部样品中一致 |
| analyte 梯度 | 0.5 nM–10 μM，先做 12–16 点两倍梯度 |
| association | 120 s |
| dissociation | 300–600 s |
| target 密度 | 先从约 100–500 RU 低密度表面开始 |
| 质控 | 空白流道、buffer blank、double reference、至少 2 次独立重复 |
| 拟合 | 先检查稳态 KD；只有无明显传质/再结合且残差支持时才报告 1:1 kon/koff |

若高浓度出现非特异表面结合，应优先降低固定密度、提高盐或加入少量 BSA，
而不是删去高浓度点后直接拟合。USP4 和 USP11 DUSP 应使用相同缓冲液、
固定策略和浓度范围平行反筛。

## MST 起始方案

建议标记 USP15 DUSP，而不是改造无 Cys binder。可使用 His-tag
RED-tris-NTA 标记或经验证不影响界面的其他 target 标记方式。

| 项目 | 建议起始条件 |
|---|---|
| labeled USP15 DUSP | 10–50 nM，以荧光和 capillary scan 结果定值 |
| binder 梯度 | 0.3 nM–10 μM，16 点两倍梯度 |
| 缓冲液 | 10 mM HEPES pH 7.4、150 mM NaCl、0.01–0.05% Tween-20 |
| 抗吸附 | 先做 capillary scan；必要时加入 0.05% BSA |
| 重复 | 每个候选至少技术重复，入选者进行 3 次独立实验 |
| 反筛 | USP4、USP11 使用相同标记策略和浓度范围 |

需要同时检查原始荧光、温度跃迁、热泳曲线和毛细管吸附。只有信号随浓度
呈可重复饱和趋势且不同标记/测定方式不改变排序时，才把拟合 KD 用于候选
比较。首轮优先顺序仍建议 Rank 1、2、3；若资源允许扩展至五个，再加入
Rank 9 和 Rank 10。

## 报告 MD 结果时必须补齐的字段

后续每个体系应记录 OpenMM 版本、CUDA/驱动、GPU、XML 文件校验和、完整
质子化状态、盒矢量、原子/水/离子数、随机 seed、实际步数、失败重启记录、
轨迹校验和和分析脚本版本。没有这些字段的轨迹不进入候选比较。
