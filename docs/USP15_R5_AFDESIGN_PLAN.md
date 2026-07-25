# USP15 R5 AFDesign 定向序列优化计划

## 依据

R1–R4 共审计 200 条 OVO target-template AF2 记录，iPAE ≤ 10 和
target-aligned binder RMSD ≤ 2 Å 的记录均为 0。单纯扩大 RFdiffusion 或
ProteinMPNN 的既有分布缺少成功信号。R3 的
`P10_rank1_rfdiffusion_8_standardized_packed_1_1` 是当前唯一同时接近
iPAE 门控（11.38）并达到 binder pLDDT 门控（85.30）的设计，因此 R5
只在该 RFD1 正确界面骨架上做 AF2 梯度序列优化。

ColabDesign/AFDesign 是 OVO 服务器已预装的组件，不新增系统软件。它只
用于设计序列；最终正向验证仍由独立保留的
`af2_model_1_multimer_tt_3rec` 完成。

## 设计方法

- 输入：R3 P10 rank-1 的 76-aa RFD1 binder–USP15 DUSP 复合物。
- 靶点与热点保持不变：完整 3T9L A6–134；
  `50,52,53,55,57,61`。
- AFDesign 使用固定复合物坐标模板但移除 binder 序列模板，联合优化：
  - binder/复合物监督 distogram；
  - binder pLDDT 与内部 PAE；
  - 热点限制的界面 PAE 与界面接触。
- 设计 loss 以热点界面为主：`i_pae=5.0`、`i_con=5.0`；
  监督 distogram 仅以 `0.05` 保持折叠约束，避免其数值尺度压制界面梯度。
- 设计只使用 multimer model 3 和 model 4；model 1 完全保留给 OVO
  target-template 独立验证。
- 去除 Cys。
- GPU 阶段在单张 V100 上串行。

## 分阶段执行

1. 一个低迭代 smoke：10 个 soft iteration、2 个 semi-greedy
   iteration，验证链、热点、权重、Cys 排除和输出 PDB。
2. smoke 技术成功且 AFDesign 预测仍覆盖至少 4 个热点后，运行 3 个
   完整种子：
   - 80 个 soft iteration；
   - 16 个 semi-greedy iteration；
   - 每个 mutation step 10 次尝试；
   - 设计 recycles=1，model 3/4。
3. 每个种子生成一个无 Cys 的离散序列，以原 RFD1 正确界面骨架作为
   OVO 验证输入。
4. 只运行 `af2_model_1_multimer_tt_3rec`，门控不变：
   - iPAE ≤ 10；
   - target-aligned binder RMSD ≤ 2 Å；
   - binder pLDDT ≥ 80。

## 决策

- 至少 3 个完整种子通过：进入 USP4/USP11 反筛和序列去冗余。
- 1–2 个通过：只增加不同设计种子，直到获得 3 个去冗余通过者或该
  优化分布的预设预算耗尽。
- 0 个通过：停止 R5，不按 AFDesign 内部 loss 或单项 AF2 指标晋级。
- 不改变靶点、热点、AF2 测试或阈值，不使用 PyRosetta，不把 AFDesign
  训练内指标当成候选门控。

## Model-in-the-loop 有界诊断

如果 model 3/4 smoke 无法形成热点界面，可额外运行一次 20-step 诊断：

- 使用 model 1 做设计期梯度，但不给 binder 坐标或序列模板；
- 仍只提供完整 USP15 target template 与六个热点；
- 不做 semi-greedy 离散放大；
- 仅判断热点 iPAE/i-contact 是否呈持续改善趋势。

该诊断与最终 model 1 测试不独立，必须在报告中标记
`model_in_the_loop=true`。内部指标不能计为候选；若后续序列进入原
OVO 验证，仍必须满足全部固定门控和 USP4/USP11 反筛。
