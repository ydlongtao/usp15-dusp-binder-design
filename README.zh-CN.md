# USP15 DUSP Binder 设计项目

[详细项目记录](README.md) ·
[在线 HTML 报告](https://ydlongtao.github.io/usp15-dusp-binder-design/) ·
[仓库内报告文件](docs/USP15_R10_complete_report.html) ·
[10 个候选结构图](docs/figures/USP15_R10) ·
[MIT License](LICENSE)

这是一个针对人 USP15 N 端 DUSP 结构域的开源计算型微型蛋白 binder
设计项目。仓库保存从 RFdiffusion 骨架生成、LigandMPNN 序列设计、
AlphaFold2 结构评估、USP4/USP11 计算反筛，到最终候选导出的参数、
脚本、失败记录和结果报告。

## 当前结果

R10 固定面板的筛选过程为：

```text
52 个输入设计
  → 41 个通过 USP15 正向筛选
  → 36 个通过非 PyRosetta 界面审计
  → 28 个通过 USP4/USP11 同姿势计算反筛
  → 24 个序列簇
  → 10 个最终代表候选
```

最终候选均为 76 aa、无 Cys 的计算型微型蛋白。每个候选都提供：

- FASTA 序列和候选指标；
- 最佳 USP15 预测复合物 PDB；
- 全复合物和热点界面的统一视角结构图；
- USP4/USP11 计算反筛结果；
- ProteinQC、序列组成和溶解度相关指标；
- 后续 OpenMM 分子动力学与 SPR/MST 实验验证起始方案。

主要结果入口：

- [R10 结果说明](docs/USP15_R10_RESULTS.md)
- [自包含 HTML 完整报告](docs/USP15_R10_complete_report.html)
- [10 张高分辨率结构图](docs/figures/USP15_R10)
- [10 个配套复合物 PDB](docs/structures/USP15_R10)
- [MD 与 SPR/MST 验证计划](docs/USP15_R10_STRUCTURE_MD_AND_SPR_MST_PLAN.md)
- [机器可读 OpenMM 参数](config/usp15_r10_openmm_md.json)
- [10 个 MD 制备体系审计](docs/USP15_R10_MD_PREPARED_SYSTEM_AUDIT.json)
- [V100 CUDA smoke 审计](docs/USP15_R10_MD_CUDA_SMOKE_AUDIT.json)

## 靶点与设计约束

- 靶点：human USP15，UniProt `Q9Y4E8`
- 结构：3T9L chain A residues 6–134
- 界面参考：6DJ9 DUSP–UbV 复合物
- 热点：`A50,A52,A53,A55,A57,A61`
- 反筛：human USP4 和 USP11 DUSP
- RFdiffusion：仅使用 RFD1
- 不使用 PyRosetta
- backbone 硬门控：
  - `N_contact_hotspots >= 8`
  - `N_hotspots_on_interface >= 4`
- AF2 正向门控：
  - iPAE ≤ 10
  - target-aligned binder RMSD ≤ 2 Å
  - binder pLDDT ≥ 80

阈值失败时停止，不自动放宽门控，也不把“接近阈值”的设计标记为候选。

## 结果应如何解释

R10 使用经过阳性对照校准的 AF2 `model_2_ptm`
interface-template 协议，是 geometry-conditioned compatibility test，
不是独立的 sequence-only fold-and-dock 证明。

因此：

- 结构图是 AF2 预测，不是实验解析结构；
- 结构图不是分子动力学轨迹帧；
- OpenMM 执行队列已建立，10/10 体系和真实 CUDA smoke 已通过审计，
  首条 100 ns 生产轨迹正在运行；
- 计算反筛不等于已经获得实验选择性；
- 当前结果不能直接证明 USP15 抑制、天然互作阻断或细胞活性。

在宣称结合或选择性之前，需要使用 SPR、MST、BLI 或其他正交方法进行
实验验证，并在相同条件下测试 USP4 和 USP11。

## 仓库结构

```text
config/       设计矩阵、固定门控和机器可读 MD 参数
docker/       特定运行时兼容镜像定义
docs/         阶段计划、结果、HTML 报告、结构图和候选 PDB
scripts/      数据准备、运行、审计、汇总和报告生成脚本
README.md     完整计算历程和服务器端执行说明
Agent.md      自动化代理必须遵守的科学与运行约束
```

大型模型权重、Nextflow `work/`、原始服务器结果、数据库和登录凭据不会
进入仓库。

## 基本运行要求

经过测试的主要环境：

- Linux GPU 服务器
- Docker 与 NVIDIA Container Toolkit
- Nextflow 25.10.4
- NVIDIA V100 32 GB
- OVO 1.0.2 相关镜像

服务器只有一张 GPU 时，RFdiffusion、LigandMPNN、AF2 和 Boltz 等
GPU 阶段必须串行运行。运行路径通过环境变量传入，不要把真实服务器
路径、IP、用户名、token 或密码写入脚本。

```bash
export USP15_CAMPAIGN_DIR=/path/to/USP15_DUSP_R1
export OVO_HOME_DIR=/path/to/initialized/ovo
export OVO_ENV_DIR=/path/to/conda/envs/ovo
```

具体的 R1–R10 执行记录、失败恢复方式和脚本命令见
[详细 README](README.md)。

## 分子动力学与实验验证

已实现的 MD 队列使用 OpenMM 8.5.2，并通过 AmberTools `tleap`
构建 AMBER ff19SB/OPC 显式溶剂截角八面体。每个复合物运行 3 个独立的
100 ns 重复；分析 binder RMSD、界面接触占有率、热点保持、界面氢键、
buried SASA、相对 MM/GBSA 和重复间一致性。

10/10 个制备体系已经通过拓扑、SHA-256、OPC 四点水、盒矢量和
“生产期无位置约束”审计。真实 V100 CUDA smoke 已进一步确认 25,000
steps、5 个 ×10 ps 帧及 3,336 个蛋白轨迹原子；首条 100 ns 轨迹正在运行，
因此目前仍没有完整重复可用于稳定性或 MM/GBSA 结论。

SPR/MST 条件只作为方法开发起点。首轮实验应优先：

1. 检查候选的可溶表达、单体比例和聚集；
2. 测量 USP15 DUSP 的完整浓度依赖曲线；
3. 在匹配条件下平行测试 USP4 和 USP11 DUSP；
4. 对通过者复核全长 USP15；
5. 最后再开展机制或细胞实验。

## 安全与贡献

提交代码前请运行语法检查、JSON 校验和凭据扫描。不要提交服务器地址、
登录凭据、私有路径、模型权重或未经审计的大型结果文件。所有新增候选
都必须保留输入结构、残基编号、模型版本、参数、指标和淘汰原因。

## 许可证

本仓库原创代码和文档采用 [MIT License](LICENSE)。

第三方软件、模型权重和来源结构继续适用各自许可证或使用条款，包括但
不限于 RFdiffusion、LigandMPNN、AlphaFold/ColabDesign、Boltz、
OpenMM、OVO 和 RCSB PDB 条目。本项目许可证不会替代第三方条款。
