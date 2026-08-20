# ESM3 接入 OVO 的分阶段实施方案

## 目标

在不改变 USP15 DUSP 靶点、热点、AF2 门控或 USP4/USP11 反筛规则的前提下，评估并接入 ESM3 作为 OVO 的正交序列/结构评估模块。ESM3 初期不替代 RFdiffusion、LigandMPNN、AF2，也不作为结合亲和力或选择性的直接证据。

## 当前 OVO 工作流中的定位

1. RFdiffusion RFD1：生成界面 backbone。
2. LigandMPNN/ProteinMPNN：设计 binder 序列。
3. FastRelax：进行结构松弛。
4. AF2/ESMFold：复折叠和结构质量检查。
5. ProteinQC：序列、二级结构、溶解性和复杂度检查。
6. Nextflow/插件：承载可复现的自定义计算步骤。

OVO 文档中已有 ESMFold 复折叠选项，但 ESMFold 与 ESM3 不是同一个模型。ESM3 应通过独立适配器接入，而不是把 ESMFold 的输出名称直接替换为 ESM3。

## 分阶段方案

### 阶段 A：独立 ESM3 smoke test

目的：确认本地运行时、模型权重、GPU/CPU 资源和输出格式可用。

- 使用官方 `esm3-sm-open-v1` 权重。
- 对一条现有 USP15 binder FASTA 执行 sequence-only 生成。
- 对同一输入执行 structure-track 生成（若本地运行时和显存允许）。
- 记录 Python、PyTorch、ESM 包版本、设备、显存、运行时间、输出文件和失败原因。
- 不导入 OVO，不运行批量设计，不改变任何现有筛选阈值。

阶段 A 的通过条件：模型成功加载；sequence-track 成功；structure-track 在当前设备上成功，或有明确、可复现的资源/运行时阻断记录。

### 阶段 B：OVO Nextflow 评估步骤

建立独立的 `esm3_eval.nf` 和容器/脚本，输入候选 FASTA、backbone PDB 和 AF2 复合物，输出 JSON、FASTA、PDB 及 CSV 指标。ESM3 结果先作为可选 descriptor 展示和排序特征，不作为硬门控。

建议输出：

- ESM3 生成序列和结构；
- 模型配置、随机种子、设备和软件版本；
- 序列/结构轨迹可复现性指标；
- 与 AF2/ESMFold 结构的一致性指标；
- 热点接触保留和界面几何诊断。

### 阶段 C：ESM3 候选序列生成

在阶段 B 稳定后，允许 ESM3 针对固定的 USP15 DUSP backbone 或界面条件生成候选序列。所有候选仍必须经过原有 ProteinQC、AF2 正向门控和 USP4/USP11 反筛。ESM3 内部评分不能替代 iPAE、target-aligned binder RMSD、binder pLDDT、实验 KD 或 Rosetta/物理能量指标。

阶段 C 的第一步采用固定单链 binder backbone 的序列生成试运行。由于当前 ESM3 适配器不接受 USP15-binder 复合物作为联合条件，试运行结果只代表 backbone-conditioned sequence proposal；在通过 AF2 复合物重建和原有反筛前，不得解释为 USP15 结合候选。

执行脚本：`scripts/run_esm3_phase_c_generate.py`。该脚本固定随机种子、温度和采样步数，输出 FASTA、JSON 审计和 CSV 清单；默认不允许 Cys 的后处理判断，但不替代 ProteinQC。

### 阶段 D：校准与长期维护

使用已知阳性/阴性对照比较 ESM3、AF2 和 ESMFold 的一致性，建立版本化阈值和失败分类。只有完成独立校准后，才考虑把 ESM3 描述符用于候选排序；不因 ESM3 单一分数放宽原有门控。

## 许可和数据安全

优先使用本地权重，不上传 USP15 设计序列到第三方 API。ESM3 开放模型受非商业许可约束，仓库只提交脚本、配置、结果摘要和许可证说明，不提交模型权重、访问令牌、缓存或私有序列数据。

## 阶段 A 执行记录

执行脚本：`scripts/run_esm3_phase_a_smoke.py`。

机器可读结果：`results/esm3_phase_a/phase_a_summary.json`（运行后生成；该目录在 `.gitignore` 中，不上传原始模型输出）。
