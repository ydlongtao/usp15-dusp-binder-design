# USP15 R4–R5 结果与当前阻塞

## 结论

R4 和 R5 均技术完成，但在不变的 OVO
`af2_model_1_multimer_tt_3rec` 正向门控下没有通过者：

- R4：12/12 记录完成，0 通过；
- R5：4/4 记录完成，0 通过；
- 因此没有启动 USP4/USP11 反筛，也没有新增候选。

这些结果不能解释为目标表面不可结合。6DJ9 已经给出了实验结构明确的
UbV 结合界面；更直接的结论是，当前 target-template 验证协议没有重现
这个界面，包括已知 6DJ9 阳性对照。

## R4：6DJ9 四构象集

### 输入审计

6DJ9 不对称单元的 A–K、B–L、C–J 和 D–H 四组 DUSP–UbV 构象均移植到
完整 3T9L A6–134 靶点。三维最小排斥平移后：

| 构象 | DUSP CA 拟合 RMSD (Å) | 去碰撞位移 (Å) | `<2 Å` 原子对 | 5 Å 热点 |
|---|---:|---:|---:|---|
| A–K | 0.993 | 0.930 | 0 | 50,52,53,55,57 |
| B–L | 1.079 | 0.780 | 0 | 50,52,53,55,57 |
| C–J | 1.123 | 1.029 | 0 | 50,52,53,55,57 |
| D–H | 1.485 | 1.197 | 0 | 50,52,53,55,57 |

独立 CA 硬门控全部通过：

- `N_contact_hotspots`：26、31、29、31；
- `N_hotspots_on_interface`：四组均为 5。

ProteinMPNN 每构象生成 3 条序列。12/12 均为 76 aa、无 Cys、保留
15 个晶体界面位置，且全局互不重复。

### OVO 结果

12 条序列均未保持输入界面：

- iPAE：23.89–26.81；
- target-aligned binder RMSD：47.68–54.47 Å；
- binder pLDDT：42.32–59.66。

R4 按预设停止规则终止，没有放大或反筛。

## R5：AFDesign 定向序列优化

R5 选用历史最接近 iPAE/pLDDT 门控的 RFD1 P10 rank-1 复合物。
服务器已有 `ovo-colabdesign` 镜像和 AF2 权重，没有修改 OVO Python
环境，也没有使用 PyRosetta。

依次审计了：

1. model 3/4、不使用 binder template；
2. model 3/4、使用固定 binder 结构模板；
3. model 3/4、提高热点 `i_pae`/`i_con` loss 权重；
4. 一个 20-step、无 binder template 的 model-1-in-the-loop 诊断。

四个 smoke 均生成 76-aa、无 Cys且互不重复的离散序列。AFDesign
内部预测没有形成稳定热点界面；model-1-in-the-loop 的最佳瞬时
normalized iPAE 为约 0.81，但第 20 步回到 0.89，内部 RMSD 为约
56 Å。

原 OVO 验证结果：

| Seed | iPAE | Binder RMSD (Å) | Binder pLDDT | 通过 |
|---:|---:|---:|---:|---|
| 9000 | 25.04 | 48.74 | 76.69 | 否 |
| 9001 | 11.57 | 50.01 | 84.30 | 否 |
| 9002 | 26.13 | 54.35 | 78.38 | 否 |
| 9100 | 25.92 | 60.88 | 76.47 | 否 |

Seed 9001 只接近 iPAE 和 pLDDT，binder 仍迁移到错误表面，不能计为候选。

## 累积证据

R1–R4 的横向审计覆盖 200 条已汇总的 OVO AF2 记录：

- iPAE ≤ 10：0；
- target-aligned binder RMSD ≤ 2 Å：0；
- binder pLDDT ≥ 80：13；
- 同时满足两项门控：0。

加上 R5 的 4 条记录后仍没有正向通过者。尤其重要的是：

- 6DJ9 晶体 UbV 阳性对照未通过；
- 完整无 Cys UbV/ubiquitin 支架未通过；
- RFD1 设计、partial diffusion、LigandMPNN、ProteinMPNN、
  四个实验构象和 AFDesign 均未通过；
- 失败共同表现为 binder 离开设计界面，而不是仅差一个数值阈值。

## 当前阻塞与允许的下一步

在“不更换靶点、不改变 `tt_3rec` 验证模型、不放宽三项门控”的组合约束
下，继续扩大相同分布缺少科学依据，也不能诚实地产生 3 个候选。

解除阻塞需要用户授权至少一种验证层面的改变，例如：

- 增加 binder/interface-template AF2 诊断，先证明验证器能重现 6DJ9；
- 允许一个独立结构预测器（例如已安装的 Boltz）作为正交验证；
- 重新定义阳性控制失败时的验证决策树。

这些都不等同于降低门控；目的应是先建立能够识别已知阳性结合构象的
验证协议。获得明确授权前，不把任何近门控结果晋级为候选。
