# USP15 R10 OpenMM 迁移至 Minerva

## 原服务器暂停点

原服务器上的 `usp15-r10-md` 队列已按用户要求停止，GPU 已释放。10/10
ff19SB/OPC 制备体系和 V100 CUDA smoke 审计均已通过。`rank01/seed0`
完成了最小化、NVT 和第一段 NPT，第二段 NPT 被中断；正式无约束生产
采样尚未开始。该部分状态已迁移到 Minerva，原服务器保持暂停且未被
重新连接或重启。

迁移必须保留：

- `inputs/`、`prepared/`、`scripts/`、`reports/`；
- 已通过的 `smoke/` 及其 `audit.json`；
- `runs/` 中的已完成阶段状态；
- `logs/`、`failed_attempts/` 和失败 smoke，作为审计记录；
- `build/usp15-openmm-8.5.2.docker.tar.gz` 及其 SHA-256。

## Minerva 调度与容器约定

Minerva 使用 IBM LSF，不使用 Slurm。正式 GPU 队列为 `gpu`，上限 144
小时；短 smoke 可用 `gpuexpress`。GPU 通过 `-gpu num=1` 申请，并用
`-R span[hosts=1]` 把 CPU/GPU 资源限制在一个节点。LSF 会设置
`CUDA_VISIBLE_DEVICES`，脚本不得覆盖它。

Minerva 支持 Apptainer/Singularity。计算节点默认无外网，因此先上传
OCI/Docker archive，再转换为 SIF。转换脚本会根据 `oci-layout` 或
`manifest.json` 自动选择 `oci-archive:` 或 `docker-archive:`。若归档
外层为 gzip，脚本会先在 `APPTAINER_TMPDIR` 中临时展开原始 tar，完成后
自动清理：

```bash
bash scripts/minerva/convert_docker_archive_to_sif.sh \
  /path/to/md_openmm/build/usp15-openmm-8.5.2.docker.tar.gz \
  /path/to/containers/usp15-openmm-8.5.2.sif
```

2026-07-27 的实际转换使用 Apptainer 1.4.5，生成的 SIF 为
1,210,654,720 bytes，SHA-256 为
`e270764ce7307c3f0b4bcfe339deb425d2093e9e436ea8f7f4240f9852ba8d92`。
容器内 OpenMM 版本核验为 8.5.2，临时展开归档和 rootfs 均在构建完成后
自动清理。

## 上传后校验

在本地生成 SHA-256 清单：

```bash
python scripts/create_minerva_transfer_manifest.py \
  --md-dir /local/path/md_openmm
```

上传后在 Minerva 验证：

```bash
python scripts/create_minerva_transfer_manifest.py \
  --md-dir /minerva/path/md_openmm \
  --verify
```

使用 SSH/rsync 上传时：

```bash
export LOCAL_MD_DIR=/local/path/md_openmm
export MINERVA_SSH=your_minerva_ssh_alias
export MINERVA_DEST=/minerva/absolute/project/path
bash scripts/minerva/upload_rsync.sh
```

## Minerva 专用 smoke

必须在 Minerva GPU 节点重新运行 smoke，不能仅依赖原 V100 的结果：

```bash
export LSF_PROJECT=acc_yourlab
export USP15_MD_DIR=/minerva/path/md_openmm
export OPENMM_SIF=/minerva/path/containers/usp15-openmm-8.5.2.sif
export MINERVA_GPU_MODEL=h100nvl
bash scripts/minerva/submit_smoke.sh
```

只有 `smoke_minerva/audit.json` 为 `passed` 才能提交正式数组。

首次 Minerva smoke（LSF `258355644`）在最小化后创建第二个 CUDA
Context 时失败。独立诊断确认 H100、驱动、CUDA Driver API、SIF 和单个
完整体系 Context 均正常；问题是前一个 OpenMM Context 仍存活时创建
下一个 Context。脚本现会在每个阶段结束后显式释放 Context 并运行垃圾
回收，未改变力场、温压、步长、约束或采样长度。首次失败目录和日志均
保留。

修复后的 smoke（LSF `258356236`）通过全部审计：

- OpenMM 8.5.2，CUDA 平台；
- 25,000 个生产步，实际 0.05 ns；
- 5 个间隔 10 ps 的轨迹帧；
- 3336 个蛋白原子；
- `production_restraints=false`。

smoke 只证明运行链路正确，其 0.05-ns 结构指标不能用于宣称 MD 稳定性
或结合结论。

## 30 个串行重复

```bash
export LSF_PROJECT=acc_yourlab
export USP15_MD_DIR=/minerva/path/md_openmm
export OPENMM_SIF=/minerva/path/containers/usp15-openmm-8.5.2.sif
export MINERVA_GPU_MODEL=h100nvl
bash scripts/minerva/submit_replica_array.sh
```

数组固定为 30 项、最大并发 1。每项对应一个 `rank × seed`，运行一个
100-ns NPT 重复，然后执行 20-ns burn-in 后的 RMSD、接触占有率、buried
SASA、氢键及相对 MM/GBSA。需要局部重提时可设置
`MINERVA_ARRAY_RANGE`，例如 `MINERVA_ARRAY_RANGE=1-3`。

2026-07-27 首次正式数组 LSF `258356480` 的第 1 项完成了 100-ns 轨迹和
结构分析，但相对 MM/GBSA 因容器内 `AMBERHOME` 未设置而失败；完整轨迹、
分析结果和失败日志均保留。已将 `APPTAINERENV_AMBERHOME=/opt/conda`
加入 Minerva 作业环境，并用容器内 `MMPBSA.py -h` 验证修复。旧数组的
未运行元素已取消，新的严格串行数组为 LSF `258390152`，范围
`1-30%1`，资源为 `select[h100nvl]` 和单 GPU exclusive-process。一次终端
输出延迟导致的重复提交 `258390152` 在启动前已取消，唯一保留的数组为
`258390007`。新数组
会复用第 1 项已完成的轨迹，仅重跑修复后的 MM/GBSA；其余项目从各自
的 checkpoint 或初始状态开始。

## 已验证的运行配置

- 已复用用户主动建立的 Minerva SSH 会话，不需要切换网络；
- 原服务器仍位于虚拟网络中，迁移执行期间不与其建立连接；
- LSF project account、工作目录和 scratch 目录均通过环境变量传入，
  不写入仓库；
- 正式 GPU 队列使用 `gpu`，smoke/诊断使用 `gpuexpress`；
- H100 NVL 资源语法为 `select[h100nvl]`；
- 全量迁移采用 `rsync`，上传后 660 个文件的 SHA-256 清单验证通过。

## 参考

- [Minerva GPGPU 与 GPU 资源申请](https://labs.icahn.mssm.edu/minervalab/documentation/gpgpu/)
- [Minerva LSF Job Scheduler](https://labs.icahn.mssm.edu/minervalab/documentation/lsf-job-scheduler/)
- [Minerva Apptainer/Singularity](https://labs.icahn.mssm.edu/minervalab/documentation/running-container-apptainer-singularity/)
- [Apptainer 从 Docker archive 构建 SIF](https://apptainer.org/docs/user/latest/docker_and_oci.html#containers-in-docker-archive-files)
