# USP15 R10 OpenMM 迁移至 Minerva

## 当前暂停点

原服务器上的 `usp15-r10-md` 队列已按用户要求停止，GPU 已释放。10/10
ff19SB/OPC 制备体系和 V100 CUDA smoke 审计均已通过。`rank01/seed0`
只完成了最小化和 NVT，正在进行的 NPT-1 被中断；正式无约束生产采样尚未
开始，因此当前仍为 0/30 个完整重复。

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
Docker archive，再在登录节点转换为 SIF：

```bash
bash scripts/minerva/convert_docker_archive_to_sif.sh \
  /path/to/md_openmm/build/usp15-openmm-8.5.2.docker.tar.gz \
  /path/to/containers/usp15-openmm-8.5.2.sif
```

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

## Minerva 专用 smoke

必须在 Minerva GPU 节点重新运行 smoke，不能仅依赖原 V100 的结果：

```bash
export LSF_PROJECT=acc_yourlab
export USP15_MD_DIR=/minerva/path/md_openmm
export OPENMM_SIF=/minerva/path/containers/usp15-openmm-8.5.2.sif
export MINERVA_GPU_MODEL=a100
bash scripts/minerva/submit_smoke.sh
```

只有 `smoke_minerva/audit.json` 为 `passed` 才能提交正式数组。

## 30 个串行重复

```bash
export LSF_PROJECT=acc_yourlab
export USP15_MD_DIR=/minerva/path/md_openmm
export OPENMM_SIF=/minerva/path/containers/usp15-openmm-8.5.2.sif
export MINERVA_GPU_MODEL=a100
bash scripts/minerva/submit_replica_array.sh
```

数组固定为 30 项、最大并发 1。每项对应一个 `rank × seed`，运行一个
100-ns NPT 重复，然后执行 20-ns burn-in 后的 RMSD、接触占有率、buried
SASA、氢键及相对 MM/GBSA。需要局部重提时可设置
`MINERVA_ARRAY_RANGE`，例如 `MINERVA_ARRAY_RANGE=1-3`。

## 尚需用户提供

- Minerva 用户名或已配置的 SSH host alias；
- LSF project account，例如 `acc_xxx`；
- Minerva 目标存储路径；
- 希望使用的 GPU 型号；默认模板为 `a100`；
- 是否通过 `rsync/scp` 或 Globus 上传。

## 参考

- [Minerva GPGPU 与 GPU 资源申请](https://labs.icahn.mssm.edu/minervalab/documentation/gpgpu/)
- [Minerva LSF Job Scheduler](https://labs.icahn.mssm.edu/minervalab/documentation/lsf-job-scheduler/)
- [Minerva Apptainer/Singularity](https://labs.icahn.mssm.edu/minervalab/documentation/running-container-apptainer-singularity/)
- [Apptainer 从 Docker archive 构建 SIF](https://apptainer.org/docs/user/latest/docker_and_oci.html#containers-in-docker-archive-files)
