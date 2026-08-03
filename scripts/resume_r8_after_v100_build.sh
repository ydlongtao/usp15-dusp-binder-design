#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
image="${R8_BOLTZ_IMAGE:-ovo-boltz-v100:2.2.1}"
build_tmux="${R8_BUILD_TMUX:-usp15-r8-boltz-v100-build}"
queue_tmux="${R8_QUEUE_TMUX:-usp15-r8-queue}"
r8_dir="${campaign_dir}/r8"

while tmux has-session -t "${build_tmux}" 2>/dev/null; do
    sleep 30
done
if ! docker image inspect "${image}" >/dev/null 2>&1; then
    echo "V100-compatible Boltz image build did not produce ${image}"
    exit 1
fi

exec 8>"${r8_dir}/r8_gpu_pipeline.lock"
if ! flock -n 8; then
    echo "Another R8 GPU pipeline holds the project GPU lock"
    exit 1
fi
docker run --rm --gpus all "${image}" python3 - \
    > "${r8_dir}/boltz_v100_gpu_smoke.log" 2>&1 <<'PY'
import importlib.metadata as metadata

import torch

assert metadata.version("boltz") == "2.2.1"
assert torch.__version__.startswith("2.5.1"), torch.__version__
assert "sm_70" in torch.cuda.get_arch_list(), torch.cuda.get_arch_list()
left = torch.randn((512, 512), device="cuda")
right = torch.randn((512, 512), device="cuda")
result = left @ right
torch.cuda.synchronize()
assert torch.isfinite(result).all()
print(
    {
        "torch": torch.__version__,
        "cuda": torch.version.cuda,
        "device": torch.cuda.get_device_name(0),
        "arch": torch.cuda.get_arch_list(),
        "tensor_smoke": True,
    }
)
PY
flock -u 8

if tmux has-session -t "${queue_tmux}" 2>/dev/null; then
    echo "R8 queue is already running"
    exit 0
fi
tmux new-session -d -s "${queue_tmux}" \
    "export USP15_CAMPAIGN_DIR=${campaign_dir} OVO_HOME_DIR=${ovo_home_dir} OVO_ENV_DIR=${ovo_env_dir} R8_BOLTZ_IMAGE=${image}; exec ${campaign_dir}/scripts/run_r8_pipeline_queue.sh > ${r8_dir}/pipeline_queue.log 2>&1"
echo "V100 GPU smoke passed; R8 queue restarted"
