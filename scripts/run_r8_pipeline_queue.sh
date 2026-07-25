#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
export R8_BOLTZ_IMAGE="${R8_BOLTZ_IMAGE:-ovo-boltz-v100:2.2.1}"
r8_dir="${campaign_dir}/r8"
report_dir="${r8_dir}/boltz2_calibration/reports"
verified_asset_dir="${R8_VERIFIED_ASSET_DIR:-${r8_dir}/boltz2_assets_verified}"
checkpoint_dir="${R8_CHECKPOINT_DIR:-${r8_dir}/boltz2_checkpoint_verified}"
cache_dir="${R8_BOLTZ_CACHE_DIR:-${r8_dir}/boltz2_cache}"
checkpoint_tmux="${R8_CHECKPOINT_TMUX:-usp15-r8-checkpoint-range}"
checkpoint_report="${R8_CHECKPOINT_REPORT:-${report_dir}/checkpoint_full_download.json}"
mols_path="${verified_asset_dir}/mols.tar"
checkpoint_path="${checkpoint_dir}/boltz2_conf.ckpt"
affinity_dir="${R8_AFFINITY_DIR:-${r8_dir}/boltz2_affinity_verified}"
affinity_path="${affinity_dir}/boltz2_aff.ckpt"
affinity_report="${R8_AFFINITY_REPORT:-${report_dir}/affinity_full_download.json}"

mkdir -p \
    "${report_dir}" \
    "${verified_asset_dir}" \
    "${affinity_dir}" \
    "${cache_dir}"
exec 9>"${r8_dir}/r8_pipeline_queue.lock"
if ! flock -n 9; then
    echo "Another R8 pipeline queue holds the lock"
    exit 1
fi

while tmux has-session -t "${checkpoint_tmux}" 2>/dev/null; do
    sleep 30
done
if ! "${python_bin}" - "${checkpoint_report}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
raise SystemExit(
    0 if path.is_file() and json.loads(path.read_text())["verified"] is True else 1
)
PY
then
    echo "Checkpoint download did not pass exact size and LFS SHA-256 verification"
    exit 1
fi

if [[ ! -f "${mols_path}" ]]; then
    truncate -s 1855662080 "${mols_path}"
fi
"${python_bin}" "${campaign_dir}/scripts/repair_sparse_hf_asset.py" \
    --file "${mols_path}" \
    --url "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar" \
    --expected-size 1855662080 \
    --expected-sha256 39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7 \
    --report "${report_dir}/mols_full_download.json" \
    --workers 8 \
    --chunk-size 33554432 \
    --retries 4 \
    > "${r8_dir}/mols_full_parallel_resume.log" 2>&1

if ! "${python_bin}" - "${report_dir}/mols_full_download.json" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
raise SystemExit(
    0 if path.is_file() and json.loads(path.read_text())["verified"] is True else 1
)
PY
then
    echo "mols.tar did not pass exact size and LFS SHA-256 verification"
    exit 1
fi

if [[ ! -f "${affinity_path}" ]]; then
    truncate -s 2062139170 "${affinity_path}"
fi
"${python_bin}" "${campaign_dir}/scripts/repair_sparse_hf_asset.py" \
    --file "${affinity_path}" \
    --url "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt" \
    --expected-size 2062139170 \
    --expected-sha256 dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e \
    --report "${affinity_report}" \
    --workers 8 \
    --chunk-size 33554432 \
    --retries 4 \
    > "${r8_dir}/affinity_full_parallel.log" 2>&1

if ! "${python_bin}" - "${affinity_report}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
raise SystemExit(
    0 if path.is_file() and json.loads(path.read_text())["verified"] is True else 1
)
PY
then
    echo "boltz2_aff.ckpt did not pass exact size and LFS SHA-256 verification"
    exit 1
fi

cp -f "${mols_path}" "${cache_dir}/mols.tar"
cp -f "${checkpoint_path}" "${cache_dir}/boltz2_conf.ckpt"
cp -f "${affinity_path}" "${cache_dir}/boltz2_aff.ckpt"

exec 8>"${r8_dir}/r8_gpu_pipeline.lock"
if ! flock -n 8; then
    echo "Another R8 GPU pipeline holds the project GPU lock"
    exit 1
fi

export R8_BOLTZ_CACHE_DIR="${cache_dir}"
"${campaign_dir}/scripts/run_r8_boltz2_smoke.sh"
"${python_bin}" "${campaign_dir}/scripts/summarize_r8_boltz2.py" \
    --phase-dir "${r8_dir}/boltz2_calibration" \
    --exact-native "${campaign_dir}/r3/native_6dj9_diagnostic/input/ubv15d_native_wt_diagnostic.pdb" \
    --complete-target "${campaign_dir}/r3/native_6dj9_diagnostic/preparation/controls/ubv15d_wt_control.pdb" \
    --seeds 0 \
    --output "${report_dir}/seed0_smoke_summary.json"

if ! "${python_bin}" - "${report_dir}/seed0_smoke_summary.json" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1]))
raise SystemExit(0 if summary["records"][0]["both_controls_pass"] is True else 1)
PY
then
    echo "Boltz-2 seed-0 smoke failed one or more unchanged positive-control gates"
    exit 1
fi

"${campaign_dir}/scripts/run_r8_boltz2_calibration.sh"
if ! "${python_bin}" - "${report_dir}/calibration_summary.json" <<'PY'
import json
import sys

raise SystemExit(0 if json.load(open(sys.argv[1]))["calibrated"] is True else 1)
PY
then
    echo "Boltz-2 seeds 0,1,2 did not meet the fixed calibration rule"
    exit 1
fi

"${campaign_dir}/scripts/run_r8_af2_candidate_screen.sh"
"${campaign_dir}/scripts/run_r8_boltz2_candidate_screen.sh"

touch "${r8_dir}/r8_positive_pipeline.completed"
echo "R8 calibrated positive-screen pipeline completed"
