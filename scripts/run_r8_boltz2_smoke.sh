#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R8_BOLTZ_IMAGE:-ovo-boltz}"
phase_dir="${R8_BOLTZ_DIR:-${campaign_dir}/r8/boltz2_calibration}"
input_dir="${phase_dir}/input"
cache_dir="${R8_BOLTZ_CACHE_DIR:-${ovo_home_dir}/reference_files/boltz_models}"
seed=0
run_dir="${phase_dir}/seed_${seed}"
predictions_dir="${run_dir}/boltz_results_$(basename "${input_dir}")/predictions"

mkdir -p "${input_dir}" "${cache_dir}" "${run_dir}" "${phase_dir}/reports"
exec 9>"${phase_dir}/r8_boltz2.lock"
if ! flock -n 9; then
    echo "Another R8 Boltz-2 driver holds the lock"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/prepare_r8_boltz_controls.py" \
    --exact-native "${campaign_dir}/r3/native_6dj9_diagnostic/input/ubv15d_native_wt_diagnostic.pdb" \
    --complete-target "${campaign_dir}/r3/native_6dj9_diagnostic/preparation/controls/ubv15d_wt_control.pdb" \
    --output-dir "${input_dir}" \
    --report "${phase_dir}/reports/input_preparation.json"

if [[ "$(find "${predictions_dir}" -type f -name '*_model_0.pdb' 2>/dev/null | wc -l)" -ne 2 ]]; then
    docker run --rm --gpus all --shm-size 8g \
        -e NUMBA_CACHE_DIR=/tmp \
        -v "${campaign_dir}:${campaign_dir}" \
        -v "${ovo_home_dir}:${ovo_home_dir}" \
        -w "${campaign_dir}" \
        "${image}" \
        boltz predict "${input_dir}" \
        --out_dir "${run_dir}" \
        --cache "${cache_dir}" \
        --accelerator gpu \
        --model boltz2 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 1 \
        --seed "${seed}" \
        --use_msa_server \
        --output_format pdb \
        --write_full_pae \
        --override \
        > "${phase_dir}/seed_${seed}.log" 2>&1
fi

prediction_count="$(find "${predictions_dir}" -type f -name '*_model_0.pdb' | wc -l)"
if [[ "${prediction_count}" -ne 2 ]]; then
    echo "Expected two Boltz-2 seed-0 PDB predictions, found ${prediction_count}"
    exit 1
fi

touch "${phase_dir}/r8_boltz2_smoke.completed"
echo "R8 Boltz-2 seed-0 smoke completed"
