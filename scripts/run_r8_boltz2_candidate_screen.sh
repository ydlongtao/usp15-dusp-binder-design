#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R8_BOLTZ_IMAGE:-ovo-boltz}"
phase_dir="${R8_BOLTZ_CANDIDATE_DIR:-${campaign_dir}/r8/boltz2_candidates}"
cache_dir="${R8_BOLTZ_CACHE_DIR:-${campaign_dir}/r8/boltz2_cache}"
input_dir="${phase_dir}/input"
input_report="${phase_dir}/reports/input_preparation.json"
af2_summary="${campaign_dir}/r8/af2_candidates/reports/summary.json"
panel_dir="${campaign_dir}/r7/rescreen_panel/pdb"
calibration_summary="${campaign_dir}/r8/boltz2_calibration/reports/calibration_summary.json"

mkdir -p "${input_dir}" "${cache_dir}" "${phase_dir}/reports"
exec 9>"${phase_dir}/r8_boltz2_candidates.lock"
if ! flock -n 9; then
    echo "Another R8 Boltz-2 candidate driver holds the lock"
    exit 1
fi

for summary in "${calibration_summary}" "${af2_summary}"; do
    if [[ ! -f "${summary}" ]]; then
        echo "Missing prerequisite summary: ${summary}"
        exit 1
    fi
done
if ! "${python_bin}" - "${calibration_summary}" <<'PY'
import json
import sys

raise SystemExit(0 if json.load(open(sys.argv[1]))["calibrated"] is True else 1)
PY
then
    echo "Boltz-2 controls have not met the fixed R8 calibration rule"
    exit 1
fi

rm -f "${input_dir}"/*.yaml
"${python_bin}" "${campaign_dir}/scripts/prepare_r8_boltz_candidate_inputs.py" \
    --af2-summary "${af2_summary}" \
    --pdb-dir "${panel_dir}" \
    --output-dir "${input_dir}" \
    --report "${input_report}"
expected_inputs="$(find "${input_dir}" -maxdepth 1 -type f -name '*.yaml' | wc -l)"

for seed in 0 1 2; do
    run_dir="${phase_dir}/seed_${seed}"
    mkdir -p "${run_dir}"
    prediction_count="$(
        find "${run_dir}/predictions" -type f -name '*_model_0.pdb' 2>/dev/null \
            | wc -l
    )"
    if [[ "${prediction_count}" -ne "${expected_inputs}" ]]; then
        docker run --rm --gpus all \
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
    prediction_count="$(
        find "${run_dir}/predictions" -type f -name '*_model_0.pdb' | wc -l
    )"
    if [[ "${prediction_count}" -ne "${expected_inputs}" ]]; then
        echo "Seed ${seed}: expected ${expected_inputs} predictions, found ${prediction_count}"
        exit 1
    fi
done

"${python_bin}" "${campaign_dir}/scripts/summarize_r8_boltz2_candidates.py" \
    --phase-dir "${phase_dir}" \
    --input-report "${input_report}" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r8_boltz2_candidates.completed"
echo "R8 Boltz-2 candidate screen completed"
