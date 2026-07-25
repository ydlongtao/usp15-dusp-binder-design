#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R7_COLABDESIGN_IMAGE:-ovo-colabdesign}"
phase_dir="${R8_AF2_CANDIDATE_DIR:-${campaign_dir}/r8/af2_candidates}"
panel_dir="${campaign_dir}/r7/rescreen_panel"
input_dir="${panel_dir}/pdb"
output_dir="${phase_dir}/output/ptm_model_2"
jsonl="${output_dir}.jsonl"
params_dir="${ovo_home_dir}/reference_files/alphafold_models"
calibration_summary="${campaign_dir}/r8/boltz2_calibration/reports/calibration_summary.json"

mkdir -p "${phase_dir}/output" "${phase_dir}/logs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r8_af2_candidates.lock"
if ! flock -n 9; then
    echo "Another R8 AF2 candidate driver holds the lock"
    exit 1
fi

if [[ ! -f "${calibration_summary}" ]] || ! "${python_bin}" - "${calibration_summary}" <<'PY'
import json
import sys

raise SystemExit(0 if json.load(open(sys.argv[1]))["calibrated"] is True else 1)
PY
then
    echo "Boltz-2 controls have not met the fixed R8 calibration rule"
    exit 1
fi

expected_inputs="$(find "${input_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)"
if [[ "${expected_inputs}" -lt 1 ]]; then
    echo "No candidate PDBs in ${input_dir}"
    exit 1
fi
expected_records="$((expected_inputs * 3))"

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
    mkdir -p "${output_dir}"
    docker run --rm --gpus all \
        -v "${campaign_dir}:${campaign_dir}" \
        -v "${ovo_home_dir}:${ovo_home_dir}" \
        -w "${campaign_dir}" \
        "${image}" \
        python3 "${campaign_dir}/scripts/run_r7_af2_ensemble_eval.py" \
        "${input_dir}" \
        "${output_dir}" \
        --params "${params_dir}" \
        --architecture ptm \
        --model-number 2 \
        --seeds 0,1,2 \
        --num-recycles 3 \
        --dropout \
        --use-binder-template \
        --use-interface-template \
        > "${phase_dir}/logs/ptm_model_2.log" 2>&1
fi

if [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
    echo "Expected ${expected_records} AF2 records"
    exit 1
fi
"${python_bin}" "${campaign_dir}/scripts/summarize_r8_af2_candidates.py" \
    --jsonl "${jsonl}" \
    --expected-inputs "${expected_inputs}" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r8_af2_candidates.completed"
echo "R8 AF2 candidate screen completed"
