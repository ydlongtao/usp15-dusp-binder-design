#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R10_COLABDESIGN_IMAGE:-ovo-colabdesign}"
phase_dir="${R10_POSITIVE_DIR:-${campaign_dir}/r10/positive_screen}"
input_dir="${campaign_dir}/r7/rescreen_panel/pdb"
output_dir="${phase_dir}/output/ptm_model_2_ct"
jsonl="${output_dir}.jsonl"
params_dir="${ovo_home_dir}/reference_files/alphafold_models"
calibration_summary="${campaign_dir}/r7/calibration/reports/summary.json"

mkdir -p "${output_dir}" "${phase_dir}/logs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r10_positive_screen.lock"
if ! flock -n 9; then
    echo "Another R10 positive-screen driver holds the lock"
    exit 1
fi

if [[ ! -f "${calibration_summary}" ]] || ! "${python_bin}" - "${calibration_summary}" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1]))
model = summary["model_results"].get("model_2_ptm", {})
raise SystemExit(
    0
    if model.get("calibrated") is True
    and model.get("passing_seed_count") == 3
    else 1
)
PY
then
    echo "R7 model-2 pTM interface-template controls are not calibrated 3/3"
    exit 1
fi

expected_inputs="$(find "${input_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)"
if [[ "${expected_inputs}" -ne 52 ]]; then
    echo "Expected the fixed 52-member panel, found ${expected_inputs}"
    exit 1
fi
expected_records="$((expected_inputs * 3))"

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
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
        > "${phase_dir}/logs/ptm_model_2_ct.log" 2>&1
fi

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
    echo "Expected ${expected_records} R10 positive-screen records"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r8_af2_candidates.py" \
    --jsonl "${jsonl}" \
    --expected-inputs "${expected_inputs}" \
    --phase-name "R10 geometry-conditioned model-2 pTM candidate screen" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r10_positive_screen.completed"
echo "R10 geometry-conditioned positive screen completed"
