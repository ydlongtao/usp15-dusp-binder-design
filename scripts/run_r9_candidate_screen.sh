#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R9_COLABDESIGN_IMAGE:-ovo-colabdesign}"
phase_dir="${R9_POSITIVE_DIR:-${campaign_dir}/r9/positive_screen}"
input_dir="${campaign_dir}/r7/rescreen_panel/pdb"
output_dir="${phase_dir}/output/ptm_model_2_tt"
jsonl="${output_dir}.jsonl"
params_dir="${ovo_home_dir}/reference_files/alphafold_models"
calibration_summary="${campaign_dir}/r9/tt_calibration/reports/summary.json"

mkdir -p "${output_dir}" "${phase_dir}/logs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r9_positive_screen.lock"
if ! flock -n 9; then
    echo "Another R9 positive-screen driver holds the lock"
    exit 1
fi

if [[ ! -f "${calibration_summary}" ]] || ! "${python_bin}" - "${calibration_summary}" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1]))
raise SystemExit(
    0
    if summary["calibrated"] is True
    and summary["model"] == "model_2_ptm"
    and summary["template_mode"] == "tt"
    else 1
)
PY
then
    echo "R9 model-2 pTM target-template controls are not calibrated"
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
        > "${phase_dir}/logs/ptm_model_2_tt.log" 2>&1
fi

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
    echo "Expected ${expected_records} R9 positive-screen records"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r9_candidates.py" \
    --jsonl "${jsonl}" \
    --expected-inputs "${expected_inputs}" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r9_positive_screen.completed"
echo "R9 positive screen completed"
