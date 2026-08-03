#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R9_COLABDESIGN_IMAGE:-ovo-colabdesign}"
phase_dir="${R9_CALIBRATION_DIR:-${campaign_dir}/r9/tt_calibration}"
input_dir="${phase_dir}/input"
output_dir="${phase_dir}/output/ptm_model_2_tt"
jsonl="${output_dir}.jsonl"
params_dir="${ovo_home_dir}/reference_files/alphafold_models"
exact_control="${campaign_dir}/r3/native_6dj9_diagnostic/input/ubv15d_native_wt_diagnostic.pdb"
complete_control="${campaign_dir}/r3/native_6dj9_diagnostic/preparation/controls/ubv15d_wt_control.pdb"

mkdir -p "${input_dir}" "${output_dir}" "${phase_dir}/logs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r9_tt_calibration.lock"
if ! flock -n 9; then
    echo "Another R9 target-template calibration driver holds the lock"
    exit 1
fi

for path in "${exact_control}" "${complete_control}" "${params_dir}"; do
    if [[ ! -e "${path}" ]]; then
        echo "Missing R9 calibration input: ${path}"
        exit 1
    fi
done
cp -f "${exact_control}" "${input_dir}/exact_native_6dj9.pdb"
cp -f "${complete_control}" "${input_dir}/complete_3t9l_6dj9_ubv.pdb"

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 6 ]]; then
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

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 6 ]]; then
    echo "Expected six R9 target-template calibration records"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r9_tt_calibration.py" \
    --jsonl "${jsonl}" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r9_tt_calibration.completed"
echo "R9 target-template calibration completed"
