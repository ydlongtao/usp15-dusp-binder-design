#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R7_COLABDESIGN_IMAGE:-ovo-colabdesign}"
phase_dir="${R7_CALIBRATION_DIR:-${campaign_dir}/r7/calibration}"
input_dir="${phase_dir}/input"
output_root="${phase_dir}/output"
params_dir="${ovo_home_dir}/reference_files/alphafold_models"
exact_control="${campaign_dir}/r3/native_6dj9_diagnostic/input/ubv15d_native_wt_diagnostic.pdb"
complete_control="${campaign_dir}/r3/native_6dj9_diagnostic/preparation/controls/ubv15d_wt_control.pdb"

models=(
    "ptm:1"
    "ptm:2"
    "multimer:1"
    "multimer:2"
    "multimer:3"
    "multimer:4"
    "multimer:5"
)

mkdir -p "${input_dir}" "${output_root}" "${phase_dir}/logs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r7_calibration.lock"
if ! flock -n 9; then
    echo "Another R7 calibration driver holds the lock"
    exit 1
fi

for path in "${exact_control}" "${complete_control}" "${params_dir}"; do
    if [[ ! -e "${path}" ]]; then
        echo "Missing R7 input: ${path}"
        exit 1
    fi
done
cp -f "${exact_control}" "${input_dir}/exact_native_6dj9.pdb"
cp -f "${complete_control}" "${input_dir}/complete_3t9l_6dj9_ubv.pdb"

for spec in "${models[@]}"; do
    architecture="${spec%%:*}"
    model_number="${spec##*:}"
    model_id="${architecture}_model_${model_number}"
    output_dir="${output_root}/${model_id}"
    jsonl="${output_dir}.jsonl"
    if [[ -s "${jsonl}" ]] && [[ "$(wc -l < "${jsonl}")" -eq 6 ]]; then
        echo "Reusing completed ${model_id}"
        continue
    fi
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
        --architecture "${architecture}" \
        --model-number "${model_number}" \
        --seeds 0,1,2 \
        --num-recycles 3 \
        --dropout \
        --use-binder-template \
        --use-interface-template \
        > "${phase_dir}/logs/${model_id}.log" 2>&1

    if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 6 ]]; then
        echo "Expected six records for ${model_id}"
        exit 1
    fi
done

"${python_bin}" "${campaign_dir}/scripts/summarize_r7_calibration.py" \
    --phase-dir "${phase_dir}" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r7_calibration.completed"
echo "R7 calibration completed"
