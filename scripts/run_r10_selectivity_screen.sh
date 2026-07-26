#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
image="${R10_COLABDESIGN_IMAGE:-ovo-colabdesign}"
phase_dir="${R10_SELECTIVITY_DIR:-${campaign_dir}/r10/selectivity}"
positive_dir="${campaign_dir}/r10/positive_screen"
positive_summary="${positive_dir}/reports/summary.json"
positive_jsonl="${positive_dir}/output/ptm_model_2_ct.jsonl"
panel_dir="${campaign_dir}/r7/rescreen_panel/pdb"
homolog_dir="${campaign_dir}/r10/homologs/targets"
input_dir="${phase_dir}/input"
input_report="${phase_dir}/reports/input_preparation.json"
interface_summary="${phase_dir}/reports/interface_summary.json"
params_dir="${ovo_home_dir}/reference_files/alphafold_models"

mkdir -p "${input_dir}" "${phase_dir}/output" "${phase_dir}/logs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r10_selectivity.lock"
if ! flock -n 9; then
    echo "Another R10 selectivity driver holds the lock"
    exit 1
fi

for path in \
    "${positive_summary}" \
    "${positive_jsonl}" \
    "${homolog_dir}/USP4_5CTR_DUSP_aligned_chainB.pdb" \
    "${homolog_dir}/USP11_4MEL_DUSP_aligned_chainB.pdb"
do
    if [[ ! -f "${path}" ]]; then
        echo "Missing R10 selectivity prerequisite: ${path}"
        exit 1
    fi
done

"${python_bin}" "${campaign_dir}/scripts/audit_r10_interfaces.py" \
    --positive-summary "${positive_summary}" \
    --panel-csv "${campaign_dir}/r7/rescreen_panel/reports/panel.csv" \
    --panel-dir "${panel_dir}" \
    --json "${interface_summary}" \
    --csv "${phase_dir}/reports/interface_metrics.csv"

"${python_bin}" "${campaign_dir}/scripts/prepare_r10_selectivity_inputs.py" \
    --positive-summary "${positive_summary}" \
    --interface-summary "${interface_summary}" \
    --panel-dir "${panel_dir}" \
    --usp4-target "${homolog_dir}/USP4_5CTR_DUSP_aligned_chainB.pdb" \
    --usp11-target "${homolog_dir}/USP11_4MEL_DUSP_aligned_chainB.pdb" \
    --output-dir "${input_dir}" \
    --report "${input_report}"

expected_inputs="$(
    "${python_bin}" - "${input_report}" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1]))["positive_count"])
PY
)"
if [[ "${expected_inputs}" -lt 1 ]]; then
    echo "No R10 positive candidates for selectivity"
    exit 1
fi
expected_records="$((expected_inputs * 3))"

for target in USP4 USP11; do
    output_dir="${phase_dir}/output/${target}/ptm_model_2_ct"
    jsonl="${output_dir}.jsonl"
    mkdir -p "${output_dir}"
    if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
        docker run --rm --gpus all \
            -v "${campaign_dir}:${campaign_dir}" \
            -v "${ovo_home_dir}:${ovo_home_dir}" \
            -w "${campaign_dir}" \
            "${image}" \
            python3 "${campaign_dir}/scripts/run_r7_af2_ensemble_eval.py" \
            "${input_dir}/${target}" \
            "${output_dir}" \
            --params "${params_dir}" \
            --architecture ptm \
            --model-number 2 \
            --seeds 0,1,2 \
            --num-recycles 3 \
            --dropout \
            --use-binder-template \
            --use-interface-template \
            > "${phase_dir}/logs/${target}_ptm_model_2_ct.log" 2>&1
    fi
    if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_records}" ]]; then
        echo "${target}: expected ${expected_records} selectivity records"
        exit 1
    fi
done

"${python_bin}" "${campaign_dir}/scripts/summarize_r10_selectivity.py" \
    --positive-jsonl "${positive_jsonl}" \
    --positive-summary "${positive_summary}" \
    --input-report "${input_report}" \
    --usp4-jsonl "${phase_dir}/output/USP4/ptm_model_2_ct.jsonl" \
    --usp11-jsonl "${phase_dir}/output/USP11/ptm_model_2_ct.jsonl" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r10_selectivity.completed"
echo "R10 same-pose selectivity screen completed"
