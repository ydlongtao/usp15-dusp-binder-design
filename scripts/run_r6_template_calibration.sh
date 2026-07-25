#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_dir="${R6_CALIBRATION_DIR:-${campaign_dir}/r6/template_calibration}"
source_pdb="${R6_NATIVE_INPUT:-${campaign_dir}/r3/native_6dj9_diagnostic/input/ubv15d_native_wt_diagnostic.pdb}"
input_dir="${phase_dir}/input"
input_pdb="${input_dir}/$(basename "${source_pdb}")"
control_label="${R6_CONTROL_LABEL:-exact_6dj9_native}"

tests=(
    af2_model_1_ptm_tbt_3rec
    af2_model_1_multimer_tbt_3rec
    af2_model_1_ptm_ct_3rec
    af2_model_1_multimer_ct_3rec
)

mkdir -p "${input_dir}" "${phase_dir}/runs" "${phase_dir}/reports"
exec 9>"${phase_dir}/r6_template_calibration.lock"
if ! flock -n 9; then
    echo "Another R6 template-calibration driver holds the lock"
    exit 1
fi

if [[ ! -s "${source_pdb}" ]]; then
    echo "Missing exact-native positive-control PDB: ${source_pdb}"
    exit 1
fi
cp -f "${source_pdb}" "${input_pdb}"

for test_name in "${tests[@]}"; do
    run_dir="${phase_dir}/runs/${test_name}"
    output_dir="${run_dir}/output"
    jsonl="${output_dir}/contig1_batch1/${test_name}.jsonl"
    mkdir -p "${run_dir}" "${run_dir}/work"

    if [[ -s "${jsonl}" ]] && [[ "$(wc -l < "${jsonl}")" -eq 1 ]]; then
        echo "Reusing completed ${test_name}"
        continue
    fi

    (
        cd "${run_dir}"
        "${nextflow_bin}" \
            -log "${run_dir}/nextflow.log" \
            run \
            -with-trace "${run_dir}/trace.txt" \
            -with-report "${run_dir}/report.html" \
            -work-dir "${run_dir}/work" \
            "${ovo_package_dir}/pipelines/refolding" \
            --publish_dir "${output_dir}" \
            --reference_files_dir "${ovo_home_dir}/reference_files" \
            --shared_modules "ovo:${ovo_package_dir}" \
            -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
            -config "${ovo_package_dir}/pipelines/refolding/nextflow.config" \
            -profile docker \
            -config "${ovo_home_dir}/nextflow_local.config" \
            --max_memory 512GB \
            -ansi-log false \
            --input_designs "${input_pdb}" \
            --tests "${test_name}" \
            --design_type binder \
            --batch_size 20 \
            -resume \
            > "${run_dir}/nextflow.stdout.log" 2>&1
    )

    if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 1 ]]; then
        echo "Expected one AF2 record for ${test_name}"
        exit 1
    fi
done

"${python_bin}" "${campaign_dir}/scripts/summarize_r6_template_calibration.py" \
    --phase-dir "${phase_dir}" \
    --control-label "${control_label}" \
    --json "${phase_dir}/reports/summary.json" \
    --csv "${phase_dir}/reports/metrics.csv"

touch "${phase_dir}/r6_template_calibration.completed"
echo "R6 template calibration completed"
