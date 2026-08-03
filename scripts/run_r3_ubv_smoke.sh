#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_dir="${R3_UBV_DIR:-${campaign_dir}/r3/ubv_positive_smoke}"
source_pdb="${campaign_dir}/inputs/6DJ9.pdb"
input_dir="${phase_dir}/input"
af2_dir="${phase_dir}/af2"
report_dir="${phase_dir}/reports"
af2_test="af2_model_1_multimer_tt_3rec"
jsonl="${af2_dir}/output/contig1_batch1/${af2_test}.jsonl"

mkdir -p \
    "${input_dir}" \
    "${af2_dir}" \
    "${report_dir}" \
    "${phase_dir}/work/af2"
exec 9>"${phase_dir}/r3_ubv_smoke.lock"
if ! flock -n 9; then
    echo "Another R3 UbV smoke driver holds the lock"
    exit 1
fi

if [[ ! -s "${source_pdb}" ]]; then
    echo "Missing 6DJ9 input: ${source_pdb}"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/prepare_r3_ubv_controls.py" \
    --source "${source_pdb}" \
    --target-reference "${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
    --output-dir "${input_dir}" \
    --report "${report_dir}/preparation.json" \
    --fasta "${report_dir}/controls.fasta"

stage_completed() {
    local trace_file="$1"
    [[ -s "${trace_file}" ]] && awk -F '\t' \
        'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' \
        "${trace_file}"
}

if ! stage_completed "${af2_dir}/trace.txt"; then
    (
        cd "${af2_dir}"
        "${nextflow_bin}" \
            -log "${af2_dir}/nextflow.log" \
            run \
            -with-trace "${af2_dir}/trace.txt" \
            -with-report "${af2_dir}/report.html" \
            -work-dir "${phase_dir}/work/af2" \
            "${ovo_package_dir}/pipelines/refolding" \
            --publish_dir "${af2_dir}/output" \
            --reference_files_dir "${ovo_home_dir}/reference_files" \
            --shared_modules "ovo:${ovo_package_dir}" \
            -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
            -config "${ovo_package_dir}/pipelines/refolding/nextflow.config" \
            -profile docker \
            -config "${ovo_home_dir}/nextflow_local.config" \
            --max_memory 512GB \
            -ansi-log false \
            --input_designs "${input_dir}/" \
            --native_pdb "${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
            --tests "${af2_test}" \
            --design_type binder \
            --batch_size 20 \
            -resume \
            > "${af2_dir}/nextflow.stdout.log" 2>&1
    )
fi

if [[ ! -s "${jsonl}" ]]; then
    echo "Missing R3 AF2 JSONL: ${jsonl}"
    exit 1
fi
if [[ "$(wc -l < "${jsonl}")" -ne 4 ]]; then
    echo "Expected four R3 AF2 records"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r3_ubv_smoke.py" \
    --jsonl "${jsonl}" \
    --preparation-report "${report_dir}/preparation.json" \
    --csv-output "${report_dir}/af2_metrics.csv" \
    --json-output "${report_dir}/summary.json"

touch "${phase_dir}/r3_ubv_smoke.completed"
echo "R3 UbV positive-control smoke completed"
