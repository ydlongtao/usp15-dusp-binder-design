#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_dir="${R3_UBIQUITIN_DIR:-${campaign_dir}/r3/stable_ubiquitin_interface}"
preparation_dir="${phase_dir}/preparation"
sequence_dir="${phase_dir}/sequence"
af2_input_dir="${phase_dir}/af2_input"
af2_dir="${phase_dir}/af2"
report_dir="${phase_dir}/reports"
sequence_pdb_dir="${sequence_dir}/output/batch1/ligandmpnn/standardized_pdb"
af2_test="af2_model_1_multimer_tt_3rec"
jsonl="${af2_dir}/output/contig1_batch1/${af2_test}.jsonl"

mkdir -p \
    "${preparation_dir}" \
    "${sequence_dir}/input" \
    "${af2_input_dir}" \
    "${af2_dir}" \
    "${report_dir}" \
    "${phase_dir}/work/sequence" \
    "${phase_dir}/work/af2"
exec 9>"${phase_dir}/r3_ubiquitin_interface.lock"
if ! flock -n 9; then
    echo "Another R3 ubiquitin-interface driver holds the lock"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/prepare_r3_ubiquitin_scaffold.py" \
    --ubiquitin "${campaign_dir}/inputs/1UBQ.pdb" \
    --complex "${campaign_dir}/inputs/6DJ9.pdb" \
    --target-reference "${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
    --wt-output "${preparation_dir}/ubiquitin_wt_control.pdb" \
    --design-output "${sequence_dir}/input/ubiquitin_interface_design.pdb" \
    --report "${report_dir}/preparation.json"

stage_completed() {
    local trace_file="$1"
    [[ -s "${trace_file}" ]] && awk -F '\t' \
        'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' \
        "${trace_file}"
}

if ! stage_completed "${sequence_dir}/trace.txt"; then
    (
        cd "${sequence_dir}"
        "${nextflow_bin}" \
            -log "${sequence_dir}/nextflow.log" \
            run \
            -with-trace "${sequence_dir}/trace.txt" \
            -with-report "${sequence_dir}/report.html" \
            -work-dir "${phase_dir}/work/sequence" \
            "${ovo_package_dir}/pipelines/ligandmpnn-sequence-design" \
            --publish_dir "${sequence_dir}/output" \
            --shared_modules "ovo:${ovo_package_dir}" \
            -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
            -config "${ovo_package_dir}/pipelines/ligandmpnn-sequence-design/nextflow.config" \
            -profile docker \
            -config "${ovo_home_dir}/nextflow_local.config" \
            --max_memory 512GB \
            -ansi-log false \
            --pdb_path "${sequence_dir}/input" \
            --num_seq_per_target 3 \
            --run_parameters "--temperature 0.1 --omit_AA C" \
            -resume \
            > "${sequence_dir}/nextflow.stdout.log" 2>&1
    )
fi
if [[ "$(find -L "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne 3 ]]; then
    echo "Expected three partial-interface LigandMPNN outputs"
    exit 1
fi

find "${af2_input_dir}" -maxdepth 1 -type f -name '*.pdb' -delete
cp "${preparation_dir}/ubiquitin_wt_control.pdb" "${af2_input_dir}/"
cp "${sequence_pdb_dir}"/*.pdb "${af2_input_dir}/"
if [[ "$(find "${af2_input_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne 4 ]]; then
    echo "Expected one control plus three interface designs"
    exit 1
fi

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
            --input_designs "${af2_input_dir}/" \
            --native_pdb "${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
            --tests "${af2_test}" \
            --design_type binder \
            --batch_size 20 \
            -resume \
            > "${af2_dir}/nextflow.stdout.log" 2>&1
    )
fi
if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 4 ]]; then
    echo "Expected four stable-scaffold AF2 records"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r3_ubiquitin_interface.py" \
    --jsonl "${jsonl}" \
    --af2-input-dir "${af2_input_dir}" \
    --preparation-report "${report_dir}/preparation.json" \
    --csv-output "${report_dir}/af2_metrics.csv" \
    --json-output "${report_dir}/summary.json" \
    --fasta-output "${report_dir}/sequences.fasta"

touch "${phase_dir}/r3_ubiquitin_interface.completed"
echo "R3 stable ubiquitin interface stage completed"
