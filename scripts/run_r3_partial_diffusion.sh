#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_dir="${R3_PARTIAL_DIR:-${campaign_dir}/r3/partial_diffusion}"
matrix="${R3_PARTIAL_MATRIX:-${campaign_dir}/config/r3_partial_diffusion.tsv}"
source_pdb="${campaign_dir}/r3/stable_ubiquitin_interface/preparation/ubiquitin_wt_control.pdb"
selected_dir="${phase_dir}/selected_backbones"
sequence_dir="${phase_dir}/sequence"
af2_dir="${phase_dir}/af2"
report_dir="${phase_dir}/reports"
sequence_pdb_dir="${sequence_dir}/output/batch1/ligandmpnn/standardized_pdb"
af2_test="af2_model_1_multimer_tt_3rec"
jsonl="${af2_dir}/output/contig1_batch1/${af2_test}.jsonl"

mkdir -p \
    "${phase_dir}/conditions" \
    "${selected_dir}" \
    "${sequence_dir}" \
    "${af2_dir}" \
    "${report_dir}" \
    "${phase_dir}/work"
exec 9>"${phase_dir}/r3_partial_diffusion.lock"
if ! flock -n 9; then
    echo "Another R3 partial-diffusion driver holds the lock"
    exit 1
fi
if [[ ! -s "${source_pdb}" ]]; then
    echo "Missing stable-scaffold source PDB: ${source_pdb}"
    exit 1
fi

stage_completed() {
    local trace_file="$1"
    [[ -s "${trace_file}" ]] && awk -F '\t' \
        'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' \
        "${trace_file}"
}

while IFS=$'\t' read -r condition_id partial_t num_designs; do
    if [[ "${condition_id}" == "condition_id" ]]; then
        continue
    fi
    condition_dir="${phase_dir}/conditions/${condition_id}"
    backbone_dir="${condition_dir}/backbone"
    metrics_dir="${condition_dir}/metrics"
    work_dir="${phase_dir}/work/${condition_id}"
    mkdir -p "${backbone_dir}" "${metrics_dir}" "${work_dir}/backbone" "${work_dir}/metrics"

    if ! stage_completed "${backbone_dir}/trace.txt"; then
        (
            cd "${backbone_dir}"
            "${nextflow_bin}" \
                -log "${backbone_dir}/nextflow.log" \
                run \
                -with-trace "${backbone_dir}/trace.txt" \
                -with-report "${backbone_dir}/report.html" \
                -work-dir "${work_dir}/backbone" \
                "${ovo_package_dir}/pipelines/rfdiffusion-backbone" \
                --publish_dir "${backbone_dir}/output" \
                --reference_files_dir "${ovo_home_dir}/reference_files" \
                --shared_modules "ovo:${ovo_package_dir}" \
                -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
                -config "${ovo_package_dir}/pipelines/rfdiffusion-backbone/nextflow.config" \
                -profile docker \
                -config "${ovo_home_dir}/nextflow_local.config" \
                --max_memory 512GB \
                -ansi-log false \
                --contig "76-76/0 B6-134" \
                --input_pdb "${source_pdb}" \
                --num_designs "${num_designs}" \
                --hotspot "B50,B52,B53,B55,B57,B61" \
                --run_parameters "diffuser.T=50 diffuser.partial_T=${partial_t}" \
                -resume \
                > "${backbone_dir}/nextflow.stdout.log" 2>&1
        )
    fi
    standardized_dir="${backbone_dir}/output/rfdiffusion/rfdiffusion_standardized_pdb"
    if [[ "$(find -L "${standardized_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne "${num_designs}" ]]; then
        echo "${condition_id}: unexpected standardized backbone count"
        exit 1
    fi

    if ! stage_completed "${metrics_dir}/trace.txt"; then
        (
            cd "${metrics_dir}"
            "${nextflow_bin}" \
                -log "${metrics_dir}/nextflow.log" \
                run \
                -with-trace "${metrics_dir}/trace.txt" \
                -with-report "${metrics_dir}/report.html" \
                -work-dir "${work_dir}/metrics" \
                "${ovo_package_dir}/pipelines/backbone-metrics" \
                --publish_dir "${metrics_dir}/output" \
                --shared_modules "ovo:${ovo_package_dir}" \
                -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
                -config "${ovo_package_dir}/pipelines/backbone-metrics/nextflow.config" \
                -profile docker \
                -config "${ovo_home_dir}/nextflow_local.config" \
                --max_memory 512GB \
                -ansi-log false \
                --output_dir "${condition_id}" \
                --pdb_dir "${standardized_dir}" \
                --hotspot "B50,B52,B53,B55,B57,B61" \
                --filters "N_contact_hotspots>=8,N_hotspots_on_interface>=4" \
                -resume \
                > "${metrics_dir}/nextflow.stdout.log" 2>&1
        )
    fi
    "${python_bin}" "${campaign_dir}/scripts/summarize_backbone_metrics.py" \
        "${metrics_dir}/output/${condition_id}/backbone_metrics.csv" \
        --filtered-dir "${metrics_dir}/output/${condition_id}/backbones_filtered" \
        --json "${report_dir}/${condition_id}_backbone_summary.json"
    "${python_bin}" "${campaign_dir}/scripts/validate_backbones.py" \
        "${standardized_dir}" \
        --json "${report_dir}/${condition_id}_independent_backbones.json" \
        --csv "${report_dir}/${condition_id}_independent_backbones.csv"
done < "${matrix}"

find "${selected_dir}" -maxdepth 1 -type f -name '*.pdb' -delete
"${python_bin}" "${campaign_dir}/scripts/select_r3_partial_backbones.py" \
    --phase-dir "${phase_dir}" \
    --matrix "${matrix}" \
    --output-dir "${selected_dir}" \
    --report "${report_dir}/selection.json" \
    --per-condition 3
selected_count="$(find "${selected_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)"
if [[ "${selected_count}" -eq 0 ]]; then
    touch "${phase_dir}/r3_partial_diffusion.no_backbone_pass"
    echo "No partial-diffusion backbone passed the unchanged hard filters"
    exit 0
fi

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
            --pdb_path "${selected_dir}" \
            --num_seq_per_target 3 \
            --run_parameters "--temperature 0.1 --omit_AA C" \
            -resume \
            > "${sequence_dir}/nextflow.stdout.log" 2>&1
    )
fi
expected_sequences=$((selected_count * 3))
if [[ "$(find -L "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne "${expected_sequences}" ]]; then
    echo "Unexpected partial-diffusion sequence count"
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
            --input_designs "${sequence_pdb_dir}/" \
            --native_pdb "${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
            --tests "${af2_test}" \
            --design_type binder \
            --batch_size 30 \
            -resume \
            > "${af2_dir}/nextflow.stdout.log" 2>&1
    )
fi
if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne "${expected_sequences}" ]]; then
    echo "Unexpected partial-diffusion AF2 record count"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r3_partial_diffusion.py" \
    --jsonl "${jsonl}" \
    --sequence-pdb-dir "${sequence_pdb_dir}" \
    --selection-report "${report_dir}/selection.json" \
    --csv-output "${report_dir}/af2_metrics.csv" \
    --json-output "${report_dir}/summary.json" \
    --fasta-output "${report_dir}/sequences.fasta"

touch "${phase_dir}/r3_partial_diffusion.completed"
echo "R3 partial-diffusion pilot completed"
