#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_b_dir="${R2_PHASE_B_DIR:-${campaign_dir}/r2/phase_b}"
phase_b_matrix="${R2_PHASE_B_MATRIX:-${campaign_dir}/config/r2_phase_b.tsv}"
phase_dir="${R2_PHASE_C_DIR:-${campaign_dir}/r2/phase_c}"
phase_matrix="${R2_PHASE_C_MATRIX:-${phase_dir}/r2_phase_c.tsv}"
af2_test="af2_model_1_multimer_tt_3rec"

if [[ ! -f "${phase_b_dir}/phase_b_backbones.completed" ]]; then
    echo "Phase B backbone queue is not complete"
    exit 1
fi

mkdir -p "${phase_dir}/runs" "${phase_dir}/reports" "${phase_dir}/work"
exec 9>"${phase_dir}/phase_c.lock"
if ! flock -n 9; then
    echo "Another R2 Phase C driver already holds ${phase_dir}/phase_c.lock"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/build_r2_phase_c_matrix.py" \
    --phase-b-dir "${phase_b_dir}" \
    --phase-b-matrix "${phase_b_matrix}" \
    --output "${phase_matrix}" \
    --report "${phase_dir}/reports/phase_c_matrix_summary.json"

stage_completed() {
    local trace_file="$1"
    [[ -s "${trace_file}" ]] && awk -F '\t' \
        'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' \
        "${trace_file}"
}

verify_no_binder_cys() {
    local pdb_file="$1"
    if awk '
        $1 == "ATOM" &&
        substr($0, 22, 1) == "A" &&
        substr($0, 18, 3) == "CYS" {found=1}
        END {exit found ? 0 : 1}
    ' "${pdb_file}"; then
        echo "Binder chain contains Cys: ${pdb_file}"
        return 1
    fi
}

run_condition() {
    local run_id="$1"
    local condition_id="$2"
    local candidate_id="$3"
    local backbone_file="$4"
    local temperature="$5"
    local sequence_count="$6"
    local omit_amino_acids="$7"

    local source_pdb="${phase_b_dir}/conditions/${condition_id}/selected_backbones/${backbone_file}"
    local run_dir="${phase_dir}/runs/${run_id}"
    local input_dir="${run_dir}/input"
    local sequence_dir="${run_dir}/sequence"
    local af2_dir="${run_dir}/af2"
    local sequence_pdb_dir="${sequence_dir}/output/batch1/ligandmpnn/standardized_pdb"
    local af2_jsonl="${af2_dir}/output/contig1_batch1/${af2_test}.jsonl"

    if [[ ! -f "${source_pdb}" ]]; then
        echo "Missing selected Phase B backbone: ${source_pdb}"
        return 1
    fi
    mkdir -p \
        "${input_dir}" \
        "${sequence_dir}" \
        "${af2_dir}" \
        "${phase_dir}/work/${run_id}/sequence" \
        "${phase_dir}/work/${run_id}/af2"
    cp -f --remove-destination \
        "${source_pdb}" \
        "${input_dir}/${candidate_id}.pdb"

    {
        echo "run_id=${run_id}"
        echo "condition_id=${condition_id}"
        echo "candidate_id=${candidate_id}"
        echo "backbone_file=${backbone_file}"
        echo "temperature=${temperature}"
        echo "sequences_per_backbone=${sequence_count}"
        echo "omit_amino_acids=${omit_amino_acids}"
        echo "af2_test=${af2_test}"
    } > "${run_dir}/parameters.txt"

    if ! stage_completed "${sequence_dir}/trace.txt"; then
        (
            cd "${sequence_dir}"
            "${nextflow_bin}" \
                -log "${sequence_dir}/nextflow.log" \
                run \
                -with-trace "${sequence_dir}/trace.txt" \
                -with-report "${sequence_dir}/report.html" \
                -work-dir "${phase_dir}/work/${run_id}/sequence" \
                "${ovo_package_dir}/pipelines/ligandmpnn-sequence-design" \
                --publish_dir "${sequence_dir}/output" \
                --shared_modules "ovo:${ovo_package_dir}" \
                -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
                -config "${ovo_package_dir}/pipelines/ligandmpnn-sequence-design/nextflow.config" \
                -profile docker \
                -config "${ovo_home_dir}/nextflow_local.config" \
                --max_memory 512GB \
                -ansi-log false \
                --pdb_path "${input_dir}" \
                --num_seq_per_target "${sequence_count}" \
                --run_parameters "--temperature ${temperature} --omit_AA ${omit_amino_acids}" \
                -resume \
                > "${sequence_dir}/nextflow.stdout.log" 2>&1
        )
    fi

    local observed_sequence_count
    observed_sequence_count="$(
        find -L "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l
    )"
    if [[ "${observed_sequence_count}" -ne "${sequence_count}" ]]; then
        echo "Expected ${sequence_count} sequence PDBs for ${run_id}; found ${observed_sequence_count}"
        return 1
    fi
    while IFS= read -r sequence_pdb; do
        verify_no_binder_cys "${sequence_pdb}"
    done < <(
        find -L "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' | sort
    )

    if ! stage_completed "${af2_dir}/trace.txt"; then
        (
            cd "${af2_dir}"
            "${nextflow_bin}" \
                -log "${af2_dir}/nextflow.log" \
                run \
                -with-trace "${af2_dir}/trace.txt" \
                -with-report "${af2_dir}/report.html" \
                -work-dir "${phase_dir}/work/${run_id}/af2" \
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
                --batch_size 20 \
                -resume \
                > "${af2_dir}/nextflow.stdout.log" 2>&1
        )
    fi

    if [[ ! -s "${af2_jsonl}" ]]; then
        echo "Missing AF2 metrics JSONL for ${run_id}: ${af2_jsonl}"
        return 1
    fi
    local observed_af2_count
    observed_af2_count="$(wc -l < "${af2_jsonl}")"
    if [[ "${observed_af2_count}" -ne "${sequence_count}" ]]; then
        echo "Expected ${sequence_count} AF2 records for ${run_id}; found ${observed_af2_count}"
        return 1
    fi
}

while IFS=$'\t' read -r \
    run_id \
    condition_id \
    candidate_id \
    backbone_file \
    temperature \
    sequences_per_backbone \
    omit_amino_acids
do
    if [[ "${run_id}" == "run_id" ]]; then
        continue
    fi
    echo "Starting ${run_id}: ${candidate_id}, temperature ${temperature}"
    run_condition \
        "${run_id}" \
        "${condition_id}" \
        "${candidate_id}" \
        "${backbone_file}" \
        "${temperature}" \
        "${sequences_per_backbone}" \
        "${omit_amino_acids}"
done < "${phase_matrix}"

"${python_bin}" "${campaign_dir}/scripts/summarize_r2_phase_c.py" \
    --phase-dir "${phase_dir}" \
    --matrix "${phase_matrix}" \
    --csv-output "${phase_dir}/reports/af2_metrics.csv" \
    --json-output "${phase_dir}/reports/phase_c_summary.json"

touch "${phase_dir}/phase_c.completed"
echo "R2 Phase C completed; thresholds were not relaxed"
