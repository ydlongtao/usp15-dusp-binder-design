#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <source-pool>"
    exit 2
fi

source_pool="$1"
campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
source_filtered_dir="${campaign_dir}/pilots/${source_pool}/metrics/output/${source_pool}/backbones_filtered"
smoke_dir="${campaign_dir}/smoke/ligandmpnn_af2/${source_pool}"
input_dir="${smoke_dir}/input"
sequence_dir="${smoke_dir}/sequence"
af2_dir="${smoke_dir}/af2"

mkdir -p \
    "${input_dir}" \
    "${sequence_dir}" \
    "${af2_dir}" \
    "${campaign_dir}/work/smoke/sequence" \
    "${campaign_dir}/work/smoke/af2"

source_pdb="$(find -L "${source_filtered_dir}" -maxdepth 1 -type f -name '*.pdb' | sort | head -n 1)"
if [[ -z "${source_pdb}" ]]; then
    echo "No backbone passed the hard filters in ${source_pool}"
    exit 1
fi
ln -sfn "${source_pdb}" "${input_dir}/smoke_backbone.pdb"

sequence_done=false
if [[ -s "${sequence_dir}/trace.txt" ]] && awk -F '\t' 'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' "${sequence_dir}/trace.txt"; then
    sequence_done=true
fi

if [[ "${sequence_done}" == "false" ]]; then
    cd "${sequence_dir}"
    "${nextflow_bin}" \
        -log "${sequence_dir}/nextflow.log" \
        run \
        -with-trace "${sequence_dir}/trace.txt" \
        -with-report "${sequence_dir}/report.html" \
        -work-dir "${campaign_dir}/work/smoke/sequence" \
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
        --num_seq_per_target 3 \
        --run_parameters "--temperature 0.1 --omit_AA C" \
        -resume \
        > "${sequence_dir}/nextflow.stdout.log" 2>&1
fi

sequence_pdb_dir="${sequence_dir}/output/batch1/ligandmpnn/standardized_pdb"
if [[ "$(find "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne 3 ]]; then
    echo "Expected 3 LigandMPNN sequence PDBs in ${sequence_pdb_dir}"
    exit 1
fi

af2_done=false
if [[ -s "${af2_dir}/trace.txt" ]] && awk -F '\t' 'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' "${af2_dir}/trace.txt"; then
    af2_done=true
fi

if [[ "${af2_done}" == "false" ]]; then
    cd "${af2_dir}"
    "${nextflow_bin}" \
        -log "${af2_dir}/nextflow.log" \
        run \
        -with-trace "${af2_dir}/trace.txt" \
        -with-report "${af2_dir}/report.html" \
        -work-dir "${campaign_dir}/work/smoke/af2" \
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
        --tests "af2_model_1_multimer_tt_3rec" \
        --design_type binder \
        --batch_size 20 \
        -resume \
        > "${af2_dir}/nextflow.stdout.log" 2>&1
fi

awk -F '\t' 'NR > 1 {print $4, $5, $6, $8}' "${sequence_dir}/trace.txt"
awk -F '\t' 'NR > 1 {print $4, $5, $6, $8}' "${af2_dir}/trace.txt"
