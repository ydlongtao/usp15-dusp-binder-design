#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <pool>"
    exit 2
fi

pool="$1"
campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
pool_table="${campaign_dir}/config/pools.tsv"
input_pdb="${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb"
pilot_dir="${campaign_dir}/pilots/${pool}"
backbone_dir="${pilot_dir}/backbone"
metrics_dir="${pilot_dir}/metrics"

pool_row="$(awk -F '\t' -v wanted="${pool}" 'NR > 1 && $1 == wanted {print; exit}' "${pool_table}")"
if [[ -z "${pool_row}" ]]; then
    echo "Unknown pool: ${pool}"
    exit 2
fi
IFS=$'\t' read -r pool_name length_range contig model_weights <<< "${pool_row}"

mkdir -p \
    "${backbone_dir}" \
    "${metrics_dir}" \
    "${campaign_dir}/work/pilots/backbone" \
    "${campaign_dir}/work/pilots/metrics" \
    "${campaign_dir}/reports"

backbone_done=false
if [[ -s "${backbone_dir}/trace.txt" ]] && awk -F '\t' 'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' "${backbone_dir}/trace.txt"; then
    backbone_done=true
fi

if [[ "${backbone_done}" == "false" ]]; then
    run_parameters="diffuser.T=50"
    if [[ "${model_weights}" == "Complex_beta" ]]; then
        run_parameters="${run_parameters} inference.ckpt_override_path=rfdiffusion_models/Complex_beta_ckpt.pt"
    fi

    cd "${backbone_dir}"
    "${nextflow_bin}" \
        -log "${backbone_dir}/nextflow.log" \
        run \
        -with-trace "${backbone_dir}/trace.txt" \
        -with-report "${backbone_dir}/report.html" \
        -work-dir "${campaign_dir}/work/pilots/backbone" \
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
        --contig "${contig}" \
        --input_pdb "${input_pdb}" \
        --num_designs 100 \
        --hotspot "A50,A52,A53,A55,A57,A61" \
        --run_parameters "${run_parameters}" \
        -resume \
        > "${backbone_dir}/nextflow.stdout.log" 2>&1
fi

standardized_dir="${backbone_dir}/output/rfdiffusion/rfdiffusion_standardized_pdb"
if [[ "$(find -L "${standardized_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne 100 ]]; then
    echo "Expected 100 standardized backbones in ${standardized_dir}"
    exit 1
fi

metrics_done=false
if [[ -s "${metrics_dir}/trace.txt" ]] && awk -F '\t' 'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' "${metrics_dir}/trace.txt"; then
    metrics_done=true
fi

if [[ "${metrics_done}" == "false" ]]; then
    cd "${metrics_dir}"
    "${nextflow_bin}" \
        -log "${metrics_dir}/nextflow.log" \
        run \
        -with-trace "${metrics_dir}/trace.txt" \
        -with-report "${metrics_dir}/report.html" \
        -work-dir "${campaign_dir}/work/pilots/metrics" \
        "${ovo_package_dir}/pipelines/backbone-metrics" \
        --publish_dir "${metrics_dir}/output" \
        --shared_modules "ovo:${ovo_package_dir}" \
        -config "${ovo_package_dir}/pipelines/nextflow_default.config" \
        -config "${ovo_package_dir}/pipelines/backbone-metrics/nextflow.config" \
        -profile docker \
        -config "${ovo_home_dir}/nextflow_local.config" \
        --max_memory 512GB \
        -ansi-log false \
        --output_dir "${pool}" \
        --pdb_dir "${standardized_dir}" \
        --hotspot "B50,B52,B53,B55,B57,B61" \
        --filters "N_contact_hotspots>=8,N_hotspots_on_interface>=4" \
        -resume \
        > "${metrics_dir}/nextflow.stdout.log" 2>&1
fi

metrics_csv="${metrics_dir}/output/${pool}/backbone_metrics.csv"
filtered_dir="${metrics_dir}/output/${pool}/backbones_filtered"
python3 "${campaign_dir}/scripts/summarize_backbone_metrics.py" \
    "${metrics_csv}" \
    --filtered-dir "${filtered_dir}" \
    --json "${campaign_dir}/reports/${pool}.pilot_summary.json"
