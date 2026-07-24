#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 <pool> [attempt]"
    exit 2
fi

pool="$1"
attempt="${2:-1}"
if ! [[ "${attempt}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Attempt must be a positive integer: ${attempt}"
    exit 2
fi
campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
pool_table="${campaign_dir}/config/pools.tsv"
input_pdb="${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb"
if [[ "${attempt}" == "1" ]]; then
    run_dir="${campaign_dir}/previews/${pool}"
else
    run_dir="${campaign_dir}/previews/${pool}/attempt_${attempt}"
fi

pool_row="$(awk -F '\t' -v wanted="${pool}" 'NR > 1 && $1 == wanted {print; exit}' "${pool_table}")"
if [[ -z "${pool_row}" ]]; then
    echo "Unknown pool: ${pool}"
    exit 2
fi

IFS=$'\t' read -r pool_name length_range contig model_weights <<< "${pool_row}"
if [[ "${pool_name}" != "${pool}" ]]; then
    echo "Pool resolution failed: ${pool_name}"
    exit 2
fi

mkdir -p "${run_dir}" "${campaign_dir}/work/previews"

if [[ -s "${run_dir}/trace.txt" ]] && awk -F '\t' 'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' "${run_dir}/trace.txt"; then
    echo "Preview already completed: ${pool}"
    exit 0
fi

run_parameters="diffuser.T=15"
if [[ "${model_weights}" == "Complex_beta" ]]; then
    run_parameters="${run_parameters} inference.ckpt_override_path=rfdiffusion_models/Complex_beta_ckpt.pt"
fi

cd "${run_dir}"
"${nextflow_bin}" \
    -log "${run_dir}/nextflow.log" \
    run \
    -with-trace "${run_dir}/trace.txt" \
    -with-report "${run_dir}/report.html" \
    -work-dir "${campaign_dir}/work/previews" \
    "${ovo_package_dir}/pipelines/rfdiffusion-backbone" \
    --publish_dir "${run_dir}/output" \
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
    --num_designs 1 \
    --hotspot "A50,A52,A53,A55,A57,A61" \
    --run_parameters "${run_parameters}" \
    -resume \
    > "${run_dir}/nextflow.stdout.log" 2>&1

awk -F '\t' 'NR > 1 {print $4, $5, $6, $8}' "${run_dir}/trace.txt"
