#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <condition-id>"
    exit 2
fi

condition_id="$1"
campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_matrix="${R2_PHASE_B_MATRIX:-${campaign_dir}/config/r2_phase_b.tsv}"
phase_dir="${R2_PHASE_B_DIR:-${campaign_dir}/r2/phase_b}"
input_pdb="${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb"

condition_row="$(
    awk -F '\t' -v wanted="${condition_id}" \
        'NR > 1 && $1 == wanted {print; exit}' \
        "${phase_matrix}"
)"
if [[ -z "${condition_row}" ]]; then
    echo "Unknown Phase B condition: ${condition_id}"
    exit 2
fi

IFS=$'\t' read -r \
    parsed_condition \
    mode \
    length_range \
    contig \
    model_weights \
    generation_hotspots \
    configured_backbones \
    scaffold_min_length \
    scaffold_max_length <<< "${condition_row}"

num_designs="${R2_PHASE_B_NUM_DESIGNS:-${configured_backbones}}"
condition_dir="${phase_dir}/conditions/${condition_id}"
backbone_dir="${condition_dir}/backbone"
metrics_dir="${condition_dir}/metrics"
standardized_dir="${backbone_dir}/output/rfdiffusion/rfdiffusion_standardized_pdb"
resource_dir="${phase_dir}/resources"
model_dir="${ovo_home_dir}/reference_files/rfdiffusion_models"

mkdir -p \
    "${backbone_dir}" \
    "${metrics_dir}" \
    "${condition_dir}/reports" \
    "${phase_dir}/work/${condition_id}/backbone" \
    "${phase_dir}/work/${condition_id}/metrics"

stage_completed() {
    local trace_file="$1"
    [[ -s "${trace_file}" ]] && awk -F '\t' \
        'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' \
        "${trace_file}"
}

run_parameters="diffuser.T=50"
if [[ "${model_weights}" == "Complex_beta" ]]; then
    run_parameters+=" inference.ckpt_override_path=rfdiffusion_models/Complex_beta_ckpt.pt"
elif [[ "${model_weights}" != "Complex_base" ]]; then
    echo "Unsupported model weights: ${model_weights}"
    exit 1
fi

if [[ "${mode}" == "scaffold" ]]; then
    scaffold_set="scaffolds_${scaffold_min_length}_${scaffold_max_length}"
    scaffold_dir="${resource_dir}/${scaffold_set}"
    target_ss="${resource_dir}/target_folds/USP15_DUSP_3T9L_A6-134_ss.pt"
    target_adj="${resource_dir}/target_folds/USP15_DUSP_3T9L_A6-134_adj.pt"
    for required in "${scaffold_dir}" "${target_ss}" "${target_adj}"; do
        if [[ ! -e "${required}" ]]; then
            echo "Missing scaffold-guided resource: ${required}"
            exit 1
        fi
    done
    scaffold_checkpoint="${model_dir}/Complex_Fold_base_ckpt.pt"
    if [[ ! -s "${scaffold_checkpoint}" ]]; then
        echo "Missing official scaffold checkpoint: ${scaffold_checkpoint}"
        exit 1
    fi
    run_parameters+=" scaffoldguided.scaffoldguided_enable=True"
    run_parameters+=" scaffoldguided.target_pdb=True"
    run_parameters+=" scaffoldguided.target_path=${input_pdb}"
    run_parameters+=" scaffoldguided.target_ss=${target_ss}"
    run_parameters+=" scaffoldguided.target_adj=${target_adj}"
    run_parameters+=" scaffoldguided.scaffold_dir=${scaffold_dir}"
    run_parameters+=" scaffoldguided.mask_loops=False"
    run_parameters+=" denoiser.noise_scale_ca=0"
    run_parameters+=" denoiser.noise_scale_frame=0"
elif [[ "${mode}" != "standard" ]]; then
    echo "Unsupported Phase B mode: ${mode}"
    exit 1
fi

{
    echo "condition_id=${condition_id}"
    echo "mode=${mode}"
    echo "length_range=${length_range}"
    echo "contig=${contig}"
    echo "model_weights=${model_weights}"
    echo "generation_hotspots=${generation_hotspots}"
    echo "evaluation_hotspots=A50,A52,A53,A55,A57,A61"
    echo "num_designs=${num_designs}"
    echo "run_parameters=${run_parameters}"
} > "${condition_dir}/parameters.txt"

if [[ "${mode}" == "standard" ]] && ! stage_completed "${backbone_dir}/trace.txt"; then
    (
        cd "${backbone_dir}"
        "${nextflow_bin}" \
            -log "${backbone_dir}/nextflow.log" \
            run \
            -with-trace "${backbone_dir}/trace.txt" \
            -with-report "${backbone_dir}/report.html" \
            -work-dir "${phase_dir}/work/${condition_id}/backbone" \
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
            --num_designs "${num_designs}" \
            --hotspot "${generation_hotspots}" \
            --run_parameters "${run_parameters}" \
            -resume \
            > "${backbone_dir}/nextflow.stdout.log" 2>&1
    )
fi

if [[ "${mode}" == "scaffold" ]]; then
    raw_dir="${backbone_dir}/direct_raw"
    published_root="${backbone_dir}/output/rfdiffusion"
    pdb_dir="${published_root}/rfdiffusion_pdb"
    raw_trb_dir="${published_root}/rfdiffusion_trb_raw"
    trb_dir="${published_root}/rfdiffusion_trb"
    traj_dir="${published_root}/rfdiffusion_traj"
    mkdir -p \
        "${raw_dir}" \
        "${pdb_dir}" \
        "${raw_trb_dir}" \
        "${trb_dir}" \
        "${traj_dir}" \
        "${standardized_dir}"

    existing_pdbs="$(find "${raw_dir}" -maxdepth 1 -type f -name 'rfdiffusion_*.pdb' | wc -l)"
    existing_trbs="$(find "${raw_dir}" -maxdepth 1 -type f -name 'rfdiffusion_*.trb' | wc -l)"
    if [[ "${existing_pdbs}" -ne "${existing_trbs}" ]]; then
        echo "Mismatched direct scaffold outputs: ${existing_pdbs} PDB, ${existing_trbs} TRB"
        exit 1
    fi
    if [[ "${existing_pdbs}" -gt "${num_designs}" ]]; then
        echo "Direct scaffold output count exceeds requested designs"
        exit 1
    fi

    remaining_designs="$((num_designs - existing_pdbs))"
    if [[ "${remaining_designs}" -gt 0 ]]; then
        {
            echo -e "stage\tstatus\tstart_utc\tdesign_startnum\tnum_designs"
            echo -e "direct_scaffold_rfdiffusion\tRUNNING\t$(date -u +%FT%TZ)\t${existing_pdbs}\t${remaining_designs}"
        } > "${backbone_dir}/direct_trace.tsv"
        set +e
        docker run --rm \
            --gpus all \
            -v "${ovo_home_dir}:${ovo_home_dir}" \
            -w "${backbone_dir}" \
            ovo-rfdiffusion \
            bash -lc "
                set -euo pipefail
                export PYTHONPATH=/opt/RFdiffusion:/opt/RFdiffusion/env/SE3Transformer
                export HYDRA_FULL_ERROR=1
                python3 /opt/RFdiffusion/scripts/run_inference.py \
                    inference.output_prefix=${raw_dir}/rfdiffusion \
                    inference.model_directory_path=${model_dir} \
                    inference.schedule_directory_path=${backbone_dir}/schedules \
                    inference.num_designs=${remaining_designs} \
                    inference.design_startnum=${existing_pdbs} \
                    inference.write_trajectory=false \
                    ppi.hotspot_res=[${generation_hotspots}] \
                    diffuser.T=50 \
                    scaffoldguided.scaffoldguided_enable=True \
                    scaffoldguided.target_pdb=True \
                    scaffoldguided.target_path=${input_pdb} \
                    scaffoldguided.target_ss=${target_ss} \
                    scaffoldguided.target_adj=${target_adj} \
                    scaffoldguided.scaffold_dir=${scaffold_dir} \
                    scaffoldguided.mask_loops=False \
                    denoiser.noise_scale_ca=0 \
                    denoiser.noise_scale_frame=0
            " > "${backbone_dir}/direct_rfdiffusion.stdout.log" 2>&1
        direct_exit="$?"
        set -e
        echo -e "direct_scaffold_rfdiffusion\t$([[ "${direct_exit}" -eq 0 ]] && echo COMPLETED || echo FAILED)\t$(date -u +%FT%TZ)\t${existing_pdbs}\t${remaining_designs}" \
            >> "${backbone_dir}/direct_trace.tsv"
        if [[ "${direct_exit}" -ne 0 ]]; then
            exit "${direct_exit}"
        fi
    fi

    observed_raw_pdbs="$(find "${raw_dir}" -maxdepth 1 -type f -name 'rfdiffusion_*.pdb' | wc -l)"
    observed_raw_trbs="$(find "${raw_dir}" -maxdepth 1 -type f -name 'rfdiffusion_*.trb' | wc -l)"
    if [[ "${observed_raw_pdbs}" -ne "${num_designs}" || "${observed_raw_trbs}" -ne "${num_designs}" ]]; then
        echo "Expected ${num_designs} paired direct scaffold outputs; found ${observed_raw_pdbs} PDB and ${observed_raw_trbs} TRB"
        exit 1
    fi

    cp "${raw_dir}"/*.pdb "${pdb_dir}/"
    cp "${raw_dir}"/*.trb "${raw_trb_dir}/"
    "${python_bin}" "${campaign_dir}/scripts/normalize_scaffold_trb.py" \
        "${raw_trb_dir}" \
        "${trb_dir}" \
        --report "${condition_dir}/reports/scaffold_trb_normalization.json"
    "${python_bin}" \
        "${ovo_package_dir}/pipelines/rfdiffusion-backbone/bin/standardize_pdb.py" \
        "${pdb_dir}" \
        "${trb_dir}" \
        "${standardized_dir}"
    touch "${backbone_dir}/direct_scaffold.completed"
fi

observed_backbones="$(
    find -L "${standardized_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l
)"
if [[ "${observed_backbones}" -ne "${num_designs}" ]]; then
    echo "Expected ${num_designs} standardized backbones for ${condition_id}; found ${observed_backbones}"
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
            -work-dir "${phase_dir}/work/${condition_id}/metrics" \
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

metrics_csv="${metrics_dir}/output/${condition_id}/backbone_metrics.csv"
"${python_bin}" "${campaign_dir}/scripts/filter_r2_phase_b_metrics.py" \
    "${metrics_csv}" \
    --pdb-dir "${standardized_dir}" \
    --condition "${condition_id}" \
    --filtered-csv "${condition_dir}/reports/backbone_metrics_r2.csv" \
    --selected-dir "${condition_dir}/selected_backbones" \
    --summary-json "${condition_dir}/reports/phase_b_backbone_summary.json"

echo "R2 Phase B condition ${condition_id} completed"
