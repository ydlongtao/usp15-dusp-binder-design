#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_dir="${R4_POSE_ENSEMBLE_DIR:-${campaign_dir}/r4/6dj9_pose_ensemble}"
source_complex="${R4_6DJ9_PDB:-${campaign_dir}/inputs/6DJ9.pdb}"
target_reference="${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb"
input_dir="${phase_dir}/input"
sequence_work="${phase_dir}/sequence_work"
sequence_dir="${phase_dir}/sequence"
sequence_pdb_dir="${sequence_dir}/standardized_pdb"
af2_dir="${phase_dir}/af2"
report_dir="${phase_dir}/reports"
af2_test="af2_model_1_multimer_tt_3rec"
jsonl="${af2_dir}/output/contig1_batch1/${af2_test}.jsonl"

mkdir -p \
    "${input_dir}" \
    "${sequence_work}/input" \
    "${sequence_dir}" \
    "${af2_dir}" \
    "${report_dir}" \
    "${phase_dir}/work/af2"
exec 9>"${phase_dir}/r4_pose_ensemble.lock"
if ! flock -n 9; then
    echo "Another R4 crystal-pose ensemble driver holds the lock"
    exit 1
fi

if [[ ! -f "${report_dir}/preparation.json" ]]; then
    "${python_bin}" "${campaign_dir}/scripts/prepare_r4_6dj9_pose_ensemble.py" \
        --complex "${source_complex}" \
        --target-reference "${target_reference}" \
        --output-dir "${input_dir}" \
        --report "${report_dir}/preparation.json"
fi
if [[ "$(find "${input_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne 4 ]]; then
    echo "Expected four validated 6DJ9 pose inputs"
    exit 1
fi
"${python_bin}" "${campaign_dir}/scripts/validate_backbones.py" \
    "${input_dir}" \
    --json "${report_dir}/backbone_gate.json" \
    --csv "${report_dir}/backbone_gate.csv" \
    --require-all-pass

if [[ ! -f "${sequence_dir}/proteinmpnn.completed" ]]; then
    if find "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' \
        -print -quit 2>/dev/null | grep -q .; then
        echo "Unmarked R4 ProteinMPNN outputs exist; preserve and audit before resume"
        exit 1
    fi
    cp "${input_dir}"/*.pdb "${sequence_work}/input/"
    (
        cd "${sequence_work}"
        "${python_bin}" \
            "${ovo_package_dir}/pipelines/ligandmpnn-sequence-design/bin/prepare_json.py" \
            --pdb_dir input \
            --pdb_ids_json pdb_ids.json \
            --redesigned_residues_json redesigned_residues_multi.json \
            --remark_json remark_multi.json
        docker run --rm --gpus all \
            --user "$(id -u):$(id -g)" \
            -v "${sequence_work}:/work" \
            -w /work \
            ovo-ligandmpnn \
            bash -lc '
                ln -s /opt/LigandMPNN/model_params ./model_params
                python /opt/LigandMPNN/run.py \
                    --model_type protein_mpnn \
                    --pdb_path_multi pdb_ids.json \
                    --redesigned_residues_multi redesigned_residues_multi.json \
                    --out_folder output \
                    --number_of_batches 3 \
                    --pack_side_chains 1 \
                    --number_of_packs_per_design 1 \
                    --repack_everything 1 \
                    --temperature 0.1 \
                    --omit_AA C
            ' > "${sequence_dir}/proteinmpnn.stdout.log" 2>&1
        bash \
            "${ovo_package_dir}/pipelines/ligandmpnn-sequence-design/bin/copy_remarks.sh" \
            remark_multi.json \
            output/packed/ \
            "${sequence_pdb_dir}/"
    )
    "${python_bin}" "${campaign_dir}/scripts/validate_r4_pose_sequences.py" \
        --input-dir "${input_dir}" \
        --sequence-pdb-dir "${sequence_pdb_dir}" \
        --report "${report_dir}/sequence_validation.json"
    touch "${sequence_dir}/proteinmpnn.completed"
fi
if [[ "$(find "${sequence_pdb_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)" -ne 12 ]]; then
    echo "Expected twelve fixed-interface ProteinMPNN outputs"
    exit 1
fi

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
            --input_designs "${sequence_pdb_dir}/" \
            --native_pdb "${target_reference}" \
            --tests "${af2_test}" \
            --design_type binder \
            --batch_size 20 \
            -resume \
            > "${af2_dir}/nextflow.stdout.log" 2>&1
    )
fi
if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 12 ]]; then
    echo "Expected twelve R4 crystal-pose AF2 records"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/summarize_r3_partial_diffusion.py" \
    --jsonl "${jsonl}" \
    --sequence-pdb-dir "${sequence_pdb_dir}" \
    --selection-report "${report_dir}/preparation.json" \
    --sequence-model protein_mpnn \
    --phase-label "R4 6DJ9 crystal-pose ensemble" \
    --csv-output "${report_dir}/af2_metrics.csv" \
    --json-output "${report_dir}/summary.json" \
    --fasta-output "${report_dir}/sequences.fasta"

touch "${phase_dir}/r4_pose_ensemble.completed"
echo "R4 6DJ9 crystal-pose ensemble completed"
