#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR to the OVO conda environment}"
ovo_package_dir="${OVO_PACKAGE_DIR:-${ovo_env_dir}/lib/python3.13/site-packages/ovo}"
nextflow_bin="${NEXTFLOW_BIN:-${ovo_env_dir}/bin/nextflow}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
phase_dir="${R3_NATIVE_DIAGNOSTIC_DIR:-${campaign_dir}/r3/native_6dj9_diagnostic}"
preparation_dir="${phase_dir}/preparation"
input_pdb="${phase_dir}/input/ubv15d_native_wt_diagnostic.pdb"
af2_dir="${phase_dir}/af2"
af2_test="af2_model_1_multimer_tt_3rec"
jsonl="${af2_dir}/output/contig1_batch1/${af2_test}.jsonl"

mkdir -p \
    "${preparation_dir}/controls" \
    "${phase_dir}/input" \
    "${af2_dir}" \
    "${phase_dir}/work/af2" \
    "${phase_dir}/reports"
exec 9>"${phase_dir}/r3_native_diagnostic.lock"
if ! flock -n 9; then
    echo "Another R3 native diagnostic driver holds the lock"
    exit 1
fi

"${python_bin}" "${campaign_dir}/scripts/prepare_r3_ubv_controls.py" \
    --source "${campaign_dir}/inputs/6DJ9.pdb" \
    --target-reference "${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb" \
    --output-dir "${preparation_dir}/controls" \
    --report "${preparation_dir}/preparation.json" \
    --fasta "${preparation_dir}/controls.fasta" \
    --native-diagnostic-output "${input_pdb}"

if [[ ! -s "${af2_dir}/trace.txt" ]] || ! awk -F '\t' \
    'NR > 1 && $5 == "COMPLETED" && $6 == 0 {ok=1} END {exit !ok}' \
    "${af2_dir}/trace.txt"; then
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
            --input_designs "${input_pdb}" \
            --tests "${af2_test}" \
            --design_type binder \
            --batch_size 20 \
            -resume \
            > "${af2_dir}/nextflow.stdout.log" 2>&1
    )
fi

if [[ ! -s "${jsonl}" ]] || [[ "$(wc -l < "${jsonl}")" -ne 1 ]]; then
    echo "Expected one native diagnostic AF2 record"
    exit 1
fi

"${python_bin}" - "${jsonl}" "${phase_dir}/reports/summary.json" <<'PY'
import json
import sys
from pathlib import Path

record = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8").strip())
record["passed_all_af2_gates"] = (
    float(record["ipae"]) <= 10.0
    and float(record["target_aligned_binder_rmsd"]) <= 2.0
    and float(record["binder_plddt"]) >= 80.0
)
summary = {
    "phase": "R3 exact 6DJ9 native-complex diagnostic",
    "candidate_eligible": False,
    "thresholds_unchanged": {
        "ipae_max": 10.0,
        "target_aligned_binder_rmsd_A_max": 2.0,
        "binder_plddt_min": 80.0,
    },
    "record": record,
    "diagnosis": (
        "3t9l_pose_transfer_requires_optimization"
        if record["passed_all_af2_gates"]
        else "tt_protocol_does_not_recover_native_ubv_sequence"
    ),
}
Path(sys.argv[2]).write_text(
    json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
    encoding="utf-8",
)
print(json.dumps(summary, ensure_ascii=False))
PY

touch "${phase_dir}/r3_native_diagnostic.completed"
echo "R3 exact-native diagnostic completed"
