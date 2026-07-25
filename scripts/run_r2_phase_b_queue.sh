#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
phase_dir="${R2_PHASE_B_DIR:-${campaign_dir}/r2/phase_b}"
phase_matrix="${R2_PHASE_B_MATRIX:-${campaign_dir}/config/r2_phase_b.tsv}"
queue_log="${phase_dir}/queue_status.tsv"

mkdir -p "${phase_dir}"
exec 9>"${phase_dir}/phase_b_queue.lock"
if ! flock -n 9; then
    echo "Another R2 Phase B queue already holds ${phase_dir}/phase_b_queue.lock"
    exit 1
fi

if [[ ! -s "${queue_log}" ]]; then
    echo -e "condition_id\tstatus\ttimestamp_utc" > "${queue_log}"
fi

mapfile -t condition_ids < <(
    awk -F '\t' 'NR > 1 && $1 != "" {print $1}' "${phase_matrix}"
)
if [[ "${#condition_ids[@]}" -ne 6 ]]; then
    echo "Expected six R2 Phase B conditions; found ${#condition_ids[@]}"
    exit 1
fi

for condition_id in "${condition_ids[@]}"; do
    summary="${phase_dir}/conditions/${condition_id}/reports/phase_b_backbone_summary.json"
    if [[ -s "${summary}" ]] && \
        "${PYTHON_BIN:-/usr/bin/python3}" - "${summary}" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1], encoding="utf-8"))
raise SystemExit(0 if data.get("generated_backbones") == 50 else 1)
PY
    then
        echo -e "${condition_id}\tALREADY_COMPLETED\t$(date -u +%FT%TZ)" >> "${queue_log}"
        continue
    fi

    echo -e "${condition_id}\tRUNNING\t$(date -u +%FT%TZ)" >> "${queue_log}"
    if "${campaign_dir}/scripts/run_r2_phase_b_backbones.sh" "${condition_id}"; then
        echo -e "${condition_id}\tCOMPLETED\t$(date -u +%FT%TZ)" >> "${queue_log}"
    else
        queue_exit="$?"
        echo -e "${condition_id}\tFAILED_${queue_exit}\t$(date -u +%FT%TZ)" >> "${queue_log}"
        exit "${queue_exit}"
    fi
done

touch "${phase_dir}/phase_b_backbones.completed"
echo "R2 Phase B backbone queue completed"
