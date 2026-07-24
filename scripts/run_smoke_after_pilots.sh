#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
pilot_queue_session="${USP15_PILOT_QUEUE_SESSION:-usp15_pilot_queue}"
pools=(
    USP15_R1_short_base
    USP15_R1_short_beta
    USP15_R1_long_base
    USP15_R1_long_beta
)

while tmux has-session -t "${pilot_queue_session}" 2>/dev/null; do
    sleep 30
done

for pool in "${pools[@]}"; do
    if [[ ! -s "${campaign_dir}/reports/${pool}.pilot_summary.json" ]]; then
        echo "Missing pilot summary: ${pool}"
        exit 1
    fi
done

source_pool=""
for pool in "${pools[@]}"; do
    filtered_dir="${campaign_dir}/pilots/${pool}/metrics/output/${pool}/backbones_filtered"
    if find -L "${filtered_dir}" -maxdepth 1 -type f -name '*.pdb' -print -quit 2>/dev/null | grep -q .; then
        source_pool="${pool}"
        break
    fi
done

if [[ -z "${source_pool}" ]]; then
    echo "No pilot backbone passed the strict hard filters"
    exit 1
fi

echo "Using ${source_pool} for LigandMPNN + AF2 smoke"
"${campaign_dir}/scripts/run_sequence_af2_smoke.sh" "${source_pool}"
