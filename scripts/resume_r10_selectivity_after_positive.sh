#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
positive_session="${R10_POSITIVE_SESSION:-usp15-r10-positive}"
summary="${campaign_dir}/r10/positive_screen/reports/summary.json"
poll_seconds="${R10_POLL_SECONDS:-20}"

while tmux has-session -t "${positive_session}" 2>/dev/null; do
    sleep "${poll_seconds}"
done

if [[ ! -f "${summary}" ]]; then
    echo "R10 positive screen did not produce a summary"
    exit 1
fi

"${campaign_dir}/scripts/run_r10_selectivity_screen.sh"
echo "R10 positive and selectivity screens completed"
