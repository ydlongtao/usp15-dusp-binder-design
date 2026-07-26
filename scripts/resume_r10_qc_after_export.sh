#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
selectivity_session="${R10_SELECTIVITY_SESSION:-usp15-r10-selectivity-watch}"
download_session="${R10_ESM_DOWNLOAD_SESSION:-usp15-r10-esm-download}"
poll_seconds="${R10_POLL_SECONDS:-20}"

while tmux has-session -t "${selectivity_session}" 2>/dev/null \
    || tmux has-session -t "${download_session}" 2>/dev/null
do
    sleep "${poll_seconds}"
done

if [[ ! -f "${campaign_dir}/r10/r10_final_export.completed" ]]; then
    echo "R10 final export is not complete"
    exit 1
fi

"${campaign_dir}/scripts/run_r10_protein_qc.sh"
echo "R10 final export and ProteinQC completed"
