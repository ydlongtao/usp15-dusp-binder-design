#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
active_session="${USP15_ACTIVE_PILOT_SESSION:-usp15_pilot_short_base}"

while tmux has-session -t "${active_session}" 2>/dev/null; do
    sleep 30
done

if [[ ! -s "${campaign_dir}/reports/USP15_R1_short_base.pilot_summary.json" ]]; then
    echo "short_base pilot did not produce its validated summary"
    exit 1
fi

for pool in \
    USP15_R1_short_beta \
    USP15_R1_long_base \
    USP15_R1_long_beta
do
    echo "Starting ${pool}"
    "${campaign_dir}/scripts/run_backbone_pilot.sh" "${pool}"
done

echo "All four backbone pilots completed"
