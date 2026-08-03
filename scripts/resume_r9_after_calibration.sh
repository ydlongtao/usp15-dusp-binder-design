#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
calibration_session="${R9_CALIBRATION_SESSION:-usp15-r9-calibration}"
summary="${campaign_dir}/r9/tt_calibration/reports/summary.json"
poll_seconds="${R9_POLL_SECONDS:-20}"

while tmux has-session -t "${calibration_session}" 2>/dev/null; do
    sleep "${poll_seconds}"
done

if [[ ! -f "${summary}" ]] || ! "${python_bin}" - "${summary}" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1]))
raise SystemExit(0 if summary["calibrated"] is True else 1)
PY
then
    echo "R9 calibration did not pass; positive screen not started"
    exit 1
fi

"${campaign_dir}/scripts/run_r9_candidate_screen.sh"
echo "R9 calibrated positive screen completed"
