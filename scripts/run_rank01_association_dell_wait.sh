#!/usr/bin/env bash
set -euo pipefail
BASE=/DATABANK/users/hflt/usp15_md_simulation
PREP=$BASE/inputs/rank01_association_prepared
SCRIPT=$BASE/scripts/run_rank01_association.py
OUT=$BASE/results/rank01_association_dell
PY=$BASE/envs/md/bin/python
mkdir -p "$OUT" "$BASE/logs"
echo "$(date -Is) waiting for Dell V100 free memory >= 20000 MiB" >> "$BASE/logs/rank01_association_wait.log"
while true; do
  free_mb=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ' || true)
  pids=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | sed '/^$/d' || true)
  if [ -n "$free_mb" ] && [ "$free_mb" -ge 20000 ]; then break; fi
  echo "$(date -Is) waiting; free=${free_mb:-unknown} MiB; other_pids=${pids:-none}" >> "$BASE/logs/rank01_association_wait.log"
  sleep 60
done
echo "$(date -Is) free memory threshold met; starting rank01 association seeds 0-2 serially" >> "$BASE/logs/rank01_association_wait.log"
for seed in 0 1 2; do
  "$PY" "$SCRIPT" --prepared-dir "$PREP" --output-dir "$OUT/seed$seed" --seed $((210001+seed)) --production-ns 100 --separation-a 35 --platform OpenCL >> "$BASE/logs/rank01_association_seed${seed}.log" 2>&1
done
echo "$(date -Is) rank01 association seeds 0-2 completed" >> "$BASE/logs/rank01_association_wait.log"
