#!/usr/bin/env bash
set -euo pipefail

: "${USP15_MD_DIR:?Set USP15_MD_DIR to the R10 md_openmm directory}"

image=${OPENMM_IMAGE:-usp15-openmm:8.5.2}
allow_shared_gpu=${USP15_ALLOW_SHARED_GPU:-0}
minimum_free_gpu_mb=${USP15_MIN_FREE_GPU_MB:-20000}
scripts_dir="$USP15_MD_DIR/scripts"
inputs_dir="$USP15_MD_DIR/inputs"
prepared_root="$USP15_MD_DIR/prepared"
runs_root="$USP15_MD_DIR/runs"
analysis_root="$USP15_MD_DIR/analysis"
reports_dir="$USP15_MD_DIR/reports"
logs_dir="$USP15_MD_DIR/logs"
mkdir -p "$prepared_root" "$runs_root" "$analysis_root" "$reports_dir" "$logs_dir"

exec 9>"$USP15_MD_DIR/.queue.lock"
if ! flock -n 9; then
  echo "Another USP15 MD queue already holds $USP15_MD_DIR/.queue.lock" >&2
  exit 3
fi

container_base=(
  docker run --rm
  --user "$(id -u):$(id -g)"
  --shm-size=8g
  -v "$USP15_MD_DIR:/campaign"
  "$image"
)

container_gpu=(
  docker run --rm
  --gpus "device=0"
  --user "$(id -u):$(id -g)"
  --shm-size=8g
  -v "$USP15_MD_DIR:/campaign"
  "$image"
)

wait_for_free_gpu() {
  while true; do
    mapfile -t gpu_pids < <(
      nvidia-smi \
        --query-compute-apps=pid \
        --format=csv,noheader,nounits 2>/dev/null |
        awk 'NF {print $1}'
    )
    if [[ ${#gpu_pids[@]} -eq 0 ]]; then
      return 0
    fi
    if [[ "$allow_shared_gpu" == "1" ]]; then
      free_gpu_mb=$(
        nvidia-smi \
          --query-gpu=memory.free \
          --format=csv,noheader,nounits 2>/dev/null |
          awk 'NR == 1 {print int($1)}'
      )
      if [[ "$free_gpu_mb" =~ ^[0-9]+$ ]] &&
         (( free_gpu_mb >= minimum_free_gpu_mb )); then
        {
          printf '%s\tshared_gpu_authorized\tfree_mb=%s\tpids=' \
            "$(date -Is)" "$free_gpu_mb"
          printf '%s,' "${gpu_pids[@]}"
          echo
        } >>"$logs_dir/gpu_wait.log"
        return 0
      fi
    fi
    {
      printf '%s\twaiting_for_gpu\t' "$(date -Is)"
      printf '%s,' "${gpu_pids[@]}"
      echo
    } >>"$logs_dir/gpu_wait.log"
    sleep 60
  done
}

prepare_rank() {
  local rank=$1
  local input="$inputs_dir/USP15_R10_rank${rank}.pdb"
  local output="$prepared_root/rank${rank}"
  if [[ -s "$output/preparation.json" ]]; then
    return 0
  fi
  mkdir -p "$output"
  "${container_base[@]}" \
    python /campaign/scripts/prepare_r10_openmm_system.py \
      --input-pdb "/campaign/inputs/USP15_R10_rank${rank}.pdb" \
      --output-dir "/campaign/prepared/rank${rank}" \
      >"$logs_dir/prepare_rank${rank}.log" 2>&1
}

run_replica() {
  local rank=$1
  local seed=$2
  local output="$runs_root/rank${rank}/seed${seed}"
  mkdir -p "$output"
  if [[ -s "$output/status.json" ]] &&
     grep -q '"status": "completed"' "$output/status.json"; then
    return 0
  fi
  wait_for_free_gpu
  "${container_gpu[@]}" \
    python /campaign/scripts/run_r10_openmm_replica.py \
      --prepared-dir "/campaign/prepared/rank${rank}" \
      --output-dir "/campaign/runs/rank${rank}/seed${seed}" \
      --seed "$((100000 + 1000 * 10#$rank + seed))" \
      --production-ns 100 \
      >"$logs_dir/run_rank${rank}_seed${seed}.log" 2>&1
}

analyze_replica() {
  local rank=$1
  local seed=$2
  local output="$analysis_root/rank${rank}/seed${seed}"
  mkdir -p "$output"
  if [[ ! -s "$output/summary.json" ]]; then
    "${container_base[@]}" \
      python /campaign/scripts/analyze_r10_openmm_replica.py \
        --topology "/campaign/prepared/rank${rank}/protein_protonated.pdb" \
        --trajectory "/campaign/runs/rank${rank}/seed${seed}/production_protein.xtc" \
        --output-dir "/campaign/analysis/rank${rank}/seed${seed}" \
        >"$logs_dir/analyze_rank${rank}_seed${seed}.log" 2>&1
  fi
  if [[ ! -s "$output/mmpbsa/FINAL_RESULTS_MMPBSA.dat" ]]; then
    "${container_base[@]}" \
      bash /campaign/scripts/run_r10_mmpbsa.sh \
        "/campaign/prepared/rank${rank}" \
        "/campaign/runs/rank${rank}/seed${seed}" \
        "/campaign/analysis/rank${rank}/seed${seed}/mmpbsa"
  fi
}

for rank in $(seq -w 1 10); do
  prepare_rank "$rank"
done

if [[ ! -s "$reports_dir/prepared_system_audit.json" ]]; then
  "${container_base[@]}" \
    python /campaign/scripts/audit_r10_openmm_prepared.py \
      --prepared-root /campaign/prepared \
      --input-root /campaign/inputs \
      --output /campaign/reports/prepared_system_audit.json \
      >"$logs_dir/prepared_system_audit.log" 2>&1
fi

# A real CUDA smoke must complete before the formal 3 x 100 ns queue begins.
smoke_dir="$USP15_MD_DIR/smoke/rank01_seed0"
if [[ ! -s "$smoke_dir/status.json" ]] ||
   ! grep -q '"status": "completed"' "$smoke_dir/status.json"; then
  mkdir -p "$smoke_dir"
  wait_for_free_gpu
  "${container_gpu[@]}" \
    python /campaign/scripts/run_r10_openmm_replica.py \
      --prepared-dir /campaign/prepared/rank01 \
      --output-dir /campaign/smoke/rank01_seed0 \
      --seed 101000 \
      --smoke \
      >"$logs_dir/smoke_rank01_seed0.log" 2>&1
fi

if [[ ! -s "$USP15_MD_DIR/smoke/analysis/summary.json" ]]; then
  mkdir -p "$USP15_MD_DIR/smoke/analysis"
  "${container_base[@]}" \
    python /campaign/scripts/analyze_r10_openmm_replica.py \
      --topology /campaign/prepared/rank01/protein_protonated.pdb \
      --trajectory /campaign/smoke/rank01_seed0/production_protein.xtc \
      --output-dir /campaign/smoke/analysis \
      --burn-in-ns 0 \
      >"$logs_dir/smoke_analysis.log" 2>&1
fi

if [[ ! -s "$USP15_MD_DIR/smoke/audit.json" ]]; then
  "${container_base[@]}" \
    python /campaign/scripts/audit_r10_openmm_smoke.py \
      --status /campaign/smoke/rank01_seed0/status.json \
      --analysis /campaign/smoke/analysis/summary.json \
      --preparation /campaign/prepared/rank01/preparation.json \
      --output /campaign/smoke/audit.json \
      >"$logs_dir/smoke_audit.log" 2>&1
fi

for rank in $(seq -w 1 10); do
  for seed in 0 1 2; do
    run_replica "$rank" "$seed"
    analyze_replica "$rank" "$seed"
  done
done

touch "$reports_dir/r10_3x100ns.completed"
