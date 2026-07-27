#!/usr/bin/env bash
set -euo pipefail

: "${LSF_PROJECT:?Set LSF_PROJECT, for example acc_yourlab}"
: "${USP15_MD_DIR:?Set an absolute Minerva path to md_openmm}"
: "${OPENMM_SIF:?Set an absolute Minerva path to the OpenMM SIF image}"

if ! grep -q '"status": "passed"' "$USP15_MD_DIR/smoke_minerva/audit.json"; then
  echo "Run submit_smoke.sh and wait for smoke_minerva/audit.json to pass" >&2
  exit 5
fi

queue=${MINERVA_GPU_QUEUE:-gpu}
gpu_model=${MINERVA_GPU_MODEL:-a100}
walltime=${MINERVA_REPLICA_WALLTIME:-36:00}
array_range=${MINERVA_ARRAY_RANGE:-1-30}
maximum_concurrent=${MINERVA_MAX_CONCURRENT:-1}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

case "$gpu_model" in
  v100|a100|a10080g|h10080g|h100nvl|l40s|b200) ;;
  *)
    echo "Unsupported Minerva GPU resource name: $gpu_model" >&2
    exit 3
    ;;
esac
if [[ "$maximum_concurrent" != "1" ]]; then
  echo "This campaign requires MINERVA_MAX_CONCURRENT=1" >&2
  exit 4
fi

mkdir -p "$USP15_MD_DIR/lsf_logs"
export USP15_MD_DIR OPENMM_SIF
bsub \
  -J "usp15_r10_md[${array_range}]%${maximum_concurrent}" \
  -P "$LSF_PROJECT" \
  -q "$queue" \
  -n 4 \
  -R "span[hosts=1]" \
  -R "rusage[mem=8000]" \
  -R "$gpu_model" \
  -gpu "num=1" \
  -W "$walltime" \
  -oo "$USP15_MD_DIR/lsf_logs/replica.%J.%I.out" \
  -eo "$USP15_MD_DIR/lsf_logs/replica.%J.%I.err" \
  <"$script_dir/run_replica_job.lsf"
