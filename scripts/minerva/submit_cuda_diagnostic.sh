#!/usr/bin/env bash
set -euo pipefail

: "${LSF_PROJECT:?Set LSF_PROJECT, for example acc_yourlab}"
: "${USP15_MD_DIR:?Set an absolute Minerva path to md_openmm}"
: "${OPENMM_SIF:?Set an absolute Minerva path to the OpenMM SIF image}"

queue=${MINERVA_GPU_QUEUE:-gpuexpress}
gpu_model=${MINERVA_GPU_MODEL:-h100nvl}
walltime=${MINERVA_DIAGNOSTIC_WALLTIME:-00:30}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

case "$gpu_model" in
  v100|a100|a10080g|h10080g|h100nvl|l40s|b200) ;;
  *)
    echo "Unsupported Minerva GPU resource name: $gpu_model" >&2
    exit 3
    ;;
esac

mkdir -p "$USP15_MD_DIR/lsf_logs"
export USP15_MD_DIR OPENMM_SIF
bsub \
  -J usp15_r10_cuda_diagnostic \
  -P "$LSF_PROJECT" \
  -q "$queue" \
  -n 1 \
  -R "select[$gpu_model] span[hosts=1] rusage[mem=4000]" \
  -gpu "num=1" \
  -W "$walltime" \
  -oo "$USP15_MD_DIR/lsf_logs/cuda_diagnostic.%J.out" \
  -eo "$USP15_MD_DIR/lsf_logs/cuda_diagnostic.%J.err" \
  <"$script_dir/run_cuda_diagnostic_job.lsf"
