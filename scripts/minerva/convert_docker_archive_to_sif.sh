#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 OPENMM_DOCKER_TAR_GZ OUTPUT_SIF" >&2
  exit 2
fi

archive=$1
output_sif=$2
runtime=${APPTAINER_CMD:-apptainer}

if [[ ! -s "$archive" ]]; then
  echo "Docker archive not found: $archive" >&2
  exit 3
fi
if ! command -v "$runtime" >/dev/null 2>&1 &&
   type module >/dev/null 2>&1; then
  module load apptainer/1.4.5
fi
if ! command -v "$runtime" >/dev/null 2>&1; then
  if command -v singularity >/dev/null 2>&1; then
    runtime=singularity
  else
    echo "Neither apptainer nor singularity is available" >&2
    exit 4
  fi
fi

mkdir -p "$(dirname "$output_sif")"
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-"$(dirname "$output_sif")/.apptainer-tmp"}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-"$(dirname "$output_sif")/.apptainer-cache"}
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

archive_input=$archive
temporary_archive=
archive_magic=$(od -An -tx1 -N2 "$archive" | tr -d '[:space:]')
if [[ "$archive_magic" == "1f8b" ]]; then
  temporary_archive=$(mktemp "$APPTAINER_TMPDIR/openmm-archive.XXXXXX.tar")
  trap 'rm -f "$temporary_archive"' EXIT
  echo "Expanding gzip archive into scratch: $temporary_archive"
  gzip -dc "$archive" >"$temporary_archive"
  archive_input=$temporary_archive
fi

archive_type=${APPTAINER_ARCHIVE_TYPE:-}
if [[ -z "$archive_type" ]]; then
  if tar -tf "$archive_input" | grep -qx 'oci-layout'; then
    archive_type=oci-archive
  elif tar -tf "$archive_input" | grep -qx 'manifest.json'; then
    archive_type=docker-archive
  else
    echo "Archive is neither an OCI layout nor a Docker archive: $archive" >&2
    exit 5
  fi
fi
case "$archive_type" in
  oci-archive|docker-archive) ;;
  *)
    echo "Unsupported APPTAINER_ARCHIVE_TYPE: $archive_type" >&2
    exit 6
    ;;
esac

echo "Building from ${archive_type}:$archive_input"
"$runtime" build "$output_sif" "${archive_type}:$archive_input"
"$runtime" exec "$output_sif" /opt/conda/bin/python -c \
  'import openmm; print(openmm.__version__)'
sha256sum "$output_sif" >"$output_sif.sha256"
