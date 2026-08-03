#!/usr/bin/env bash
set -euo pipefail

: "${LOCAL_MD_DIR:?Set LOCAL_MD_DIR to the local md_openmm directory}"
: "${MINERVA_SSH:?Set MINERVA_SSH to a configured SSH host or user@host}"
: "${MINERVA_DEST:?Set MINERVA_DEST to an absolute Minerva destination directory}"

if [[ ! -d "$LOCAL_MD_DIR/prepared" ]] ||
   [[ ! -s "$LOCAL_MD_DIR/reports/prepared_system_audit.json" ]]; then
  echo "LOCAL_MD_DIR is not a complete audited MD directory" >&2
  exit 3
fi

ssh "$MINERVA_SSH" "mkdir -p '$MINERVA_DEST/md_openmm'"
rsync -a --partial --progress \
  --exclude '.queue.lock' \
  "$LOCAL_MD_DIR/" \
  "$MINERVA_SSH:$MINERVA_DEST/md_openmm/"

echo "Upload complete. Verify transfer_manifest.json on Minerva before use."
