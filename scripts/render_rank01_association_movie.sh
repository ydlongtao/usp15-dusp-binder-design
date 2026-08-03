#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?association result directory}
OUT=${2:-$ROOT/animation}
mkdir -p "$OUT"
cp "$(dirname "$0")/make_rank01_association_movie.pml" "$OUT/movie.pml"
cd "$ROOT"
pymol -cq "$OUT/movie.pml"
if command -v ffmpeg >/dev/null 2>&1; then
  ffmpeg -y -loglevel error -framerate 12 -i "$OUT/association_frames%04d.png" -vf "format=yuv420p" "$OUT/rank01_association.mp4"
  ffmpeg -y -loglevel error -i "$OUT/rank01_association.mp4" -vf "fps=8,scale=960:-2:flags=lanczos" "$OUT/rank01_association.gif"
fi
