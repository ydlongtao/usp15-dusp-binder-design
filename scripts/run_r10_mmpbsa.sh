#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 PREPARED_DIR REPLICA_DIR OUTPUT_DIR" >&2
  exit 2
fi

prepared_dir=$1
replica_dir=$2
output_dir=$3
mkdir -p "$output_dir"
# MMPBSA.py creates intermediate files in the current working directory.
# Keep those files on the GPFS work allocation rather than the 30-GB home NFS.
cd "$output_dir"

if [[ -s "$output_dir/FINAL_RESULTS_MMPBSA.dat" ]]; then
  exit 0
fi

input_file="$output_dir/mmpbsa.in"
{
  echo "Relative MM/GBSA diagnostic for USP15 R10"
  echo "&general"
  # Protein coordinates are written every 10 ps.  Start after the
  # predeclared 20 ns burn-in and sample every 1 ns thereafter.
  echo "  startframe=2001, endframe=99999999, interval=100,"
  echo "  receptor_mask=':77-205', ligand_mask=':1-76',"
  echo "  keep_files=0, verbose=1,"
  echo "/"
  echo "&gb"
  echo "  igb=8, saltcon=0.150,"
  echo "/"
} >"$input_file"

MMPBSA.py -O \
  -i "$input_file" \
  -o "$output_dir/FINAL_RESULTS_MMPBSA.dat" \
  -eo "$output_dir/FINAL_RESULTS_MMPBSA.csv" \
  -cp "$prepared_dir/complex.prmtop" \
  -rp "$prepared_dir/target.prmtop" \
  -lp "$prepared_dir/binder.prmtop" \
  -y "$replica_dir/production_protein.xtc" \
  >"$output_dir/mmpbsa.stdout.log" \
  2>"$output_dir/mmpbsa.stderr.log"
