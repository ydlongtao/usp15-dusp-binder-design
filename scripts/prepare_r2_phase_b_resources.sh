#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR to the campaign working directory}"
phase_dir="${R2_PHASE_B_DIR:-${campaign_dir}/r2/phase_b}"
image="${RFDIFFUSION_IMAGE:-ovo-rfdiffusion}"
resource_dir="${phase_dir}/resources"
target_pdb="${campaign_dir}/inputs/USP15_DUSP_3T9L_A6-134.pdb"

mkdir -p "${resource_dir}"

if docker run --rm "${image}" python3 -c "import pyrosetta" >/dev/null 2>&1; then
    echo "PyRosetta is installed in ${image}; refusing to use the helper implicitly"
    exit 1
fi

docker run --rm -i \
    -v "${campaign_dir}:${campaign_dir}" \
    "${image}" \
    bash -lc "
        set -euo pipefail
        mkdir -p '${resource_dir}/target_folds'
        python3 /opt/RFdiffusion/helper_scripts/make_secstruc_adj.py \
            --input_pdb '${target_pdb}' \
            --out_dir '${resource_dir}/target_folds'
        mkdir -p '${resource_dir}/all_scaffolds'
        tar -xzf /opt/RFdiffusion/examples/ppi_scaffolds_subset.tar.gz \
            -C '${resource_dir}/all_scaffolds'
    "

docker run --rm -i \
    -v "${campaign_dir}:${campaign_dir}" \
    "${image}" \
    python3 - "${resource_dir}" <<'PY'
from pathlib import Path
import csv
import shutil
import sys
import torch

resource_dir = Path(sys.argv[1])
source_dir = resource_dir / "all_scaffolds" / "ppi_scaffolds"
ranges = {
    "scaffolds_50_65": (50, 65),
    "scaffolds_60_75": (60, 75),
}
manifest_rows = []

for output_name, (minimum, maximum) in ranges.items():
    output_dir = resource_dir / output_name
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = 0
    observed_lengths = []
    for ss_path in sorted(source_dir.glob("*_ss.pt")):
        scaffold_name = ss_path.name.removesuffix("_ss.pt")
        adj_path = source_dir / f"{scaffold_name}_adj.pt"
        if not adj_path.is_file():
            continue
        length = len(torch.load(ss_path, map_location="cpu"))
        if minimum <= length <= maximum:
            shutil.copy2(ss_path, output_dir / ss_path.name)
            shutil.copy2(adj_path, output_dir / adj_path.name)
            selected += 1
            observed_lengths.append(length)
            manifest_rows.append(
                {
                    "set": output_name,
                    "scaffold": scaffold_name,
                    "length": length,
                }
            )
    if selected == 0:
        raise RuntimeError(f"No scaffolds selected for {output_name}")
    print(
        output_name,
        f"pairs={selected}",
        f"observed_min={min(observed_lengths)}",
        f"observed_max={max(observed_lengths)}",
    )

manifest_path = resource_dir / "scaffold_manifest.tsv"
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["set", "scaffold", "length"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)
PY

for required in \
    "${resource_dir}/target_folds/USP15_DUSP_3T9L_A6-134_ss.pt" \
    "${resource_dir}/target_folds/USP15_DUSP_3T9L_A6-134_adj.pt" \
    "${resource_dir}/scaffold_manifest.tsv"
do
    if [[ ! -s "${required}" ]]; then
        echo "Missing prepared resource: ${required}"
        exit 1
    fi
done

echo "R2 Phase B resources prepared without PyRosetta"
