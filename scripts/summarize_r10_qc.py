#!/usr/bin/env python3
"""Summarize R10 sequence-composition, ProteinSol, and ESM-IF diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ENTROPY_MIN = 2.5
GRAVY_MAX = 1.0
PROTEINSOL_SCALED_MIN = 0.30


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--seq-composition", required=True, type=Path)
    parser.add_argument("--proteinsol", required=True, type=Path)
    parser.add_argument("--esmif", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path) -> dict[str, dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["id"]: row for row in csv.DictReader(handle)}


def finite(row: dict, key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {key} for {row.get('id')}")
    return value


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    sequence_rows = read_rows(args.seq_composition)
    proteinsol_rows = read_rows(args.proteinsol)
    esmif_rows = read_rows(args.esmif)

    rows = []
    for record in manifest["records"]:
        candidate_id = record["id"]
        seq = sequence_rows[candidate_id]
        sol = proteinsol_rows[candidate_id]
        esm = esmif_rows[candidate_id]
        avg_entropy = finite(seq, "avg_entropy")
        gravy = finite(seq, "gravy")
        scaled_sol = finite(sol, "scaled-sol")
        esmif_probability = finite(esm, "native_seq_avg_softmax")
        low_complexity_pass = avg_entropy >= ENTROPY_MIN
        hydrophobicity_pass = gravy <= GRAVY_MAX
        solubility_pass = scaled_sol >= PROTEINSOL_SCALED_MIN
        qc_pass = low_complexity_pass and hydrophobicity_pass and solubility_pass
        rows.append(
            {
                "id": candidate_id,
                "avg_entropy": avg_entropy,
                "gravy": gravy,
                "proteinsol_scaled": scaled_sol,
                "esmif_native_seq_avg_softmax": esmif_probability,
                "low_complexity_pass": low_complexity_pass,
                "hydrophobicity_pass": hydrophobicity_pass,
                "solubility_pass": solubility_pass,
                "qc_pass": qc_pass,
            }
        )
    passing_ids = [row["id"] for row in rows if row["qc_pass"]]
    summary = {
        "phase": "R10 ProteinQC",
        "gates": {
            "avg_entropy_min": ENTROPY_MIN,
            "gravy_max": GRAVY_MAX,
            "proteinsol_scaled_min": PROTEINSOL_SCALED_MIN,
        },
        "esm_if_role": "ranking diagnostic only",
        "input_count": len(rows),
        "passing_count": len(passing_ids),
        "passing_ids": passing_ids,
        "rows": rows,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    with args.csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"input_count": len(rows), "passing_count": len(passing_ids)}))


if __name__ == "__main__":
    main()
