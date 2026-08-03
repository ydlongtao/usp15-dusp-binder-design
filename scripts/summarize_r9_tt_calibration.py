#!/usr/bin/env python3
"""Audit R9 model-2 pTM target-template positive-control calibration."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


CONTROLS = {"exact_native_6dj9", "complete_3t9l_6dj9_ubv"}
GATES = {
    "ipae_max": 10.0,
    "target_aligned_binder_rmsd_max": 2.0,
    "binder_plddt_min": 80.0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def passes(record: dict) -> bool:
    values = (
        record.get("ipae"),
        record.get("target_aligned_binder_rmsd"),
        record.get("binder_plddt"),
    )
    if any(value is None or not math.isfinite(float(value)) for value in values):
        return False
    return (
        float(record["ipae"]) <= GATES["ipae_max"]
        and float(record["target_aligned_binder_rmsd"])
        <= GATES["target_aligned_binder_rmsd_max"]
        and float(record["binder_plddt"]) >= GATES["binder_plddt_min"]
    )


def main() -> None:
    args = parse_args()
    records = [
        json.loads(line)
        for line in args.jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) != 6:
        raise ValueError(f"Expected six records, found {len(records)}")

    grouped: dict[int, dict[str, bool]] = {}
    rows = []
    for record in records:
        if record.get("architecture") != "ptm":
            raise ValueError("R9 calibration requires the pTM architecture")
        if record.get("model_name") != "model_2_ptm":
            raise ValueError(f"Unexpected model {record.get('model_name')}")
        if record.get("template_mode") != "tt":
            raise ValueError(f"Unexpected template mode {record.get('template_mode')}")
        if record.get("num_recycles") != 3 or not record.get("dropout"):
            raise ValueError("R9 recycle/dropout invariant violated")
        control = record["input_id"]
        seed = int(record["seed"])
        if control not in CONTROLS or seed not in (0, 1, 2):
            raise ValueError(f"Unexpected control/seed: {control}/{seed}")
        passed = passes(record)
        grouped.setdefault(seed, {})[control] = passed
        rows.append(
            {
                "input_id": control,
                "seed": seed,
                "ipae": float(record["ipae"]),
                "target_aligned_binder_rmsd_A": float(
                    record["target_aligned_binder_rmsd"]
                ),
                "binder_plddt": float(record["binder_plddt"]),
                "passed_all_gates": passed,
            }
        )

    seed_pass = {
        str(seed): set(grouped.get(seed, {})) == CONTROLS
        and all(grouped[seed].values())
        for seed in (0, 1, 2)
    }
    passing_seeds = [int(seed) for seed, passed in seed_pass.items() if passed]
    calibrated = len(passing_seeds) >= 2
    summary = {
        "phase": "R9 model-2 pTM target-template calibration",
        "candidate_eligible": False,
        "model": "model_2_ptm",
        "template_mode": "tt",
        "dropout": True,
        "seeds": [0, 1, 2],
        "recycles": 3,
        "gates": GATES,
        "seed_pass_both_controls": seed_pass,
        "passing_seeds": passing_seeds,
        "calibrated": calibrated,
        "decision": (
            "proceed_to_r9_candidate_screen"
            if calibrated
            else "stop_target_template_mode_not_calibrated"
        ),
        "records": sorted(rows, key=lambda row: (row["input_id"], row["seed"])),
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
        writer.writerows(summary["records"])
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
