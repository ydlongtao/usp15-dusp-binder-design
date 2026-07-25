#!/usr/bin/env python3
"""Apply the fixed AF2 gates to the R8 pTM-model-2 candidate ensemble."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


GATES = {
    "ipae_max": 10.0,
    "target_aligned_binder_rmsd_max": 2.0,
    "binder_plddt_min": 80.0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--expected-inputs", required=True, type=int)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def seed_pass(record: dict) -> bool:
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
    expected_records = args.expected_inputs * 3
    if len(records) != expected_records:
        raise ValueError(f"Expected {expected_records} records, found {len(records)}")

    grouped: dict[str, list[dict]] = {}
    for record in records:
        if record.get("model_name") != "model_2_ptm":
            raise ValueError(f"Unexpected model: {record.get('model_name')}")
        if record.get("template_mode") != "ct":
            raise ValueError(f"Unexpected template mode: {record.get('template_mode')}")
        if record.get("num_recycles") != 3 or not record.get("dropout"):
            raise ValueError("R8 AF2 candidate record violates recycle/dropout invariants")
        grouped.setdefault(record["input_id"], []).append(record)
    if len(grouped) != args.expected_inputs:
        raise ValueError(
            f"Expected {args.expected_inputs} unique inputs, found {len(grouped)}"
        )

    rows = []
    for input_id, input_records in sorted(grouped.items()):
        input_records = sorted(input_records, key=lambda record: record["seed"])
        seeds = sorted(record["seed"] for record in input_records)
        if seeds != [0, 1, 2]:
            raise ValueError(f"{input_id}: expected seeds 0,1,2, found {seeds}")
        passing_seeds = sorted(
            record["seed"] for record in input_records if seed_pass(record)
        )
        rows.append(
            {
                "id": input_id,
                "passing_seed_count": len(passing_seeds),
                "passing_seeds": ",".join(map(str, passing_seeds)),
                "af2_positive_pass": len(passing_seeds) >= 2,
                "seed0_ipae": input_records[0]["ipae"],
                "seed0_binder_rmsd": input_records[0][
                    "target_aligned_binder_rmsd"
                ],
                "seed0_binder_plddt": input_records[0]["binder_plddt"],
                "seed1_ipae": input_records[1]["ipae"],
                "seed1_binder_rmsd": input_records[1][
                    "target_aligned_binder_rmsd"
                ],
                "seed1_binder_plddt": input_records[1]["binder_plddt"],
                "seed2_ipae": input_records[2]["ipae"],
                "seed2_binder_rmsd": input_records[2][
                    "target_aligned_binder_rmsd"
                ],
                "seed2_binder_plddt": input_records[2]["binder_plddt"],
            }
        )

    passing_ids = [row["id"] for row in rows if row["af2_positive_pass"]]
    summary = {
        "phase": "R8 AF2 pTM-model-2 candidate ensemble",
        "model": "model_2_ptm",
        "template_mode": "ct",
        "dropout": True,
        "seeds": [0, 1, 2],
        "recycles": 3,
        "gates": GATES,
        "candidate_rule": "at least two of three seeds pass all fixed gates",
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
