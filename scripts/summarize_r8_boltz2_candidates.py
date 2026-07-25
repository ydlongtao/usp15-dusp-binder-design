#!/usr/bin/env python3
"""Apply fixed gates to sequence-only Boltz-2 candidate predictions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from summarize_r8_boltz2 import GATES, audit_prediction


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--input-report", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_report = json.loads(args.input_report.read_text(encoding="utf-8"))
    candidates = input_report["records"]
    rows = []
    for candidate in candidates:
        seed_metrics = []
        for seed in (0, 1, 2):
            prediction_root = (
                args.phase_dir
                / f"seed_{seed}"
                / "predictions"
                / candidate["id"]
            )
            metrics = audit_prediction(
                prediction_root,
                Path(candidate["source_pdb"]),
            )
            metrics["seed"] = seed
            seed_metrics.append(metrics)
        passing_seeds = [
            metrics["seed"] for metrics in seed_metrics if metrics["pass"]
        ]
        row = {
            "id": candidate["id"],
            "passing_seed_count": len(passing_seeds),
            "passing_seeds": ",".join(map(str, passing_seeds)),
            "boltz_positive_pass": len(passing_seeds) >= 2,
        }
        for metrics in seed_metrics:
            seed = metrics["seed"]
            row[f"seed{seed}_ipae"] = metrics["ipae"]
            row[f"seed{seed}_binder_rmsd"] = metrics[
                "target_aligned_binder_rmsd"
            ]
            row[f"seed{seed}_binder_plddt"] = metrics["binder_plddt"]
        rows.append(row)

    passing_ids = [row["id"] for row in rows if row["boltz_positive_pass"]]
    summary = {
        "phase": "R8 Boltz-2 candidate ensemble",
        "model": "boltz2",
        "templates": None,
        "forced_constraints": False,
        "seeds": [0, 1, 2],
        "recycles": 3,
        "sampling_steps": 200,
        "diffusion_samples": 1,
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
