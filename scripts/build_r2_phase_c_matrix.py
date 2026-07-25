#!/usr/bin/env python3
"""Build the R2 Phase C sequence/AF2 matrix from Phase B selections."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-b-dir", required=True, type=Path)
    parser.add_argument("--phase-b-matrix", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.phase_b_matrix.open(newline="", encoding="utf-8") as handle:
        condition_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(condition_rows) != 6:
        raise ValueError("Phase B matrix must contain exactly six conditions")

    rows: list[dict[str, str | int | float]] = []
    condition_counts: dict[str, int] = {}
    for condition_row in condition_rows:
        condition = condition_row["condition_id"]
        condition_dir = args.phase_b_dir / "conditions" / condition
        summary_path = (
            condition_dir / "reports" / "phase_b_backbone_summary.json"
        )
        if not summary_path.is_file():
            raise FileNotFoundError(summary_path)
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if int(summary.get("generated_backbones", 0)) != 50:
            raise ValueError(
                f"{condition}: expected 50 generated backbones, got "
                f"{summary.get('generated_backbones')}"
            )

        selected_ids = list(summary.get("selected_ids", []))
        if len(selected_ids) > 5:
            raise ValueError(f"{condition}: more than five selected backbones")
        condition_counts[condition] = len(selected_ids)
        for rank, selected_id in enumerate(selected_ids, start=1):
            backbone_file = f"{selected_id}.pdb"
            source_path = (
                condition_dir / "selected_backbones" / backbone_file
            )
            if not source_path.is_file():
                raise FileNotFoundError(source_path)
            candidate_id = f"{condition.lower()}_{selected_id}"
            for temperature_label, temperature in (
                ("t005", 0.05),
                ("t010", 0.10),
            ):
                rows.append(
                    {
                        "run_id": (
                            f"r2c_{condition.lower()}_{rank:02d}_"
                            f"{temperature_label}"
                        ),
                        "condition_id": condition,
                        "candidate_id": candidate_id,
                        "backbone_file": backbone_file,
                        "temperature": temperature,
                        "sequences_per_backbone": 3,
                        "omit_amino_acids": "C",
                    }
                )

    if not rows:
        raise ValueError("No Phase B backbones were selected for Phase C")

    fieldnames = [
        "run_id",
        "condition_id",
        "candidate_id",
        "backbone_file",
        "temperature",
        "sequences_per_backbone",
        "omit_amino_acids",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "status": "completed",
        "conditions": condition_counts,
        "selected_backbones": sum(condition_counts.values()),
        "matrix_runs": len(rows),
        "expected_af2_designs": sum(
            int(row["sequences_per_backbone"]) for row in rows
        ),
        "temperatures": [0.05, 0.10],
        "sequences_per_backbone_per_temperature": 3,
        "omit_amino_acids": "C",
        "af2_test": "af2_model_1_multimer_tt_3rec",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
