#!/usr/bin/env python3
"""Select up to three hard-filtered partial-diffusion backbones per condition."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--per-condition", type=int, default=3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    matrix = list(csv.DictReader(args.matrix.open(encoding="utf-8"), delimiter="\t"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for condition in matrix:
        condition_id = condition["condition_id"]
        metrics_csv = (
            args.phase_dir
            / "conditions"
            / condition_id
            / "metrics"
            / "output"
            / condition_id
            / "backbone_metrics.csv"
        )
        filtered_dir = metrics_csv.parent / "backbones_filtered"
        rows = list(csv.DictReader(metrics_csv.open(encoding="utf-8")))
        passing = [row for row in rows if truthy(row.get("passed_filters", ""))]
        passing.sort(
            key=lambda row: (
                float(row["N_hotspots_on_interface"]),
                float(row["N_contact_hotspots"]),
                float(row["N_contact_interface"]),
            ),
            reverse=True,
        )
        selected = passing[: args.per_condition]
        for rank, row in enumerate(selected, start=1):
            source_name = f"{row['id']}.pdb"
            source = filtered_dir / source_name
            if not source.exists():
                matches = list(filtered_dir.glob(f"*{Path(source_name).stem}*.pdb"))
                if len(matches) != 1:
                    raise ValueError(
                        f"Cannot resolve filtered PDB for {condition_id}: {source_name}"
                    )
                source = matches[0]
            output_name = f"{condition_id}_rank{rank}_{source.name}"
            destination = args.output_dir / output_name
            shutil.copy2(source, destination)
            records.append(
                {
                    "condition_id": condition_id,
                    "rank": rank,
                    "source": str(source),
                    "selected_pdb": output_name,
                    "N_contact_hotspots": float(row["N_contact_hotspots"]),
                    "N_hotspots_on_interface": float(
                        row["N_hotspots_on_interface"]
                    ),
                    "N_contact_interface": float(row["N_contact_interface"]),
                }
            )
    report = {
        "selected_backbones": len(records),
        "maximum_per_condition": args.per_condition,
        "hard_filters": {
            "N_contact_hotspots": ">=8",
            "N_hotspots_on_interface": ">=4",
        },
        "records": records,
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
