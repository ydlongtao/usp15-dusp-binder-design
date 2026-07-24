#!/usr/bin/env python3
"""Summarize OVO backbone metrics with no third-party dependencies."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("metrics_csv", type=Path)
    parser.add_argument("--filtered-dir", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    with args.metrics_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No metrics rows in {args.metrics_csv}")

    passed = [row for row in rows if as_bool(row.get("passed_filters", ""))]
    filtered_pdbs = list(args.filtered_dir.glob("*.pdb"))
    if len(passed) != len(filtered_pdbs):
        raise ValueError(
            f"Filtered output mismatch: CSV passed={len(passed)}, PDB links={len(filtered_pdbs)}"
        )

    summary = {
        "metrics_csv": str(args.metrics_csv),
        "generated_backbones": len(rows),
        "passed_backbone_gate": len(passed),
        "pass_rate": len(passed) / len(rows),
        "hard_filters": {
            "N_contact_hotspots": ">=8",
            "N_hotspots_on_interface": ">=4",
        },
        "N_contact_hotspots": {
            "minimum": min(float(row["N_contact_hotspots"]) for row in rows),
            "maximum": max(float(row["N_contact_hotspots"]) for row in rows),
        },
        "N_hotspots_on_interface": {
            "minimum": min(float(row["N_hotspots_on_interface"]) for row in rows),
            "maximum": max(float(row["N_hotspots_on_interface"]) for row in rows),
        },
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
