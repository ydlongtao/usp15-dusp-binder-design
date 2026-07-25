#!/usr/bin/env python3
"""Apply R2 compactness/topology filters and select five backbones."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import Counter
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("metrics_csv", type=Path)
    parser.add_argument("--pdb-dir", required=True, type=Path)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--filtered-csv", required=True, type=Path)
    parser.add_argument("--selected-dir", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--select-count", type=int, default=5)
    return parser.parse_args()


def numeric(row: dict[str, str], field: str) -> float:
    return float(row[field])


def structured_segments(dssp: str, minimum_length: int = 4) -> int:
    segments = 0
    current_type = ""
    current_length = 0
    for state in dssp:
        state_type = state if state in {"H", "E"} else "-"
        if state_type == current_type:
            current_length += 1
            continue
        if current_type in {"H", "E"} and current_length >= minimum_length:
            segments += 1
        current_type = state_type
        current_length = 1
    if current_type in {"H", "E"} and current_length >= minimum_length:
        segments += 1
    return segments


def longest_helix_fraction(dssp: str) -> float:
    if not dssp:
        return 1.0
    longest = 0
    current = 0
    for state in dssp:
        if state == "H":
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest / len(dssp)


def evaluate(row: dict[str, str]) -> dict[str, Any]:
    dssp = row["pydssp_str"]
    binder_length = len(dssp)
    rg_max = 15.5 if binder_length <= 55 else 18.0
    contact_density = numeric(row, "N_contact_interface") / binder_length
    segment_count = structured_segments(dssp)
    helix_fraction = longest_helix_fraction(dssp)
    loop_percent = numeric(row, "pydssp_loop_percent")
    gates = {
        "pass_contact_hotspots": numeric(row, "N_contact_hotspots") >= 8,
        "pass_hotspots_on_interface": (
            numeric(row, "N_hotspots_on_interface") >= 4
        ),
        "pass_radius_of_gyration": numeric(row, "radius_of_gyration") <= rg_max,
        "pass_structured_segments": segment_count >= 2,
        "pass_longest_helix_fraction": helix_fraction <= 0.65,
        "pass_loop_percent": 5.0 <= loop_percent <= 40.0,
        "pass_contact_density": contact_density >= 1.0,
    }
    result: dict[str, Any] = dict(row)
    result.update(
        {
            "binder_length": binder_length,
            "rg_max": rg_max,
            "contact_density": contact_density,
            "structured_segments": segment_count,
            "longest_helix_fraction": helix_fraction,
            **gates,
            "passed_r2_filters": all(gates.values()),
        }
    )
    return result


def selection_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -float(row["N_hotspots_on_interface"]),
        -float(row["N_contact_hotspots"]),
        -float(row["contact_density"]),
        float(row["radius_of_gyration"]),
        float(row["longest_helix_fraction"]),
        row["id"],
    )


def main() -> int:
    args = parse_args()
    with args.metrics_csv.open(newline="", encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle))
    if not source_rows:
        raise ValueError(f"No metrics rows in {args.metrics_csv}")

    rows = [evaluate(row) for row in source_rows]
    passed = [row for row in rows if row["passed_r2_filters"]]
    selected = sorted(passed, key=selection_key)[: args.select_count]

    fieldnames = list(rows[0])
    args.filtered_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.filtered_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    args.selected_dir.mkdir(parents=True, exist_ok=True)
    for stale_pdb in args.selected_dir.glob("*.pdb"):
        stale_pdb.unlink()
    for row in selected:
        source_pdb = args.pdb_dir / f"{row['id']}.pdb"
        if not source_pdb.is_file():
            raise FileNotFoundError(source_pdb)
        shutil.copy2(source_pdb, args.selected_dir / source_pdb.name)

    filter_names = [
        "pass_contact_hotspots",
        "pass_hotspots_on_interface",
        "pass_radius_of_gyration",
        "pass_structured_segments",
        "pass_longest_helix_fraction",
        "pass_loop_percent",
        "pass_contact_density",
    ]
    summary = {
        "condition": args.condition,
        "generated_backbones": len(rows),
        "passed_r2_filters": len(passed),
        "selected_for_phase_c": len(selected),
        "selected_ids": [row["id"] for row in selected],
        "filter_pass_counts": {
            name: sum(bool(row[name]) for row in rows) for name in filter_names
        },
        "length_distribution": dict(
            sorted(Counter(int(row["binder_length"]) for row in rows).items())
        ),
        "thresholds": {
            "N_contact_hotspots_min": 8,
            "N_hotspots_on_interface_min": 4,
            "radius_of_gyration_max_A": {
                "45-55": 15.5,
                "56-75": 18.0,
            },
            "structured_segments_min": 2,
            "longest_continuous_helix_fraction_max": 0.65,
            "loop_percent_range": [5.0, 40.0],
            "contact_density_min": 1.0,
        },
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
