#!/usr/bin/env python3
"""Summarize the R3 crystallographic UbV positive-control smoke."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


IPAE_MAX = 10.0
BINDER_RMSD_MAX = 2.0
BINDER_PLDDT_MIN = 80.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--preparation-report", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    preparation = json.loads(args.preparation_report.read_text(encoding="utf-8"))
    prepared = {record["id"]: record for record in preparation["variants"]}
    records = [
        json.loads(line)
        for line in args.jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows: list[dict[str, Any]] = []
    issues: list[str] = []
    for record in records:
        design_id = record["id"]
        base_id = next(
            (variant for variant in prepared if design_id.startswith(variant)),
            None,
        )
        if base_id is None:
            issues.append(f"Unexpected AF2 design id: {design_id}")
            continue
        ipae = float(record["ipae"])
        rmsd = float(record["target_aligned_binder_rmsd"])
        plddt = float(record["binder_plddt"])
        rows.append(
            {
                "id": base_id,
                "af2_design_id": design_id,
                "candidate_eligible": prepared[base_id]["candidate_eligible"],
                "binder_cys_count": prepared[base_id]["binder_cys_count"],
                "ipae": ipae,
                "target_aligned_binder_rmsd": rmsd,
                "binder_plddt": plddt,
                "pass_ipae": ipae <= IPAE_MAX,
                "pass_binder_rmsd": rmsd <= BINDER_RMSD_MAX,
                "pass_binder_plddt": plddt >= BINDER_PLDDT_MIN,
                "passed_all_af2_gates": (
                    ipae <= IPAE_MAX
                    and rmsd <= BINDER_RMSD_MAX
                    and plddt >= BINDER_PLDDT_MIN
                ),
            }
        )

    if len(rows) != 4:
        issues.append(f"Expected four AF2 records; found {len(rows)}")
    fieldnames = list(rows[0]) if rows else [
        "id",
        "af2_design_id",
        "candidate_eligible",
        "binder_cys_count",
        "ipae",
        "target_aligned_binder_rmsd",
        "binder_plddt",
        "pass_ipae",
        "pass_binder_rmsd",
        "pass_binder_plddt",
        "passed_all_af2_gates",
    ]
    args.csv_output.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    eligible_passes = [
        row
        for row in rows
        if row["candidate_eligible"] and row["passed_all_af2_gates"]
    ]
    control = next(
        (row for row in rows if row["id"] == "ubv15d_wt_control"),
        None,
    )
    if issues:
        status = "technical_incomplete"
        next_action = "repair_and_resume_r3_smoke"
    elif control is None or not control["passed_all_af2_gates"]:
        status = "positive_control_failed"
        next_action = "diagnose_evaluation_before_library_generation"
    elif len(eligible_passes) >= 3:
        status = "r3_smoke_passed"
        next_action = "proceed_to_focused_diversification_and_selectivity"
    else:
        status = "insufficient_cys_free_passes"
        next_action = "expand_conservative_cys_replacement_panel"

    summary = {
        "phase": "USP15 R3 UbV positive-control smoke",
        "status": status,
        "expected_designs": 4,
        "observed_designs": len(rows),
        "eligible_cys_free_passes": len(eligible_passes),
        "thresholds": {
            "ipae_max": IPAE_MAX,
            "target_aligned_binder_rmsd_A_max": BINDER_RMSD_MAX,
            "binder_plddt_min": BINDER_PLDDT_MIN,
        },
        "rows": rows,
        "technical_issues": issues,
        "next_action": next_action,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
