#!/usr/bin/env python3
"""Summarize bounded R5 AFDesign smoke sequences under unchanged OVO gates."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--afdesign-report-dir", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    design_reports = {
        int(path.stem.split("_")[1]): json.loads(path.read_text(encoding="utf-8"))
        for path in args.afdesign_report_dir.glob("seed_*_afdesign_report.json")
    }
    records = [
        json.loads(line)
        for line in args.jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows: list[dict[str, object]] = []
    issues: list[str] = []
    for record in records:
        design_id = record["id"]
        try:
            seed = int(design_id.split("_")[1])
        except (IndexError, ValueError):
            issues.append(f"Cannot extract seed from {design_id}")
            continue
        report = design_reports.get(seed)
        if report is None:
            issues.append(f"Missing AFDesign report for seed {seed}")
            continue
        sequence = report["sequence"]
        ipae = float(record["ipae"])
        rmsd = float(record["target_aligned_binder_rmsd"])
        plddt = float(record["binder_plddt"])
        rows.append(
            {
                "id": design_id,
                "seed": seed,
                "sequence": sequence,
                "length": len(sequence),
                "binder_cys_count": sequence.count("C"),
                "design_models": ",".join(report["design_models"]),
                "binder_structure_template_during_design": report.get(
                    "binder_structure_template_during_design", False
                ),
                "model_in_the_loop": report.get("model_in_the_loop", False),
                "ipae": ipae,
                "target_aligned_binder_rmsd": rmsd,
                "binder_plddt": plddt,
                "passed_all_af2_gates": (
                    len(sequence) == 76
                    and "C" not in sequence
                    and ipae <= 10.0
                    and rmsd <= 2.0
                    and plddt >= 80.0
                ),
            }
        )
    if len(records) != 4 or len(rows) != 4:
        issues.append(
            f"Expected four AFDesign/OVO records; found {len(records)}/{len(rows)}"
        )
    if len({row["sequence"] for row in rows}) != len(rows):
        issues.append("AFDesign smoke sequences are not unique")
    passing = [row for row in rows if row["passed_all_af2_gates"]]

    args.csv_output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else []
    with args.csv_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    args.fasta_output.write_text(
        "".join(f">{row['id']}\n{row['sequence']}\n" for row in rows),
        encoding="utf-8",
    )
    summary = {
        "phase": "R5 AFDesign fixed-RFD1 sequence optimization",
        "status": (
            "technical_incomplete"
            if issues
            else "positive_gate_passed"
            if len(passing) >= 3
            else "positive_gate_failed"
        ),
        "expected_sequences": 4,
        "observed_sequences": len(rows),
        "passing_sequences": len(passing),
        "thresholds": {
            "ipae_max": 10.0,
            "target_aligned_binder_rmsd_A_max": 2.0,
            "binder_plddt_min": 80.0,
        },
        "rows": rows,
        "technical_issues": issues,
        "next_action": (
            "repair_and_resume"
            if issues
            else "run_usp4_usp11_selectivity"
            if len(passing) >= 3
            else "stop_afdesign_branch_and_review_validation_protocol"
        ),
    }
    args.json_output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
