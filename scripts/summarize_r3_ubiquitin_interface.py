#!/usr/bin/env python3
"""Summarize stable-ubiquitin scaffold/interface redesign AF2 results."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from summarize_r3_sequence_smoke import binder_sequence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--af2-input-dir", required=True, type=Path)
    parser.add_argument("--preparation-report", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    preparation = json.loads(args.preparation_report.read_text(encoding="utf-8"))
    wildtype = preparation["sequence"]
    designed_positions = set(preparation["designed_interface_residues"])
    sequences = {
        path.stem: binder_sequence(path)
        for path in sorted(args.af2_input_dir.glob("*.pdb"))
    }
    records = [
        json.loads(line)
        for line in args.jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows = []
    issues = []
    for record in records:
        design_id = record["id"]
        sequence = sequences.get(design_id)
        if sequence is None:
            issues.append(f"Missing AF2 input PDB for {design_id}")
            continue
        candidate_eligible = design_id != "ubiquitin_wt_control"
        changed_positions = [
            index
            for index, (original, designed) in enumerate(
                zip(wildtype, sequence), start=1
            )
            if original != designed
        ]
        if len(sequence) != 76:
            issues.append(f"{design_id}: expected 76-aa binder")
        if "C" in sequence:
            issues.append(f"{design_id}: binder contains Cys")
        if candidate_eligible and not set(changed_positions).issubset(
            designed_positions
        ):
            issues.append(f"{design_id}: non-interface scaffold mutation detected")
        ipae = float(record["ipae"])
        rmsd = float(record["target_aligned_binder_rmsd"])
        plddt = float(record["binder_plddt"])
        rows.append(
            {
                "id": design_id,
                "candidate_eligible": candidate_eligible,
                "sequence": sequence,
                "changed_positions_from_1ubq": ",".join(
                    map(str, changed_positions)
                ),
                "ipae": ipae,
                "target_aligned_binder_rmsd": rmsd,
                "binder_plddt": plddt,
                "passed_all_af2_gates": (
                    ipae <= 10.0
                    and rmsd <= 2.0
                    and plddt >= 80.0
                    and "C" not in sequence
                ),
            }
        )
    eligible = [row for row in rows if row["candidate_eligible"]]
    passing = [row for row in eligible if row["passed_all_af2_gates"]]
    if len(records) != 4 or len(rows) != 4 or len(eligible) != 3:
        issues.append(
            f"Expected four AF2 records and three eligible designs; found "
            f"{len(records)}/{len(rows)}/{len(eligible)}"
        )
    if len({row["sequence"] for row in eligible}) != len(eligible):
        issues.append("Interface design generated duplicate candidate sequences")

    args.csv_output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else [
        "id",
        "candidate_eligible",
        "sequence",
        "changed_positions_from_1ubq",
        "ipae",
        "target_aligned_binder_rmsd",
        "binder_plddt",
        "passed_all_af2_gates",
    ]
    with args.csv_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    args.fasta_output.write_text(
        "".join(f">{row['id']}\n{row['sequence']}\n" for row in rows),
        encoding="utf-8",
    )

    if issues:
        status = "technical_incomplete"
        next_action = "repair_and_resume"
    elif len(passing) >= 3:
        status = "positive_gate_passed"
        next_action = "run_usp4_usp11_selectivity"
    elif passing:
        status = "partial_positive_gate_pass"
        next_action = "replicate_interface_sampling"
    else:
        status = "positive_gate_failed"
        next_action = "evaluate_alternative_stable_scaffold_pose"
    summary = {
        "phase": "R3 stable ubiquitin scaffold/interface design",
        "status": status,
        "eligible_designs": len(eligible),
        "passing_eligible_designs": len(passing),
        "thresholds": {
            "ipae_max": 10.0,
            "target_aligned_binder_rmsd_A_max": 2.0,
            "binder_plddt_min": 80.0,
        },
        "rows": rows,
        "technical_issues": issues,
        "next_action": next_action,
    }
    args.json_output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
