#!/usr/bin/env python3
"""Summarize R2 Phase A LigandMPNN and AF2 smoke results."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


AF2_TEST = "af2_model_1_multimer_tt_3rec"
IPAE_MAX = 10.0
BINDER_RMSD_MAX = 2.0
BINDER_PLDDT_MIN = 80.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    return parser.parse_args()


def read_matrix(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "run_id",
        "candidate_id",
        "source_pool",
        "backbone_file",
        "temperature",
        "sequences_per_backbone",
        "omit_amino_acids",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Invalid Phase A matrix: {path}")
    return rows


def load_af2_rows(jsonl_path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with jsonl_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                records.append(json.loads(line))
    return records


def find_af2_jsonl(run_dir: Path) -> Path | None:
    matches = sorted(
        (run_dir / "af2" / "output").glob(f"**/{AF2_TEST}.jsonl")
    )
    if len(matches) > 1:
        raise RuntimeError(f"Multiple AF2 JSONL files found in {run_dir}")
    return matches[0] if matches else None


def count_sequence_pdbs(run_dir: Path) -> int:
    standardized = (
        run_dir
        / "sequence"
        / "output"
        / "batch1"
        / "ligandmpnn"
        / "standardized_pdb"
    )
    return len(list(standardized.glob("*.pdb")))


def metric_row(
    matrix_row: dict[str, str], record: dict[str, Any]
) -> dict[str, Any]:
    ipae = float(record["ipae"])
    binder_rmsd = float(record["target_aligned_binder_rmsd"])
    binder_plddt = float(record["binder_plddt"])
    pass_ipae = ipae <= IPAE_MAX
    pass_rmsd = binder_rmsd <= BINDER_RMSD_MAX
    pass_plddt = binder_plddt >= BINDER_PLDDT_MIN
    return {
        "run_id": matrix_row["run_id"],
        "candidate_id": matrix_row["candidate_id"],
        "source_pool": matrix_row["source_pool"],
        "backbone_file": matrix_row["backbone_file"],
        "temperature": float(matrix_row["temperature"]),
        "design_id": record["id"],
        "af2_test": AF2_TEST,
        "ipae": ipae,
        "target_aligned_binder_rmsd": binder_rmsd,
        "binder_plddt": binder_plddt,
        "pass_ipae": pass_ipae,
        "pass_binder_rmsd": pass_rmsd,
        "pass_binder_plddt": pass_plddt,
        "passed_all_af2_gates": pass_ipae and pass_rmsd and pass_plddt,
    }


def best_record(
    rows: list[dict[str, Any]], field: str, reverse: bool = False
) -> dict[str, Any] | None:
    if not rows:
        return None
    selected = sorted(rows, key=lambda row: float(row[field]), reverse=reverse)[0]
    return {
        "run_id": selected["run_id"],
        "candidate_id": selected["candidate_id"],
        "temperature": selected["temperature"],
        "design_id": selected["design_id"],
        field: selected[field],
    }


def main() -> int:
    args = parse_args()
    matrix_rows = read_matrix(args.matrix)
    metric_rows: list[dict[str, Any]] = []
    technical_issues: list[str] = []
    expected_designs = 0

    for matrix_row in matrix_rows:
        expected = int(matrix_row["sequences_per_backbone"])
        expected_designs += expected
        run_dir = args.phase_dir / "runs" / matrix_row["run_id"]
        sequence_count = count_sequence_pdbs(run_dir)
        if sequence_count != expected:
            technical_issues.append(
                f"{matrix_row['run_id']}: expected {expected} sequence PDBs, "
                f"found {sequence_count}"
            )

        jsonl_path = find_af2_jsonl(run_dir)
        if jsonl_path is None:
            technical_issues.append(
                f"{matrix_row['run_id']}: missing {AF2_TEST}.jsonl"
            )
            continue

        records = load_af2_rows(jsonl_path)
        if len(records) != expected:
            technical_issues.append(
                f"{matrix_row['run_id']}: expected {expected} AF2 records, "
                f"found {len(records)}"
            )
        metric_rows.extend(metric_row(matrix_row, record) for record in records)

    fieldnames = [
        "run_id",
        "candidate_id",
        "source_pool",
        "backbone_file",
        "temperature",
        "design_id",
        "af2_test",
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
        writer.writerows(metric_rows)

    passed_rows = [row for row in metric_rows if row["passed_all_af2_gates"]]
    if technical_issues:
        status = "technical_incomplete"
        next_action = "repair_and_resume_phase_a"
    elif passed_rows:
        status = "af2_gate_passed"
        next_action = "prepare_targeted_scaling_without_starting_automatically"
    else:
        status = "af2_gate_failed"
        next_action = "proceed_to_r2_phase_b_only_after_review"

    summary = {
        "phase": "USP15 R2 Phase A",
        "status": status,
        "matrix_runs": len(matrix_rows),
        "expected_af2_designs": expected_designs,
        "observed_af2_designs": len(metric_rows),
        "passed_all_af2_gates": len(passed_rows),
        "thresholds": {
            "ipae_max": IPAE_MAX,
            "target_aligned_binder_rmsd_A_max": BINDER_RMSD_MAX,
            "binder_plddt_min": BINDER_PLDDT_MIN,
        },
        "best_metrics_independently": {
            "ipae": best_record(metric_rows, "ipae"),
            "target_aligned_binder_rmsd": best_record(
                metric_rows, "target_aligned_binder_rmsd"
            ),
            "binder_plddt": best_record(
                metric_rows, "binder_plddt", reverse=True
            ),
        },
        "passing_designs": passed_rows,
        "technical_issues": technical_issues,
        "next_action": next_action,
        "automatic_scaling_started": False,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if technical_issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
