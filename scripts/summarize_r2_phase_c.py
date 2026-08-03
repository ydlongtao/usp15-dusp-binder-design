#!/usr/bin/env python3
"""Summarize R2 Phase C LigandMPNN and AF2 results."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
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


def find_af2_jsonl(run_dir: Path) -> Path | None:
    matches = sorted(
        (run_dir / "af2" / "output").glob(f"**/{AF2_TEST}.jsonl")
    )
    if len(matches) > 1:
        raise RuntimeError(f"Multiple AF2 JSONL files found in {run_dir}")
    return matches[0] if matches else None


def count_sequence_pdbs(run_dir: Path) -> int:
    path = (
        run_dir
        / "sequence"
        / "output"
        / "batch1"
        / "ligandmpnn"
        / "standardized_pdb"
    )
    return len(list(path.glob("*.pdb")))


def main() -> int:
    args = parse_args()
    with args.matrix.open(newline="", encoding="utf-8") as handle:
        matrix_rows = list(csv.DictReader(handle, delimiter="\t"))

    metric_rows: list[dict[str, Any]] = []
    technical_issues: list[str] = []
    expected_designs = 0
    for matrix_row in matrix_rows:
        expected = int(matrix_row["sequences_per_backbone"])
        expected_designs += expected
        run_dir = args.phase_dir / "runs" / matrix_row["run_id"]
        observed_sequences = count_sequence_pdbs(run_dir)
        if observed_sequences != expected:
            technical_issues.append(
                f"{matrix_row['run_id']}: expected {expected} sequence PDBs, "
                f"found {observed_sequences}"
            )

        jsonl_path = find_af2_jsonl(run_dir)
        if jsonl_path is None:
            technical_issues.append(
                f"{matrix_row['run_id']}: missing {AF2_TEST}.jsonl"
            )
            continue
        records = [
            json.loads(line)
            for line in jsonl_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if len(records) != expected:
            technical_issues.append(
                f"{matrix_row['run_id']}: expected {expected} AF2 records, "
                f"found {len(records)}"
            )
        for record in records:
            ipae = float(record["ipae"])
            binder_rmsd = float(record["target_aligned_binder_rmsd"])
            binder_plddt = float(record["binder_plddt"])
            metric_rows.append(
                {
                    "run_id": matrix_row["run_id"],
                    "condition_id": matrix_row["condition_id"],
                    "candidate_id": matrix_row["candidate_id"],
                    "backbone_file": matrix_row["backbone_file"],
                    "temperature": float(matrix_row["temperature"]),
                    "design_id": record["id"],
                    "af2_test": AF2_TEST,
                    "ipae": ipae,
                    "target_aligned_binder_rmsd": binder_rmsd,
                    "binder_plddt": binder_plddt,
                    "pass_ipae": ipae <= IPAE_MAX,
                    "pass_binder_rmsd": binder_rmsd <= BINDER_RMSD_MAX,
                    "pass_binder_plddt": binder_plddt >= BINDER_PLDDT_MIN,
                    "passed_all_af2_gates": (
                        ipae <= IPAE_MAX
                        and binder_rmsd <= BINDER_RMSD_MAX
                        and binder_plddt >= BINDER_PLDDT_MIN
                    ),
                }
            )

    fieldnames = [
        "run_id",
        "condition_id",
        "candidate_id",
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

    passed = [row for row in metric_rows if row["passed_all_af2_gates"]]
    per_condition = Counter(row["condition_id"] for row in passed)
    if technical_issues:
        status = "technical_incomplete"
        next_action = "repair_and_resume_phase_c"
    elif passed:
        status = "af2_gate_passed"
        next_action = "proceed_to_positive_quality_and_selectivity_evaluation"
    else:
        status = "af2_gate_failed"
        next_action = "report_r2_nonconvergence_without_threshold_relaxation"

    def best(field: str, reverse: bool = False) -> dict[str, Any] | None:
        if not metric_rows:
            return None
        row = sorted(
            metric_rows,
            key=lambda item: float(item[field]),
            reverse=reverse,
        )[0]
        return {
            "run_id": row["run_id"],
            "condition_id": row["condition_id"],
            "candidate_id": row["candidate_id"],
            "temperature": row["temperature"],
            "design_id": row["design_id"],
            field: row[field],
        }

    summary = {
        "phase": "USP15 R2 Phase C",
        "status": status,
        "matrix_runs": len(matrix_rows),
        "expected_af2_designs": expected_designs,
        "observed_af2_designs": len(metric_rows),
        "passed_all_af2_gates": len(passed),
        "passing_designs_per_condition": dict(sorted(per_condition.items())),
        "thresholds": {
            "ipae_max": IPAE_MAX,
            "target_aligned_binder_rmsd_A_max": BINDER_RMSD_MAX,
            "binder_plddt_min": BINDER_PLDDT_MIN,
        },
        "best_metrics_independently": {
            "ipae": best("ipae"),
            "target_aligned_binder_rmsd": best(
                "target_aligned_binder_rmsd"
            ),
            "binder_plddt": best("binder_plddt", reverse=True),
        },
        "passing_designs": passed,
        "technical_issues": technical_issues,
        "next_action": next_action,
        "thresholds_relaxed": False,
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
