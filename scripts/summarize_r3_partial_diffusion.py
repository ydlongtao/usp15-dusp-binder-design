#!/usr/bin/env python3
"""Summarize R3 partial-diffusion sequence/AF2 results."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from summarize_r3_sequence_smoke import binder_sequence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--sequence-pdb-dir", required=True, type=Path)
    parser.add_argument("--selection-report", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    parser.add_argument("--sequence-model", default="ligand_mpnn")
    parser.add_argument("--phase-label", default="R3 low-noise partial diffusion")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selection = json.loads(args.selection_report.read_text(encoding="utf-8"))
    expected = 3 * int(selection["selected_backbones"])
    sequences = {
        path.stem: binder_sequence(path)
        for path in sorted(args.sequence_pdb_dir.glob("*.pdb"))
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
            issues.append(f"Missing sequence PDB for {design_id}")
            continue
        ipae = float(record["ipae"])
        rmsd = float(record["target_aligned_binder_rmsd"])
        plddt = float(record["binder_plddt"])
        rows.append(
            {
                "id": design_id,
                "sequence": sequence,
                "length": len(sequence),
                "binder_cys_count": sequence.count("C"),
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
    if len(sequences) != expected or len(records) != expected or len(rows) != expected:
        issues.append(
            f"Expected {expected} sequences/AF2 rows; found "
            f"{len(sequences)}/{len(records)}/{len(rows)}"
        )
    if any(len(sequence) != 76 or "C" in sequence for sequence in sequences.values()):
        issues.append("Sequence length/Cys invariant failed")
    passing = [row for row in rows if row["passed_all_af2_gates"]]

    args.csv_output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else [
        "id",
        "sequence",
        "length",
        "binder_cys_count",
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
        next_action = "replicate_best_partial_T"
    else:
        status = "positive_gate_failed"
        next_action = "stop_partial_diffusion_branch"
    summary = {
        "phase": args.phase_label,
        "status": status,
        "selected_backbones": selection["selected_backbones"],
        "expected_sequences": expected,
        "observed_sequences": len(rows),
        "passing_sequences": len(passing),
        "parameters": {
            "sequence_model": args.sequence_model,
            "ligandmpnn_temperature": 0.1,
            "sequences_per_backbone": 3,
            "omit_amino_acids": "C",
            "af2_test": "af2_model_1_multimer_tt_3rec",
        },
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
