#!/usr/bin/env python3
"""Summarize R3 LigandMPNN redesign smoke with unchanged AF2 gates."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True, type=Path)
    parser.add_argument("--sequence-pdb-dir", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    return parser.parse_args()


def binder_sequence(path: Path) -> str:
    residues: dict[tuple[int, str], str] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM") and line[21] == "A":
            key = (int(line[22:26]), line[26])
            residues.setdefault(key, AA3_TO_1[line[17:20].strip()])
    return "".join(residues[key] for key in sorted(residues))


def main() -> int:
    args = parse_args()
    sequence_paths = sorted(args.sequence_pdb_dir.glob("*.pdb"))
    sequences = {path.stem: binder_sequence(path) for path in sequence_paths}
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
                    ipae <= 10.0
                    and rmsd <= 2.0
                    and plddt >= 80.0
                    and "C" not in sequence
                ),
            }
        )
    if len(records) != 3 or len(sequence_paths) != 3 or len(rows) != 3:
        issues.append(
            f"Expected 3 sequences/AF2 records; found "
            f"{len(sequence_paths)}/{len(records)}/{len(rows)}"
        )
    if len(set(sequences.values())) != len(sequences):
        issues.append("LigandMPNN produced duplicate binder sequences")

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

    pass_count = sum(row["passed_all_af2_gates"] for row in rows)
    if issues:
        status = "technical_incomplete"
        next_action = "repair_and_resume"
    elif pass_count >= 3:
        status = "sequence_smoke_passed"
        next_action = "proceed_to_selectivity_and_focused_replication"
    elif pass_count >= 1:
        status = "partial_sequence_smoke_pass"
        next_action = "expand_ligandmpnn_sampling_around_same_scaffold"
    else:
        status = "sequence_smoke_failed"
        next_action = "optimize_scaffold_pose_or_interface_design"
    summary = {
        "phase": "R3 UbV-scaffold LigandMPNN sequence smoke",
        "status": status,
        "expected_sequences": 3,
        "observed_sequences": len(rows),
        "passing_sequences": pass_count,
        "parameters": {
            "temperature": 0.1,
            "sequences_per_backbone": 3,
            "omit_amino_acids": "C",
            "amino_acid_bias": None,
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
