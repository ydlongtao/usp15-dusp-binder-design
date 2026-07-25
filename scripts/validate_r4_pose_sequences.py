#!/usr/bin/env python3
"""Validate R4 ProteinMPNN counts, Cys exclusion, and fixed residues."""

from __future__ import annotations

import argparse
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
POSES = ("AK", "BL", "CJ", "DH")
FIXED_INTERFACE = (4, 6, 7, 8, 9, 44, 46, 48, 49, 50, 51, 72, 73, 74, 75)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--sequence-pdb-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def binder_sequence(path: Path) -> str:
    residues: dict[int, str] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[21] != "A":
            continue
        residue = int(line[22:26])
        residues.setdefault(residue, AA3_TO_1[line[17:20].strip()])
    if sorted(residues) != list(range(1, 77)):
        raise ValueError(f"{path.name}: expected binder residues A1-76")
    return "".join(residues[index] for index in range(1, 77))


def main() -> int:
    args = parse_args()
    issues: list[str] = []
    records: list[dict[str, object]] = []
    all_sequences: list[str] = []
    for pose in POSES:
        input_path = args.input_dir / f"pose_{pose}_fixed_interface.pdb"
        if not input_path.is_file():
            issues.append(f"Missing input pose {input_path.name}")
            continue
        source_sequence = binder_sequence(input_path)
        output_paths = sorted(
            args.sequence_pdb_dir.glob(f"pose_{pose}_fixed_interface*.pdb")
        )
        if len(output_paths) != 3:
            issues.append(
                f"Pose {pose}: expected 3 output PDBs, found {len(output_paths)}"
            )
        pose_sequences: list[str] = []
        for output_path in output_paths:
            sequence = binder_sequence(output_path)
            pose_sequences.append(sequence)
            all_sequences.append(sequence)
            fixed_preserved = all(
                sequence[index - 1] == source_sequence[index - 1]
                for index in FIXED_INTERFACE
            )
            if "C" in sequence:
                issues.append(f"{output_path.name}: binder contains Cys")
            if not fixed_preserved:
                issues.append(
                    f"{output_path.name}: crystallographic interface changed"
                )
            records.append(
                {
                    "pose": pose,
                    "id": output_path.stem,
                    "sequence": sequence,
                    "length": len(sequence),
                    "binder_cys_count": sequence.count("C"),
                    "fixed_interface_preserved": fixed_preserved,
                    "fixed_interface_sequence": "".join(
                        sequence[index - 1] for index in FIXED_INTERFACE
                    ),
                }
            )
        if len(set(pose_sequences)) != len(pose_sequences):
            issues.append(f"Pose {pose}: ProteinMPNN sequences are not unique")
    if len(records) != 12:
        issues.append(f"Expected 12 sequence records, found {len(records)}")
    if len(set(all_sequences)) != len(all_sequences):
        issues.append("The 12 ProteinMPNN sequences are not globally unique")

    report = {
        "phase": "R4 6DJ9 crystal-pose ensemble",
        "status": "valid" if not issues else "invalid",
        "expected_pose_count": 4,
        "expected_sequences_per_pose": 3,
        "expected_sequence_count": 12,
        "observed_sequence_count": len(records),
        "global_unique_sequence_count": len(set(all_sequences)),
        "fixed_interface_positions": list(FIXED_INTERFACE),
        "records": records,
        "issues": issues,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
