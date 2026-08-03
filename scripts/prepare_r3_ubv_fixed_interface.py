#!/usr/bin/env python3
"""Mark the crystallographic UbV interface fixed and redesign its fold core."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


FIXED_INTERFACE = {4, 6, 7, 8, 9, 44, 46, 48, 49, 50, 51, 72, 73, 74, 75}


def compact_ranges(residues: list[int]) -> str:
    ranges = []
    start = previous = residues[0]
    for residue in residues[1:]:
        if residue == previous + 1:
            previous = residue
        else:
            ranges.append((start, previous))
            start = previous = residue
    ranges.append((start, previous))
    return "/".join(
        f"A{start}" if start == end else f"A{start}-{end}"
        for start, end in ranges
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    atom_lines = [
        line
        for line in args.input.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]
    binder_residues = sorted(
        {int(line[22:26]) for line in atom_lines if line[21] == "A"}
    )
    target_residues = sorted(
        {int(line[22:26]) for line in atom_lines if line[21] == "B"}
    )
    if binder_residues != list(range(1, 77)):
        raise ValueError("Expected binder A1-76")
    if target_residues != list(range(6, 135)):
        raise ValueError("Expected target B6-134")
    redesigned = sorted(set(binder_residues) - FIXED_INTERFACE)
    inpaint_seq = compact_ranges(redesigned)
    header = [
        'REMARK   1 Input contig: "A1-76/0 B6-134/0"',
        'REMARK   1 Standardized contig: "A1-76/0 B6-134/0"',
        f'REMARK   1 Inpaint seq: "{inpaint_seq}"',
        'REMARK   1 Chains: "A B"',
        'REMARK   1 Input hotspots: "B50,B52,B53,B55,B57,B61"',
        'REMARK   1 Standardized hotspots: "B50,B52,B53,B55,B57,B61"',
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "\n".join(header + atom_lines + ["END"]) + "\n",
        encoding="utf-8",
    )
    report = {
        "fixed_crystallographic_interface_positions": sorted(FIXED_INTERFACE),
        "redesigned_core_positions": redesigned,
        "fixed_position_count": len(FIXED_INTERFACE),
        "redesigned_position_count": len(redesigned),
        "inpaint_seq": inpaint_seq,
        "temperature": 0.1,
        "sequences_per_backbone": 3,
        "omit_amino_acids": "C",
        "sequence_model": "protein_mpnn",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
