#!/usr/bin/env python3
"""Prepare and validate the USP15 DUSP target used by the OVO campaign."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
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

HOTSPOTS = (50, 52, 53, 55, 57, 61)
REQUIRED_BACKBONE = frozenset({"N", "CA", "C", "O"})


def parse_atom(line: str) -> dict[str, object]:
    return {
        "line": line.rstrip("\n"),
        "serial": int(line[6:11]),
        "atom": line[12:16].strip(),
        "altloc": line[16],
        "resname": line[17:20].strip(),
        "chain": line[21],
        "resseq": int(line[22:26]),
        "icode": line[26],
        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
    }


def rewrite_atom(atom: dict[str, object], serial: int) -> str:
    line = str(atom["line"]).ljust(80)
    # Select altloc A when it is the only modeled conformation, then normalize to blank.
    return f"{line[:6]}{serial:5d}{line[11:16]} {line[17:21]}A{line[22:80]}".rstrip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    atoms: list[dict[str, object]] = []
    altloc_seen: set[str] = set()
    with args.source.open() as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            atom = parse_atom(line)
            if atom["chain"] != "A" or not (6 <= int(atom["resseq"]) <= 134):
                continue
            if atom["resname"] not in AA3_TO_1:
                raise ValueError(f"Non-standard residue: {atom['resname']} {atom['resseq']}")
            if str(atom["icode"]).strip():
                raise ValueError(f"Insertion code at A{atom['resseq']}{atom['icode']}")
            altloc = str(atom["altloc"])
            if altloc not in {" ", "A"}:
                altloc_seen.add(altloc)
                continue
            atoms.append(atom)

    by_residue: dict[int, list[dict[str, object]]] = defaultdict(list)
    for atom in atoms:
        by_residue[int(atom["resseq"])].append(atom)

    expected = list(range(6, 135))
    observed = sorted(by_residue)
    if observed != expected:
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise ValueError(f"Residue coverage mismatch; missing={missing}, extra={extra}")

    missing_backbone: dict[str, list[str]] = {}
    sequence: list[str] = []
    residues: dict[str, str] = {}
    atom_lookup: dict[tuple[int, str], tuple[float, float, float]] = {}
    for resseq in observed:
        residue_atoms = by_residue[resseq]
        names = {str(atom["atom"]) for atom in residue_atoms}
        absent = sorted(REQUIRED_BACKBONE - names)
        if absent:
            missing_backbone[f"A{resseq}"] = absent
        resname = str(residue_atoms[0]["resname"])
        sequence.append(AA3_TO_1[resname])
        residues[f"A{resseq}"] = resname
        for atom in residue_atoms:
            atom_lookup[(resseq, str(atom["atom"]))] = atom["xyz"]  # type: ignore[assignment]

    if missing_backbone:
        raise ValueError(f"Missing backbone atoms: {missing_backbone}")

    peptide_cn_distances: dict[str, float] = {}
    chain_breaks: list[dict[str, object]] = []
    for left, right in zip(expected, expected[1:]):
        distance = math.dist(atom_lookup[(left, "C")], atom_lookup[(right, "N")])
        peptide_cn_distances[f"A{left}-A{right}"] = round(distance, 3)
        if distance > 2.0:
            chain_breaks.append({"between": [f"A{left}", f"A{right}"], "c_n_distance": round(distance, 3)})

    if chain_breaks:
        raise ValueError(f"Backbone chain breaks detected: {chain_breaks}")

    missing_hotspots = [f"A{x}" for x in HOTSPOTS if x not in by_residue]
    if missing_hotspots:
        raise ValueError(f"Missing hotspots: {missing_hotspots}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for serial, atom in enumerate(atoms, start=1):
            handle.write(rewrite_atom(atom, serial) + "\n")
        last = atoms[-1]
        handle.write(
            f"TER   {len(atoms) + 1:5d}      {last['resname']:>3} A{int(last['resseq']):4d}\n"
        )
        handle.write("END\n")

    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(args.source.resolve()),
        "output": str(args.output.resolve()),
        "source_structure": "3T9L",
        "chain": "A",
        "residue_range": [6, 134],
        "residue_count": len(observed),
        "atom_count": len(atoms),
        "sequence": "".join(sequence),
        "hotspots": [f"A{x}" for x in HOTSPOTS],
        "hotspot_residues": {f"A{x}": residues[f"A{x}"] for x in HOTSPOTS},
        "checks": {
            "only_standard_atom_records": True,
            "single_chain_A": True,
            "continuous_residue_numbering": True,
            "no_insertion_codes": True,
            "complete_backbone": True,
            "no_peptide_chain_breaks_over_2A": True,
            "all_hotspots_present": True,
            "discarded_nonblank_nonA_altlocs": sorted(altloc_seen),
        },
        "peptide_c_n_distance_A": {
            "minimum": min(peptide_cn_distances.values()),
            "maximum": max(peptide_cn_distances.values()),
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
