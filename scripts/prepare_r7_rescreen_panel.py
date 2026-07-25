#!/usr/bin/env python3
"""Build an exact-sequence-unique, candidate-eligible R7 re-screen panel."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path


AA3 = {
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
HOTSPOTS = {50, 52, 53, 55, 57, 61}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        action="append",
        required=True,
        help="Source specification LABEL=/path/to/pdb_dir; order defines duplicate priority",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    return parser.parse_args()


def parse_pdb(path: Path) -> dict:
    sequences = {"A": [], "B": []}
    ca = {"A": [], "B": []}
    seen = set()
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        if chain not in sequences:
            continue
        residue_key = (chain, line[22:27])
        if residue_key not in seen:
            seen.add(residue_key)
            sequences[chain].append(AA3.get(line[17:20].strip(), "X"))
        if line[12:16].strip() == "CA":
            ca[chain].append(
                (
                    int(line[22:26]),
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ),
                )
            )
    return {
        "binder_sequence": "".join(sequences["A"]),
        "target_sequence": "".join(sequences["B"]),
        "binder_ca": ca["A"],
        "target_ca": ca["B"],
    }


def backbone_metrics(parsed: dict) -> dict:
    contacts = []
    for binder_residue, binder_coord in parsed["binder_ca"]:
        for target_residue, target_coord in parsed["target_ca"]:
            distance = math.dist(binder_coord, target_coord)
            if distance <= 10.0:
                contacts.append((binder_residue, target_residue))
    hotspot_contacts = [
        contact for contact in contacts if contact[1] in HOTSPOTS
    ]
    touched = sorted({contact[1] for contact in hotspot_contacts})
    return {
        "N_contact_interface": len(contacts),
        "N_contact_hotspots": len(hotspot_contacts),
        "N_hotspots_on_interface": len(touched),
        "hotspots_on_interface": ",".join(map(str, touched)),
        "passes_backbone_gate": len(hotspot_contacts) >= 8 and len(touched) >= 4,
    }


def main() -> None:
    args = parse_args()
    sources = []
    for spec in args.source:
        if "=" not in spec:
            raise ValueError(f"Expected LABEL=/path, got {spec}")
        label, raw_path = spec.split("=", 1)
        sources.append((label, Path(raw_path)))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    seen_sequences = {}
    selected = []
    for source_label, source_dir in sources:
        if not source_dir.is_dir():
            raise ValueError(f"Missing source directory: {source_dir}")
        for path in sorted(source_dir.glob("*.pdb")):
            parsed = parse_pdb(path)
            binder = parsed["binder_sequence"]
            target = parsed["target_sequence"]
            metrics = backbone_metrics(parsed)
            reasons = []
            if not 45 <= len(binder) <= 80:
                reasons.append("binder_length_outside_45_80")
            if len(target) != 129:
                reasons.append("target_not_complete_3t9l_129aa")
            if "C" in binder:
                reasons.append("binder_contains_cys")
            if "X" in binder or "X" in target:
                reasons.append("nonstandard_or_unknown_residue")
            if not metrics["passes_backbone_gate"]:
                reasons.append("failed_backbone_hotspot_gate")
            duplicate_of = seen_sequences.get(binder)
            if duplicate_of:
                reasons.append(f"exact_sequence_duplicate_of:{duplicate_of}")

            output_id = f"{source_label}__{path.stem}"
            eligible = not reasons
            row = {
                "id": output_id,
                "source_label": source_label,
                "source_pdb": str(path),
                "binder_sequence": binder,
                "binder_length": len(binder),
                "target_length": len(target),
                **metrics,
                "selected": eligible,
                "elimination_reasons": ";".join(reasons),
            }
            rows.append(row)
            if eligible:
                destination = args.output_dir / f"{output_id}.pdb"
                shutil.copy2(path, destination)
                seen_sequences[binder] = output_id
                selected.append(row)

    if not selected:
        raise ValueError("No candidate-eligible R7 panel members")
    summary = {
        "phase": "R7 prior-design re-screen panel",
        "source_count": len(sources),
        "examined": len(rows),
        "selected_exact_sequence_unique": len(selected),
        "rejected": len(rows) - len(selected),
        "invariants": {
            "binder_length": "45-80",
            "complete_target_length": 129,
            "binder_cys": 0,
            "N_contact_hotspots_min": 8,
            "N_hotspots_on_interface_min": 4,
        },
        "rows": rows,
    }

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    with args.csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    with args.fasta.open("w", encoding="utf-8") as handle:
        for row in selected:
            handle.write(f">{row['id']}\n{row['binder_sequence']}\n")
    print(
        json.dumps(
            {
                "examined": len(rows),
                "selected_exact_sequence_unique": len(selected),
                "rejected": len(rows) - len(selected),
            }
        )
    )


if __name__ == "__main__":
    main()
