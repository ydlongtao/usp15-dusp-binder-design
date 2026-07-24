#!/usr/bin/env python3
"""Summarize standardized RFdiffusion backbones and enforce the preview gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


HOTSPOTS = {50, 52, 53, 55, 57, 61}


def read_ca(path: Path) -> list[tuple[str, int, tuple[float, float, float]]]:
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        atoms.append(
            (
                line[21],
                int(line[22:26]),
                (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            )
        )
    return atoms


def evaluate(path: Path) -> dict[str, object]:
    atoms = read_ca(path)
    binder = [atom for atom in atoms if atom[0] == "A"]
    target = [atom for atom in atoms if atom[0] == "B"]
    if not binder or not target:
        raise ValueError(f"Expected standardized binder A and target B in {path}")

    contacts = []
    for binder_atom in binder:
        for target_atom in target:
            distance = math.dist(binder_atom[2], target_atom[2])
            if distance <= 10.0:
                contacts.append((binder_atom[1], target_atom[1], distance))
    hotspot_contacts = [contact for contact in contacts if contact[1] in HOTSPOTS]
    touched_hotspots = sorted({contact[1] for contact in hotspot_contacts})
    return {
        "pdb": str(path),
        "binder_length": len(binder),
        "target_length": len(target),
        "N_contact_interface": len(contacts),
        "N_contact_hotspots": len(hotspot_contacts),
        "N_hotspots_on_interface": len(touched_hotspots),
        "hotspots_on_interface": ",".join(map(str, touched_hotspots)),
        "passes_backbone_gate": len(hotspot_contacts) >= 8 and len(touched_hotspots) >= 4,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_dir", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--require-all-pass", action="store_true")
    args = parser.parse_args()

    rows = [evaluate(path) for path in sorted(args.pdb_dir.glob("*.pdb"))]
    if not rows:
        raise ValueError(f"No PDB files found in {args.pdb_dir}")

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "designs": len(rows),
        "passed": sum(bool(row["passes_backbone_gate"]) for row in rows),
        "failed": sum(not bool(row["passes_backbone_gate"]) for row in rows),
        "rows": rows,
    }
    args.json.write_text(json.dumps(summary, indent=2) + "\n")
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps({key: summary[key] for key in ("designs", "passed", "failed")}))
    if args.require_all_pass and summary["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
