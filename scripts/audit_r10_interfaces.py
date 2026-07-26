#!/usr/bin/env python3
"""Apply non-PyRosetta interface area, clash, and hotspot gates to R10 passers."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley


DELTA_SASA_MIN = 600.0
CONTACT_HOTSPOTS_MIN = 8
HOTSPOTS_ON_INTERFACE_MIN = 4
SEVERE_CLASH_DISTANCE = 1.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-summary", required=True, type=Path)
    parser.add_argument("--panel-csv", required=True, type=Path)
    parser.add_argument("--panel-dir", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def interface_sasa_and_clashes(path: Path) -> tuple[float, float, float, int]:
    structure = PDBParser(QUIET=True).get_structure(path.stem, path)
    model = structure[0]
    if set(chain.id for chain in model) != {"A", "B"}:
        raise ValueError(f"{path}: expected only chains A and B")
    chain_a = model["A"]
    chain_b = model["B"]
    sr = ShrakeRupley(n_points=100)

    sr.compute(model, level="C")
    complex_a = float(chain_a.sasa)
    complex_b = float(chain_b.sasa)
    complex_total = complex_a + complex_b
    sr.compute(chain_a, level="C")
    mono_a = float(chain_a.sasa)
    sr.compute(chain_b, level="C")
    mono_b = float(chain_b.sasa)
    delta_sasa = mono_a + mono_b - complex_total

    atoms_a = np.asarray(
        [atom.coord for atom in chain_a.get_atoms() if atom.element != "H"],
        dtype=float,
    )
    atoms_b = np.asarray(
        [atom.coord for atom in chain_b.get_atoms() if atom.element != "H"],
        dtype=float,
    )
    severe_clashes = 0
    chunk_size = 256
    cutoff_squared = SEVERE_CLASH_DISTANCE**2
    for start in range(0, len(atoms_a), chunk_size):
        distances_squared = np.sum(
            (
                atoms_a[start : start + chunk_size, None, :]
                - atoms_b[None, :, :]
            )
            ** 2,
            axis=2,
        )
        severe_clashes += int(np.count_nonzero(distances_squared < cutoff_squared))
    return delta_sasa, mono_a - complex_a, mono_b - complex_b, severe_clashes


def main() -> None:
    args = parse_args()
    positive = json.loads(args.positive_summary.read_text(encoding="utf-8"))
    passing_ids = set(positive["passing_ids"])
    with args.panel_csv.open(encoding="utf-8", newline="") as handle:
        panel_rows = {row["id"]: row for row in csv.DictReader(handle)}
    if not passing_ids:
        raise ValueError("R10 positive screen has no passers")

    rows = []
    for candidate_id in sorted(passing_ids):
        panel = panel_rows[candidate_id]
        delta_sasa, binder_buried, target_buried, severe_clashes = (
            interface_sasa_and_clashes(args.panel_dir / f"{candidate_id}.pdb")
        )
        contact_hotspots = int(panel["N_contact_hotspots"])
        hotspots_on_interface = int(panel["N_hotspots_on_interface"])
        passed = (
            delta_sasa >= DELTA_SASA_MIN
            and severe_clashes == 0
            and contact_hotspots >= CONTACT_HOTSPOTS_MIN
            and hotspots_on_interface >= HOTSPOTS_ON_INTERFACE_MIN
        )
        rows.append(
            {
                "id": candidate_id,
                "interface_delta_sasa_A2": delta_sasa,
                "binder_buried_sasa_A2": binder_buried,
                "target_buried_sasa_A2": target_buried,
                "severe_clash_pair_count": severe_clashes,
                "N_contact_hotspots": contact_hotspots,
                "N_hotspots_on_interface": hotspots_on_interface,
                "interface_pass": passed,
            }
        )

    passing = [row["id"] for row in rows if row["interface_pass"]]
    summary = {
        "phase": "R10 non-PyRosetta interface audit",
        "gates": {
            "interface_delta_sasa_A2_min": DELTA_SASA_MIN,
            "severe_clash_pair_count_max": 0,
            "severe_clash_distance_A": SEVERE_CLASH_DISTANCE,
            "N_contact_hotspots_min": CONTACT_HOTSPOTS_MIN,
            "N_hotspots_on_interface_min": HOTSPOTS_ON_INTERFACE_MIN,
        },
        "input_count": len(rows),
        "passing_count": len(passing),
        "passing_ids": passing,
        "rows": rows,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    with args.csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"input_count": len(rows), "passing_count": len(passing)}))


if __name__ == "__main__":
    main()
