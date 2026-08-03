#!/usr/bin/env python3
"""Prepare sequence-only Boltz YAML controls from standardized PDB inputs."""

from __future__ import annotations

import argparse
import json
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-native", required=True, type=Path)
    parser.add_argument("--complete-target", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sequences(path: Path) -> dict[str, str]:
    chains = {"A": [], "B": []}
    seen = set()
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[21] not in chains:
            continue
        key = (line[21], line[22:27])
        if key in seen:
            continue
        seen.add(key)
        chains[line[21]].append(AA3.get(line[17:20].strip(), "X"))
    result = {chain: "".join(items) for chain, items in chains.items()}
    if not result["A"] or not result["B"] or "X" in result["A"] + result["B"]:
        raise ValueError(f"Invalid standardized protein chains in {path}")
    return result


def yaml_text(chain_sequences: dict[str, str]) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {chain_sequences['A']}\n"
        "  - protein:\n"
        "      id: B\n"
        f"      sequence: {chain_sequences['B']}\n"
    )


def main() -> None:
    args = parse_args()
    controls = {
        "exact_native_6dj9": (args.exact_native, sequences(args.exact_native)),
        "complete_3t9l_6dj9_ubv": (
            args.complete_target,
            sequences(args.complete_target),
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for control_id, (source, chain_sequences) in controls.items():
        (args.output_dir / f"{control_id}.yaml").write_text(
            yaml_text(chain_sequences),
            encoding="utf-8",
        )
        records.append(
            {
                "id": control_id,
                "source": str(source),
                "binder_chain": "A",
                "target_chain": "B",
                "binder_length": len(chain_sequences["A"]),
                "target_length": len(chain_sequences["B"]),
                "candidate_eligible": False,
            }
        )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(
            {
                "phase": "R8 Boltz-2 positive-control preparation",
                "template": None,
                "forced_constraints": False,
                "controls": records,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"controls": len(records)}))


if __name__ == "__main__":
    main()
