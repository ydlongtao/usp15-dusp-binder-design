#!/usr/bin/env python3
"""Compute three-state PyDSSP diagnostics for final R10 binder chains."""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--pydssp-module", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_pydssp(path: Path):
    spec = importlib.util.spec_from_file_location("ovo_pydssp_numpy", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    args = parse_args()
    pydssp = load_pydssp(args.pydssp_module)
    rows = []
    for path in sorted(args.input_dir.glob("*.pdb")):
        coordinates = pydssp.read_pdbtext_with_checking(
            path.read_text(encoding="utf-8", errors="ignore"),
            chain_id="A",
        )
        assignment = pydssp.assign(coordinates)
        counts = np.asarray(assignment, dtype=int).sum(axis=0)
        length = int(counts.sum())
        rows.append(
            {
                "id": path.stem,
                "length": length,
                "loop_fraction": float(counts[0] / length),
                "helix_fraction": float(counts[1] / length),
                "strand_fraction": float(counts[2] / length),
            }
        )
    if not rows:
        raise ValueError("No R10 PyDSSP inputs")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print({"input_count": len(rows)})


if __name__ == "__main__":
    main()
