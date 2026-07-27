#!/usr/bin/env python3
"""Gate formal R10 MD on a complete, correctly sized CUDA smoke."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re


def version_tuple(value: str):
    return tuple(int(part) for part in re.findall(r"\d+", value)[:3])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--analysis", required=True, type=Path)
    parser.add_argument("--preparation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    status = json.loads(args.status.read_text())
    analysis = json.loads(args.analysis.read_text())
    preparation = json.loads(args.preparation.read_text())
    expected_atoms = int(preparation["protein_atoms"])

    checks = {
        "run_completed": status.get("status") == "completed",
        "analysis_completed": analysis.get("status") == "completed",
        "platform_cuda": status.get("platform") == "CUDA",
        "openmm_at_least_8_5": (
            version_tuple(status.get("openmm_version", "0")) >= (8, 5)
        ),
        "production_unrestrained": (
            status.get("production_restraints") is False
        ),
        "protein_atom_subset_exact": (
            status.get("protein_trajectory_atoms") == expected_atoms
        ),
        "production_steps_exact": status.get("production_steps") == 25_000,
        "five_10ps_frames": analysis.get("frames") == 5,
        "observed_duration_0p05ns": (
            abs(float(analysis.get("duration_ns_observed", -1)) - 0.05)
            < 1e-9
        ),
    }
    report = {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "expected_protein_atoms": expected_atoms,
        "observed_protein_atoms": status.get("protein_trajectory_atoms"),
        "observed_frames": analysis.get("frames"),
        "observed_duration_ns": analysis.get("duration_ns_observed"),
        "wall_seconds": status.get("wall_seconds"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
