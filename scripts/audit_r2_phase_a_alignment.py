#!/usr/bin/env python3
"""Audit R2 Phase A chain mapping and target-template alignment."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


AF2_TEST = "af2_model_1_multimer_tt_3rec"
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
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--native-target", required=True, type=Path)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--target-rmsd-max", type=float, default=2.0)
    parser.add_argument("--metric-tolerance", type=float, default=0.02)
    return parser.parse_args()


def parse_pdb(path: Path) -> dict[str, dict[str, Any]]:
    chains: dict[str, dict[str, Any]] = {}
    seen_residues: set[tuple[str, str]] = set()
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip()
            residue_key = line[22:27]
            residue_name = line[17:20].strip()
            chain_data = chains.setdefault(
                chain, {"sequence": [], "ca": [], "residue_keys": []}
            )
            if (chain, residue_key) not in seen_residues:
                seen_residues.add((chain, residue_key))
                chain_data["sequence"].append(AA3.get(residue_name, "X"))
                chain_data["residue_keys"].append(residue_key)
            if line[12:16].strip() == "CA":
                chain_data["ca"].append(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ]
                )
    for chain_data in chains.values():
        chain_data["sequence"] = "".join(chain_data["sequence"])
        chain_data["ca"] = np.asarray(chain_data["ca"], dtype=float)
    return chains


def superpose(
    fixed: np.ndarray, moving: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    if fixed.shape != moving.shape or fixed.ndim != 2 or fixed.shape[1] != 3:
        raise ValueError(
            f"Coordinate shape mismatch: fixed={fixed.shape}, moving={moving.shape}"
        )
    fixed_center = fixed.mean(axis=0)
    moving_center = moving.mean(axis=0)
    fixed_zero = fixed - fixed_center
    moving_zero = moving - moving_center
    covariance = moving_zero.T @ fixed_zero
    left, _, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    translation = fixed_center - moving_center @ rotation
    fitted = moving @ rotation + translation
    rmsd = float(
        np.sqrt(np.mean(np.sum((fitted - fixed) ** 2, axis=1)))
    )
    return rotation, translation, rmsd


def load_metrics(path: Path) -> dict[str, dict[str, Any]]:
    metrics: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                record = json.loads(line)
                metrics[record["id"]] = record
    return metrics


def main() -> int:
    args = parse_args()
    native = parse_pdb(args.native_target)
    native_target_sequence = native["A"]["sequence"]
    rows: list[dict[str, Any]] = []
    issues: list[str] = []

    for run_dir in sorted((args.phase_dir / "runs").glob("r2a_*")):
        sequence_dir = (
            run_dir
            / "sequence"
            / "output"
            / "batch1"
            / "ligandmpnn"
            / "standardized_pdb"
        )
        af2_root = run_dir / "af2" / "output" / "contig1_batch1"
        prediction_dir = af2_root / AF2_TEST
        metrics_path = af2_root / f"{AF2_TEST}.jsonl"
        if not metrics_path.is_file():
            issues.append(f"{run_dir.name}: missing {metrics_path.name}")
            continue
        metrics = load_metrics(metrics_path)

        for prediction_path in sorted(prediction_dir.glob("*.pdb")):
            suffix = f"_{AF2_TEST}.pdb"
            design_id = prediction_path.name.removesuffix(suffix)
            design_path = sequence_dir / f"{design_id}.pdb"
            if not design_path.is_file() or design_id not in metrics:
                issues.append(f"{run_dir.name}/{design_id}: missing input or metrics")
                continue

            design = parse_pdb(design_path)
            prediction = parse_pdb(prediction_path)
            expected_chains = {"A", "B"}
            chain_mapping_ok = (
                set(design) == expected_chains
                and set(prediction) == expected_chains
            )
            target_sequence_ok = (
                design.get("B", {}).get("sequence") == native_target_sequence
                and prediction.get("B", {}).get("sequence")
                == native_target_sequence
            )
            binder_sequence_ok = (
                design.get("A", {}).get("sequence")
                == prediction.get("A", {}).get("sequence")
            )

            if not chain_mapping_ok:
                issues.append(f"{run_dir.name}/{design_id}: chain mapping mismatch")
                continue

            rotation, translation, target_ca_rmsd = superpose(
                design["B"]["ca"], prediction["B"]["ca"]
            )
            fitted_binder = prediction["A"]["ca"] @ rotation + translation
            binder_ca_rmsd = float(
                np.sqrt(
                    np.mean(
                        np.sum(
                            (fitted_binder - design["A"]["ca"]) ** 2,
                            axis=1,
                        )
                    )
                )
            )
            reported_binder_rmsd = float(
                metrics[design_id]["target_aligned_binder_rmsd"]
            )
            metric_delta = abs(binder_ca_rmsd - reported_binder_rmsd)
            row = {
                "run_id": run_dir.name,
                "design_id": design_id,
                "design_binder_length": len(design["A"]["sequence"]),
                "design_target_length": len(design["B"]["sequence"]),
                "prediction_binder_length": len(prediction["A"]["sequence"]),
                "prediction_target_length": len(prediction["B"]["sequence"]),
                "chain_mapping_ok": chain_mapping_ok,
                "target_sequence_ok": target_sequence_ok,
                "binder_sequence_ok": binder_sequence_ok,
                "target_ca_rmsd": target_ca_rmsd,
                "independent_binder_ca_rmsd": binder_ca_rmsd,
                "reported_binder_rmsd": reported_binder_rmsd,
                "binder_rmsd_delta": metric_delta,
                "audit_passed": (
                    chain_mapping_ok
                    and target_sequence_ok
                    and binder_sequence_ok
                    and target_ca_rmsd <= args.target_rmsd_max
                    and metric_delta <= args.metric_tolerance
                ),
            }
            rows.append(row)

    for row in rows:
        if not row["audit_passed"]:
            issues.append(
                f"{row['run_id']}/{row['design_id']}: alignment audit failed"
            )

    fieldnames = [
        "run_id",
        "design_id",
        "design_binder_length",
        "design_target_length",
        "prediction_binder_length",
        "prediction_target_length",
        "chain_mapping_ok",
        "target_sequence_ok",
        "binder_sequence_ok",
        "target_ca_rmsd",
        "independent_binder_ca_rmsd",
        "reported_binder_rmsd",
        "binder_rmsd_delta",
        "audit_passed",
    ]
    args.csv_output.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "phase": "USP15 R2 Phase A alignment audit",
        "status": "passed" if rows and not issues else "failed",
        "audited_designs": len(rows),
        "passed_designs": sum(bool(row["audit_passed"]) for row in rows),
        "target_ca_rmsd_range": (
            {
                "minimum": min(row["target_ca_rmsd"] for row in rows),
                "maximum": max(row["target_ca_rmsd"] for row in rows),
            }
            if rows
            else None
        ),
        "maximum_binder_rmsd_metric_delta": (
            max(row["binder_rmsd_delta"] for row in rows) if rows else None
        ),
        "issues": issues,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
