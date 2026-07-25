#!/usr/bin/env python3
"""Independently audit Boltz-2 USP15 control predictions against fixed gates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


GATES = {
    "ipae_max": 10.0,
    "binder_rmsd_max": 2.0,
    "binder_plddt_min": 80.0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--exact-native", required=True, type=Path)
    parser.add_argument("--complete-target", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def ca_by_chain(path: Path) -> dict[str, np.ndarray]:
    coordinates: dict[str, list[list[float]]] = {"A": [], "B": []}
    seen: set[tuple[str, str]] = set()
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        chain = line[21]
        if chain not in coordinates:
            continue
        residue_id = line[22:27]
        key = (chain, residue_id)
        if key in seen:
            continue
        seen.add(key)
        coordinates[chain].append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    result = {chain: np.asarray(items) for chain, items in coordinates.items()}
    if any(array.ndim != 2 or array.shape[1:] != (3,) for array in result.values()):
        raise ValueError(f"Missing chain-A/chain-B C-alpha coordinates in {path}")
    return result


def kabsch_transform(
    mobile: np.ndarray, reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    if mobile.shape != reference.shape:
        raise ValueError(
            f"Alignment shape mismatch: mobile {mobile.shape}, reference {reference.shape}"
        )
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u_matrix @ vt_matrix))
    rotation = u_matrix @ correction @ vt_matrix
    translation = reference_center - mobile_center @ rotation
    return rotation, translation


def rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    if mobile.shape != reference.shape:
        raise ValueError(
            f"RMSD shape mismatch: mobile {mobile.shape}, reference {reference.shape}"
        )
    return float(np.sqrt(np.mean(np.sum((mobile - reference) ** 2, axis=1))))


def one_path(root: Path, pattern: str) -> Path:
    paths = sorted(root.rglob(pattern))
    if len(paths) != 1:
        raise ValueError(f"Expected one {pattern} below {root}, found {len(paths)}")
    return paths[0]


def audit_prediction(
    prediction_root: Path,
    reference_path: Path,
) -> dict[str, object]:
    prediction_path = one_path(prediction_root, "*_model_0.pdb")
    pae_path = one_path(prediction_root, "pae_*_model_0.npz")
    plddt_path = one_path(prediction_root, "plddt_*_model_0.npz")
    confidence_path = one_path(prediction_root, "confidence_*_model_0.json")

    reference = ca_by_chain(reference_path)
    predicted = ca_by_chain(prediction_path)
    if predicted["A"].shape != reference["A"].shape:
        raise ValueError(f"Binder C-alpha count mismatch in {prediction_path}")
    if predicted["B"].shape != reference["B"].shape:
        raise ValueError(f"Target C-alpha count mismatch in {prediction_path}")

    rotation, translation = kabsch_transform(predicted["B"], reference["B"])
    aligned_binder = predicted["A"] @ rotation + translation
    binder_rmsd = rmsd(aligned_binder, reference["A"])

    pae = np.load(pae_path)["pae"]
    plddt = np.load(plddt_path)["plddt"]
    binder_length = reference["A"].shape[0]
    target_length = reference["B"].shape[0]
    total_length = binder_length + target_length
    if pae.shape != (total_length, total_length):
        raise ValueError(
            f"PAE shape {pae.shape} does not match {total_length} tokens in {pae_path}"
        )
    if plddt.shape != (total_length,):
        raise ValueError(
            f"pLDDT shape {plddt.shape} does not match {total_length} tokens "
            f"in {plddt_path}"
        )
    cross_ab = pae[:binder_length, binder_length:]
    cross_ba = pae[binder_length:, :binder_length]
    ipae = float(np.concatenate((cross_ab.ravel(), cross_ba.ravel())).mean())
    binder_plddt = float(plddt[:binder_length].mean())
    if float(np.nanmax(plddt)) <= 1.5:
        binder_plddt *= 100.0

    passed = (
        ipae <= GATES["ipae_max"]
        and binder_rmsd <= GATES["binder_rmsd_max"]
        and binder_plddt >= GATES["binder_plddt_min"]
    )
    return {
        "prediction_pdb": str(prediction_path),
        "pae_npz": str(pae_path),
        "plddt_npz": str(plddt_path),
        "confidence_json": str(confidence_path),
        "binder_length": binder_length,
        "target_length": target_length,
        "ipae": ipae,
        "target_aligned_binder_rmsd": binder_rmsd,
        "binder_plddt": binder_plddt,
        "pass": passed,
    }


def main() -> None:
    args = parse_args()
    references = {
        "exact_native_6dj9": args.exact_native,
        "complete_3t9l_6dj9_ubv": args.complete_target,
    }
    records = []
    for seed in args.seeds:
        predictions_dir = args.phase_dir / f"seed_{seed}" / "predictions"
        seed_pass = True
        control_records = []
        for control_id, reference in references.items():
            control_root = predictions_dir / control_id
            metrics = audit_prediction(control_root, reference)
            metrics["control_id"] = control_id
            metrics["reference_pdb"] = str(reference)
            control_records.append(metrics)
            seed_pass = seed_pass and bool(metrics["pass"])
        records.append(
            {
                "seed": seed,
                "controls": control_records,
                "both_controls_pass": seed_pass,
            }
        )

    passing_seeds = [record["seed"] for record in records if record["both_controls_pass"]]
    summary = {
        "phase": "R8 Boltz-2 positive-control calibration",
        "model": "boltz2",
        "templates": None,
        "forced_constraints": False,
        "recycles": 3,
        "sampling_steps": 200,
        "diffusion_samples": 1,
        "gates": GATES,
        "records": records,
        "passing_seeds": passing_seeds,
        "calibrated": len(passing_seeds) >= 2 and len(args.seeds) >= 3,
        "calibration_rule": "at least two of seeds 0,1,2 pass all gates on both controls",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
