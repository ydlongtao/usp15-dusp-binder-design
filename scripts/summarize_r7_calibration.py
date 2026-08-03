#!/usr/bin/env python3
"""Summarize R7 two-control AF2 model/dropout-seed calibration."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


IPAE_MAX = 10.0
RMSD_MAX = 2.0
PLDDT_MIN = 80.0
CONTROLS = {"exact_native_6dj9", "complete_3t9l_6dj9_ubv"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def passes(record: dict) -> bool:
    return (
        float(record["ipae"]) <= IPAE_MAX
        and float(record["target_aligned_binder_rmsd"]) <= RMSD_MAX
        and float(record["binder_plddt"]) >= PLDDT_MIN
    )


def main() -> None:
    args = parse_args()
    rows = []
    grouped: dict[str, dict[int, dict[str, bool]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    model_architecture = {}

    for path in sorted((args.phase_dir / "output").glob("*.jsonl")):
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            record = json.loads(line)
            passed = passes(record)
            model_name = record["model_name"]
            control = record["input_id"]
            if control not in CONTROLS:
                raise ValueError(f"Unexpected control {control}")
            grouped[model_name][int(record["seed"])][control] = passed
            model_architecture[model_name] = record["architecture"]
            rows.append(
                {
                    "input_id": control,
                    "architecture": record["architecture"],
                    "model_name": model_name,
                    "seed": int(record["seed"]),
                    "dropout": bool(record["dropout"]),
                    "num_recycles": int(record["num_recycles"]),
                    "template_mode": record["template_mode"],
                    "ipae": float(record["ipae"]),
                    "target_aligned_binder_rmsd_A": float(
                        record["target_aligned_binder_rmsd"]
                    ),
                    "binder_plddt": float(record["binder_plddt"]),
                    "passed_all_af2_gates": passed,
                }
            )

    expected_records = 7 * 3 * 2
    if len(rows) != expected_records:
        raise ValueError(f"Expected {expected_records} records, found {len(rows)}")

    model_results = {}
    for model_name, seeds in sorted(grouped.items()):
        seed_pass = {
            str(seed): set(control_results) == CONTROLS
            and all(control_results.values())
            for seed, control_results in sorted(seeds.items())
        }
        passing_seed_count = sum(seed_pass.values())
        model_results[model_name] = {
            "architecture": model_architecture[model_name],
            "seed_pass_both_controls": seed_pass,
            "passing_seed_count": passing_seed_count,
            "calibrated": passing_seed_count >= 2,
        }

    calibrated_by_architecture = {
        architecture: [
            model_name
            for model_name, result in model_results.items()
            if result["architecture"] == architecture and result["calibrated"]
        ]
        for architecture in ("ptm", "multimer")
    }
    calibration_passed = all(calibrated_by_architecture.values())
    selected_models = {
        architecture: (
            sorted(
                calibrated_by_architecture[architecture],
                key=lambda name: (
                    -model_results[name]["passing_seed_count"],
                    name,
                ),
            )[0]
            if calibrated_by_architecture[architecture]
            else None
        )
        for architecture in ("ptm", "multimer")
    }
    summary = {
        "phase": "R7 AF2 model/dropout-seed calibration",
        "candidate_eligible": False,
        "thresholds_unchanged": {
            "ipae_max": IPAE_MAX,
            "target_aligned_binder_rmsd_A_max": RMSD_MAX,
            "binder_plddt_min": PLDDT_MIN,
        },
        "controls": sorted(CONTROLS),
        "expected_records": expected_records,
        "observed_records": len(rows),
        "model_results": model_results,
        "calibrated_by_architecture": calibrated_by_architecture,
        "selected_models": selected_models,
        "calibration_passed": calibration_passed,
        "decision": (
            "proceed_to_r7_candidate_rescreen"
            if calibration_passed
            else "stop_no_reproducible_ptm_and_multimer_calibration"
        ),
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
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
