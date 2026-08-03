#!/usr/bin/env python3
"""Summarize the exact-native R6 AF2 template-mode calibration."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TESTS = (
    "af2_model_1_ptm_tbt_3rec",
    "af2_model_1_multimer_tbt_3rec",
    "af2_model_1_ptm_ct_3rec",
    "af2_model_1_multimer_ct_3rec",
)
IPAE_MAX = 10.0
RMSD_MAX = 2.0
PLDDT_MIN = 80.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-dir", required=True, type=Path)
    parser.add_argument("--control-label", default="exact_6dj9_native")
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def read_record(path: Path) -> dict:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    if len(lines) != 1:
        raise ValueError(f"Expected one record in {path}, found {len(lines)}")
    return json.loads(lines[0])


def passes(record: dict) -> bool:
    return (
        float(record["ipae"]) <= IPAE_MAX
        and float(record["target_aligned_binder_rmsd"]) <= RMSD_MAX
        and float(record["binder_plddt"]) >= PLDDT_MIN
    )


def main() -> None:
    args = parse_args()
    rows = []
    records = {}
    for test_name in TESTS:
        path = (
            args.phase_dir
            / "runs"
            / test_name
            / "output"
            / "contig1_batch1"
            / f"{test_name}.jsonl"
        )
        record = read_record(path)
        passed = passes(record)
        row = {
            "test": test_name,
            "ipae": float(record["ipae"]),
            "target_aligned_binder_rmsd_A": float(
                record["target_aligned_binder_rmsd"]
            ),
            "binder_plddt": float(record["binder_plddt"]),
            "passed_all_af2_gates": passed,
        }
        rows.append(row)
        records[test_name] = {**record, "passed_all_af2_gates": passed}

    mode_pass = {
        mode: all(
            records[f"af2_model_1_{model}_{mode}_3rec"]["passed_all_af2_gates"]
            for model in ("ptm", "multimer")
        )
        for mode in ("tbt", "ct")
    }
    selected_mode = "tbt" if mode_pass["tbt"] else "ct" if mode_pass["ct"] else None
    summary = {
        "phase": "R6 AF2 template calibration",
        "control_label": args.control_label,
        "candidate_eligible": False,
        "thresholds_unchanged": {
            "ipae_max": IPAE_MAX,
            "target_aligned_binder_rmsd_A_max": RMSD_MAX,
            "binder_plddt_min": PLDDT_MIN,
        },
        "mode_pass": mode_pass,
        "selected_mode": selected_mode,
        "geometry_conditioned": selected_mode == "ct",
        "records": records,
        "decision": (
            f"use_{selected_mode}_for_r6_positive_rescreen"
            if selected_mode
            else "stop_no_template_mode_recovers_positive_control"
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
