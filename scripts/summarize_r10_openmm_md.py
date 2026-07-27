#!/usr/bin/env python3
"""Aggregate completed USP15 R10 OpenMM replicas without relaxing gates."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import statistics


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--md-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_mmpbsa_mean(path):
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    patterns = [
        r"DELTA TOTAL\s+(-?\d+(?:\.\d+)?)",
        r"DELTA G binding\s*=\s*(-?\d+(?:\.\d+)?)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return float(match.group(1))
    return None


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    replicas = []

    for rank in range(1, 11):
        for seed in range(3):
            analysis_dir = (
                args.md_dir / "analysis" / f"rank{rank:02d}" / f"seed{seed}"
            )
            run_dir = (
                args.md_dir / "runs" / f"rank{rank:02d}" / f"seed{seed}"
            )
            summary_path = analysis_dir / "summary.json"
            status_path = run_dir / "status.json"
            if not summary_path.exists() or not status_path.exists():
                continue
            summary = json.loads(summary_path.read_text())
            status = json.loads(status_path.read_text())
            per_frame = read_csv(analysis_dir / "per_frame.csv")
            rmsd_le_3_fraction = (
                sum(float(row["binder_rmsd_a"]) <= 3.0 for row in per_frame)
                / len(per_frame)
            )
            mmpbsa = parse_mmpbsa_mean(
                analysis_dir / "mmpbsa" / "FINAL_RESULTS_MMPBSA.dat"
            )
            replicate_pass = (
                rmsd_le_3_fraction >= 0.70
                and summary["hotspots_with_occupancy_ge_0p5"] >= 4
                and summary["buried_sasa_a2"]["median"] >= 600.0
            )
            replicas.append(
                {
                    "rank": rank,
                    "seed": seed,
                    "production_ns": status.get(
                        "production_ns_requested"
                    ),
                    "frames": summary["frames"],
                    "binder_rmsd_mean_a": summary["binder_rmsd_a"]["mean"],
                    "binder_rmsd_q95_a": summary["binder_rmsd_a"]["q95"],
                    "binder_rmsd_le_3_fraction": rmsd_le_3_fraction,
                    "native_contact_fraction_median": summary[
                        "native_contact_fraction"
                    ]["median"],
                    "hotspots_with_occupancy_ge_0p5": summary[
                        "hotspots_with_occupancy_ge_0p5"
                    ],
                    "buried_sasa_median_a2": summary["buried_sasa_a2"][
                        "median"
                    ],
                    "persistent_interface_hbonds": summary[
                        "persistent_interface_hbonds_freq_ge_0p10"
                    ],
                    "mmpbsa_delta_total_kcal_mol": mmpbsa,
                    "replicate_stability_gate_pass": replicate_pass,
                    "trajectory": str(
                        run_dir / "production_protein.xtc"
                    ),
                }
            )

    replica_csv = args.output_dir / "replica_metrics.csv"
    if replicas:
        with replica_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(replicas[0]))
            writer.writeheader()
            writer.writerows(replicas)

    candidates = []
    for rank in range(1, 11):
        rows = [row for row in replicas if row["rank"] == rank]
        if not rows:
            continue
        passes = sum(row["replicate_stability_gate_pass"] for row in rows)
        mmpbsa_values = [
            row["mmpbsa_delta_total_kcal_mol"]
            for row in rows
            if row["mmpbsa_delta_total_kcal_mol"] is not None
        ]
        candidates.append(
            {
                "rank": rank,
                "completed_replicas": len(rows),
                "stability_gate_pass_replicas": passes,
                "three_replica_consensus": (
                    "robust" if len(rows) == 3 and passes == 3
                    else "majority" if len(rows) == 3 and passes >= 2
                    else "failed" if len(rows) == 3
                    else "incomplete"
                ),
                "binder_rmsd_mean_across_replicas_a": statistics.mean(
                    row["binder_rmsd_mean_a"] for row in rows
                ),
                "native_contact_fraction_median_across_replicas": (
                    statistics.median(
                        row["native_contact_fraction_median"]
                        for row in rows
                    )
                ),
                "buried_sasa_median_across_replicas_a2": (
                    statistics.median(
                        row["buried_sasa_median_a2"] for row in rows
                    )
                ),
                "mmpbsa_mean_kcal_mol": (
                    statistics.mean(mmpbsa_values)
                    if mmpbsa_values
                    else None
                ),
                "mmpbsa_range_kcal_mol": (
                    max(mmpbsa_values) - min(mmpbsa_values)
                    if len(mmpbsa_values) >= 2
                    else None
                ),
            }
        )

    candidate_csv = args.output_dir / "candidate_summary.csv"
    if candidates:
        with candidate_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(candidates[0]))
            writer.writeheader()
            writer.writerows(candidates)

    robust = [
        row["rank"]
        for row in candidates
        if row["three_replica_consensus"] == "robust"
    ]
    majority = [
        row["rank"]
        for row in candidates
        if row["three_replica_consensus"] == "majority"
    ]
    campaign = {
        "status": (
            "completed"
            if len(replicas) == 30 and len(candidates) == 10
            else "partial"
        ),
        "completed_replicas": len(replicas),
        "expected_replicas": 30,
        "completed_candidates": sum(
            row["completed_replicas"] == 3 for row in candidates
        ),
        "total_sampling_ns": sum(
            float(row["production_ns"] or 0.0) for row in replicas
        ),
        "robust_candidates": robust,
        "majority_candidates": majority,
        "predeclared_stability_gate": {
            "binder_rmsd_le_3a_frame_fraction_min": 0.70,
            "hotspots_occupancy_ge_0p5_min_count": 4,
            "buried_sasa_median_a2_min": 600.0,
        },
        "free_energy_interpretation": (
            "MM/GBSA is a relative endpoint diagnostic without entropy; "
            "it is not an absolute binding free energy and cannot be "
            "converted directly to KD."
        ),
        "experimental_interpretation": (
            "MD stability does not establish experimental binding, "
            "selectivity, inhibition, or cellular activity."
        ),
    }
    (args.output_dir / "campaign_summary.json").write_text(
        json.dumps(campaign, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(campaign, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
