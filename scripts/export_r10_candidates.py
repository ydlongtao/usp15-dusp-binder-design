#!/usr/bin/env python3
"""De-redundant and export complete R10 geometry-conditioned candidates."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

from Bio.Align import PairwiseAligner


IDENTITY_CLUSTER_THRESHOLD = 0.80
MAX_EXPORTS = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-summary", required=True, type=Path)
    parser.add_argument("--interface-summary", required=True, type=Path)
    parser.add_argument("--selectivity-summary", required=True, type=Path)
    parser.add_argument("--input-report", required=True, type=Path)
    parser.add_argument("--panel-dir", required=True, type=Path)
    parser.add_argument("--positive-output-dir", required=True, type=Path)
    parser.add_argument("--selectivity-output-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sequence_identity(sequence_a: str, sequence_b: str) -> float:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.1
    alignment = aligner.align(sequence_a, sequence_b)[0]
    matches = 0
    aligned_length = 0
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        block_length = min(a_end - a_start, b_end - b_start)
        matches += sum(
            sequence_a[a_start + offset] == sequence_b[b_start + offset]
            for offset in range(block_length)
        )
        aligned_length += block_length
    denominator = max(len(sequence_a), len(sequence_b), aligned_length)
    return matches / denominator


def positive_aggregates(row: dict) -> dict:
    seeds = [int(seed) for seed in str(row["passing_seeds"]).split(",") if seed != ""]
    return {
        "positive_passing_seed_count": len(seeds),
        "positive_mean_ipae": sum(float(row[f"seed{seed}_ipae"]) for seed in seeds)
        / len(seeds),
        "positive_mean_binder_rmsd": sum(
            float(row[f"seed{seed}_binder_rmsd"]) for seed in seeds
        )
        / len(seeds),
        "positive_mean_binder_plddt": sum(
            float(row[f"seed{seed}_binder_plddt"]) for seed in seeds
        )
        / len(seeds),
        "best_positive_seed": min(
            seeds,
            key=lambda seed: (
                float(row[f"seed{seed}_ipae"]),
                float(row[f"seed{seed}_binder_rmsd"]),
                -float(row[f"seed{seed}_binder_plddt"]),
            ),
        ),
    }


def copy_prediction_set(
    candidate_id: str,
    source_dir: Path,
    destination_dir: Path,
) -> list[str]:
    copied = []
    for seed in (0, 1, 2):
        source = (
            source_dir
            / f"{candidate_id}__seed{seed}__model_2_ptm_ct.pdb"
        )
        if not source.is_file():
            raise FileNotFoundError(source)
        destination_dir.mkdir(parents=True, exist_ok=True)
        destination = destination_dir / source.name
        shutil.copy2(source, destination)
        copied.append(str(destination))
    return copied


def main() -> None:
    args = parse_args()
    positive = json.loads(args.positive_summary.read_text(encoding="utf-8"))
    interface = json.loads(args.interface_summary.read_text(encoding="utf-8"))
    selectivity = json.loads(
        args.selectivity_summary.read_text(encoding="utf-8")
    )
    inputs = json.loads(args.input_report.read_text(encoding="utf-8"))

    positive_rows = {row["id"]: row for row in positive["rows"]}
    interface_rows = {row["id"]: row for row in interface["rows"]}
    details = {record["id"]: record for record in selectivity["details"]}
    input_records = {record["id"]: record for record in inputs["records"]}
    passing_ids = list(selectivity["passing_ids"])
    if not passing_ids:
        raise ValueError("No R10 selectivity passers to export")

    ranked = []
    for candidate_id in passing_ids:
        aggregates = positive_aggregates(positive_rows[candidate_id])
        deltas = []
        selective_seed_count = 0
        for homolog in ("USP4", "USP11"):
            for seed_record in details[candidate_id]["homologs"][homolog]["seeds"]:
                if seed_record["pass"]:
                    selective_seed_count += 1
                    deltas.append(float(seed_record["delta_ipae"]))
        ranked.append(
            {
                "id": candidate_id,
                "sequence": input_records[candidate_id]["sequence"],
                "length": int(input_records[candidate_id]["length"]),
                "selective_seed_count": selective_seed_count,
                "min_passing_delta_ipae": min(deltas),
                "interface_delta_sasa_A2": float(
                    interface_rows[candidate_id]["interface_delta_sasa_A2"]
                ),
                **aggregates,
            }
        )
    ranked.sort(
        key=lambda row: (
            -row["selective_seed_count"],
            -row["min_passing_delta_ipae"],
            row["positive_mean_ipae"],
            row["positive_mean_binder_rmsd"],
            -row["positive_mean_binder_plddt"],
            -row["interface_delta_sasa_A2"],
            row["id"],
        )
    )

    representatives = []
    cluster_members: dict[str, list[str]] = {}
    assignment = {}
    for row in ranked:
        matched_rep = None
        matched_identity = 0.0
        for representative in representatives:
            identity = sequence_identity(row["sequence"], representative["sequence"])
            if identity >= IDENTITY_CLUSTER_THRESHOLD:
                matched_rep = representative
                matched_identity = identity
                break
        if matched_rep is None:
            representatives.append(row)
            cluster_members[row["id"]] = [row["id"]]
            assignment[row["id"]] = {
                "representative": row["id"],
                "identity": 1.0,
            }
        else:
            cluster_members[matched_rep["id"]].append(row["id"])
            assignment[row["id"]] = {
                "representative": matched_rep["id"],
                "identity": matched_identity,
            }

    selected = representatives[:MAX_EXPORTS]
    selected_ids = {row["id"] for row in selected}
    if len(selected) < 3:
        raise ValueError(
            f"Only {len(selected)} sequence clusters passed all R10 gates"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_lines = []
    export_records = []
    for rank, row in enumerate(selected, start=1):
        candidate_id = row["id"]
        candidate_dir = args.output_dir / f"rank_{rank:02d}_{candidate_id}"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        design_input = args.panel_dir / f"{candidate_id}.pdb"
        copied_input = candidate_dir / "design_complex.pdb"
        shutil.copy2(design_input, copied_input)
        positive_predictions = copy_prediction_set(
            candidate_id,
            args.positive_output_dir,
            candidate_dir / "predictions" / "USP15",
        )
        off_target_predictions = {}
        for homolog in ("USP4", "USP11"):
            off_target_predictions[homolog] = copy_prediction_set(
                candidate_id,
                args.selectivity_output_dir
                / homolog
                / "ptm_model_2_ct",
                candidate_dir / "predictions" / homolog,
            )
        fasta_lines.extend([f">{candidate_id}", row["sequence"]])
        export_records.append(
            {
                "rank": rank,
                **row,
                "cluster_member_ids": cluster_members[candidate_id],
                "design_complex_pdb": str(copied_input),
                "positive_prediction_pdbs": positive_predictions,
                "off_target_prediction_pdbs": off_target_predictions,
                "geometry_conditioned": True,
            }
        )

    (args.output_dir / "candidates.fasta").write_text(
        "\n".join(fasta_lines) + "\n",
        encoding="utf-8",
    )
    metrics_path = args.output_dir / "candidate_metrics.csv"
    with metrics_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            key
            for key in selected[0]
            if key not in {"sequence"}
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(
            {key: value for key, value in row.items() if key in fieldnames}
            for row in selected
        )

    elimination_rows = []
    for row in ranked:
        if row["id"] in selected_ids:
            reason = "selected_representative"
        else:
            reason = (
                "sequence_redundant"
                if assignment[row["id"]]["representative"] in selected_ids
                else "representative_below_export_rank"
            )
        elimination_rows.append(
            {
                "id": row["id"],
                "cluster_representative": assignment[row["id"]][
                    "representative"
                ],
                "identity_to_representative": assignment[row["id"]]["identity"],
                "disposition": reason,
            }
        )
    with (args.output_dir / "elimination_reasons.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=elimination_rows[0].keys())
        writer.writeheader()
        writer.writerows(elimination_rows)

    manifest = {
        "phase": "R10 final geometry-conditioned candidate export",
        "candidate_count": len(export_records),
        "sequence_identity_cluster_threshold": IDENTITY_CLUSTER_THRESHOLD,
        "max_exports": MAX_EXPORTS,
        "geometry_conditioned": True,
        "experimental_validation_required": True,
        "records": export_records,
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "selectivity_passing_count": len(ranked),
                "sequence_cluster_count": len(representatives),
                "exported_count": len(export_records),
            }
        )
    )


if __name__ == "__main__":
    main()
