#!/usr/bin/env python3
"""Apply paired-seed R10 USP4/USP11 same-pose selectivity gates."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


POSITIVE_GATES = {
    "ipae_max": 10.0,
    "target_aligned_binder_rmsd_max": 2.0,
    "binder_plddt_min": 80.0,
}
SELECTIVITY_GATES = {
    "delta_ipae_min": 5.0,
    "off_target_ipae_min": 15.0,
    "off_target_binder_rmsd_min_exclusive": 4.0,
    "passing_seed_count_min": 2,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-jsonl", required=True, type=Path)
    parser.add_argument("--positive-summary", required=True, type=Path)
    parser.add_argument("--input-report", required=True, type=Path)
    parser.add_argument("--usp4-jsonl", required=True, type=Path)
    parser.add_argument("--usp11-jsonl", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def read_jsonl(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def index_records(records: list[dict], expected_ids: set[str]) -> dict:
    indexed = {}
    for record in records:
        if record.get("model_name") != "model_2_ptm":
            raise ValueError(f"Unexpected model {record.get('model_name')}")
        if record.get("template_mode") != "ct":
            raise ValueError(f"Unexpected template mode {record.get('template_mode')}")
        if record.get("num_recycles") != 3 or not record.get("dropout"):
            raise ValueError("R10 selectivity record violates AF2 invariants")
        candidate_id = record["input_id"]
        seed = int(record["seed"])
        if candidate_id in expected_ids:
            indexed[(candidate_id, seed)] = record
    expected_keys = {
        (candidate_id, seed)
        for candidate_id in expected_ids
        for seed in (0, 1, 2)
    }
    if set(indexed) != expected_keys:
        missing = sorted(expected_keys - set(indexed))
        extra = sorted(set(indexed) - expected_keys)
        raise ValueError(f"Record key mismatch; missing={missing[:5]} extra={extra[:5]}")
    return indexed


def finite(record: dict, key: str) -> float:
    value = float(record[key])
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {key} in {record.get('id')}")
    return value


def positive_seed(record: dict) -> bool:
    return (
        finite(record, "ipae") <= POSITIVE_GATES["ipae_max"]
        and finite(record, "target_aligned_binder_rmsd")
        <= POSITIVE_GATES["target_aligned_binder_rmsd_max"]
        and finite(record, "binder_plddt")
        >= POSITIVE_GATES["binder_plddt_min"]
    )


def main() -> None:
    args = parse_args()
    positive_summary = json.loads(args.positive_summary.read_text(encoding="utf-8"))
    input_report = json.loads(args.input_report.read_text(encoding="utf-8"))
    candidate_ids = set(positive_summary["passing_ids"])
    report_ids = {record["id"] for record in input_report["records"]}
    if not candidate_ids or candidate_ids != report_ids:
        raise ValueError("Positive summary and selectivity input report disagree")

    on_target = index_records(read_jsonl(args.positive_jsonl), candidate_ids)
    off_targets = {
        "USP4": index_records(read_jsonl(args.usp4_jsonl), candidate_ids),
        "USP11": index_records(read_jsonl(args.usp11_jsonl), candidate_ids),
    }

    rows = []
    details = []
    for candidate_id in sorted(candidate_ids):
        candidate_detail = {"id": candidate_id, "homologs": {}}
        row = {"id": candidate_id}
        final_pass = True
        for homolog, records in off_targets.items():
            seed_details = []
            passing_seeds = []
            for seed in (0, 1, 2):
                on_record = on_target[(candidate_id, seed)]
                off_record = records[(candidate_id, seed)]
                on_ipae = finite(on_record, "ipae")
                off_ipae = finite(off_record, "ipae")
                off_rmsd = finite(off_record, "target_aligned_binder_rmsd")
                delta_ipae = off_ipae - on_ipae
                passed = (
                    positive_seed(on_record)
                    and delta_ipae >= SELECTIVITY_GATES["delta_ipae_min"]
                    and (
                        off_ipae >= SELECTIVITY_GATES["off_target_ipae_min"]
                        or off_rmsd
                        > SELECTIVITY_GATES[
                            "off_target_binder_rmsd_min_exclusive"
                        ]
                    )
                )
                if passed:
                    passing_seeds.append(seed)
                seed_details.append(
                    {
                        "seed": seed,
                        "on_target_ipae": on_ipae,
                        "off_target_ipae": off_ipae,
                        "delta_ipae": delta_ipae,
                        "off_target_binder_rmsd": off_rmsd,
                        "off_target_binder_plddt": finite(
                            off_record, "binder_plddt"
                        ),
                        "pass": passed,
                    }
                )
            homolog_pass = (
                len(passing_seeds)
                >= SELECTIVITY_GATES["passing_seed_count_min"]
            )
            final_pass = final_pass and homolog_pass
            candidate_detail["homologs"][homolog] = {
                "passing_seeds": passing_seeds,
                "pass": homolog_pass,
                "seeds": seed_details,
            }
            row[f"{homolog.lower()}_passing_seed_count"] = len(passing_seeds)
            row[f"{homolog.lower()}_passing_seeds"] = ",".join(
                map(str, passing_seeds)
            )
            row[f"{homolog.lower()}_pass"] = homolog_pass
        candidate_detail["selectivity_pass"] = final_pass
        row["selectivity_pass"] = final_pass
        details.append(candidate_detail)
        rows.append(row)

    passing_ids = [row["id"] for row in rows if row["selectivity_pass"]]
    summary = {
        "phase": "R10 same-pose USP4/USP11 selectivity screen",
        "model": "model_2_ptm",
        "template_mode": "ct",
        "dropout": True,
        "seeds": [0, 1, 2],
        "recycles": 3,
        "positive_gates": POSITIVE_GATES,
        "selectivity_gates": SELECTIVITY_GATES,
        "input_count": len(candidate_ids),
        "passing_count": len(passing_ids),
        "passing_ids": passing_ids,
        "rows": rows,
        "details": details,
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
    print(json.dumps({"input_count": len(rows), "passing_count": len(passing_ids)}))


if __name__ == "__main__":
    main()
