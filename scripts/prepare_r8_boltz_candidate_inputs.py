#!/usr/bin/env python3
"""Create sequence-only Boltz YAMLs for AF2-passing R8 candidates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from prepare_r8_boltz_controls import sequences, yaml_text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--af2-summary", required=True, type=Path)
    parser.add_argument("--pdb-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = json.loads(args.af2_summary.read_text(encoding="utf-8"))
    passing_ids = summary.get("passing_ids", [])
    if not passing_ids:
        raise ValueError("No AF2-positive candidates for Boltz-2 screening")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for candidate_id in passing_ids:
        source = args.pdb_dir / f"{candidate_id}.pdb"
        if not source.is_file():
            raise ValueError(f"Missing candidate PDB: {source}")
        chain_sequences = sequences(source)
        if not 45 <= len(chain_sequences["A"]) <= 80:
            raise ValueError(f"{candidate_id}: invalid binder length")
        if len(chain_sequences["B"]) != 129:
            raise ValueError(f"{candidate_id}: incomplete USP15 target")
        if "C" in chain_sequences["A"]:
            raise ValueError(f"{candidate_id}: binder contains Cys")
        yaml_path = args.output_dir / f"{candidate_id}.yaml"
        yaml_path.write_text(yaml_text(chain_sequences), encoding="utf-8")
        records.append(
            {
                "id": candidate_id,
                "source_pdb": str(source),
                "yaml": str(yaml_path),
                "binder_length": len(chain_sequences["A"]),
                "target_length": len(chain_sequences["B"]),
            }
        )
    report = {
        "phase": "R8 Boltz-2 candidate input preparation",
        "templates": None,
        "forced_constraints": False,
        "input_count": len(records),
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"input_count": len(records)}))


if __name__ == "__main__":
    main()
