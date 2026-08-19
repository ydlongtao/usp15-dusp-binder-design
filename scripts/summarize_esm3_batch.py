#!/usr/bin/env python3
"""Aggregate ESM3 descriptor JSON files without promoting candidates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_dir", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    records = []
    for path in sorted(args.input_dir.glob("*/esm3_eval.json")):
        records.append(json.loads(path.read_text()))
    summary = {
        "schema_version": "esm3_ovo_batch_summary_v1",
        "candidate_count": len(records),
        "passed_count": sum(r.get("status") == "passed" for r in records),
        "sequence_track_passed": sum(r.get("checks", {}).get("sequence_track") == "passed" for r in records),
        "structure_track_passed": sum(r.get("checks", {}).get("structure_track") == "passed" for r in records),
        "records": records,
        "promotion": "not_evaluated: ESM3 is an orthogonal descriptor, not an AF2/selectivity gate",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in summary if k != "records"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
