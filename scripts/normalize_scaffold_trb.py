#!/usr/bin/env python3
"""Normalize auto-contig scaffold TRB metadata for OVO standardization."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_trb_dir", type=Path)
    parser.add_argument("normalized_trb_dir", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    raw_paths = sorted(args.raw_trb_dir.glob("*.trb"))
    if not raw_paths:
        raise ValueError(f"No TRB files found in {args.raw_trb_dir}")

    args.normalized_trb_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for raw_path in raw_paths:
        with raw_path.open("rb") as handle:
            trb = pickle.load(handle)

        sampled_mask = trb.get("sampled_mask")
        if not sampled_mask or not all(
            isinstance(item, str) and item for item in sampled_mask
        ):
            raise ValueError(f"Missing valid sampled_mask in {raw_path}")

        config = trb.get("config")
        if not isinstance(config, dict):
            raise ValueError(f"Missing config dictionary in {raw_path}")
        contigmap = config.get("contigmap")
        if not isinstance(contigmap, dict):
            raise ValueError(f"Missing config.contigmap in {raw_path}")

        original_contigs = contigmap.get("contigs")
        if original_contigs not in (None, list(sampled_mask)):
            raise ValueError(
                f"Unexpected pre-existing contigs in {raw_path}: "
                f"{original_contigs!r}"
            )

        normalized_contigs = list(sampled_mask)
        contigmap["contigs"] = normalized_contigs
        trb["ovo_scaffold_metadata_normalization"] = {
            "version": 1,
            "source_field": "sampled_mask",
            "original_config_contigs": original_contigs,
            "normalized_config_contigs": normalized_contigs,
            "coordinates_modified": False,
        }

        output_path = args.normalized_trb_dir / raw_path.name
        with output_path.open("wb") as handle:
            pickle.dump(trb, handle, protocol=pickle.HIGHEST_PROTOCOL)
        records.append(
            {
                "trb": raw_path.name,
                "original_config_contigs": original_contigs,
                "normalized_config_contigs": normalized_contigs,
                "coordinates_modified": False,
            }
        )

    report = {
        "status": "completed",
        "normalized_trb_count": len(records),
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
