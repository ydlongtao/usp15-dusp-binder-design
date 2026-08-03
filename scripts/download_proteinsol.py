#!/usr/bin/env python3
"""Download and verify the official Protein-Sol sequence software archive."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import requests


URL = (
    "https://protein-sol.manchester.ac.uk/cgi-bin/utilities/"
    "download_sequence_code.php"
)
EXPECTED_BYTES = 34958
EXPECTED_SHA256 = (
    "4df32c61fca53adcb2394a528babd1ad85cb5c551bf7bd1c56d134097fb2b1b8"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if not args.output.is_file():
        response = requests.get(URL, timeout=120)
        response.raise_for_status()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(response.content)
    observed_bytes = args.output.stat().st_size
    observed_sha256 = sha256(args.output)
    verified = (
        observed_bytes == EXPECTED_BYTES
        and observed_sha256 == EXPECTED_SHA256
    )
    report = {
        "url": URL,
        "path": str(args.output),
        "expected_bytes": EXPECTED_BYTES,
        "observed_bytes": observed_bytes,
        "expected_sha256": EXPECTED_SHA256,
        "observed_sha256": observed_sha256,
        "verified": verified,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    if not verified:
        raise SystemExit("Protein-Sol archive integrity verification failed")


if __name__ == "__main__":
    main()
