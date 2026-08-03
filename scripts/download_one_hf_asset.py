#!/usr/bin/env python3
"""Download one Hugging Face asset and verify size plus LFS SHA-256."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from huggingface_hub import hf_hub_download


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-id", required=True)
    parser.add_argument("--filename", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-size", required=True, type=int)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = Path(
        hf_hub_download(
            repo_id=args.repo_id,
            filename=args.filename,
            local_dir=args.output_dir,
        )
    )
    observed_size = path.stat().st_size
    observed_sha256 = sha256(path)
    verified = (
        observed_size == args.expected_size
        and observed_sha256 == args.expected_sha256
    )
    report = {
        "repo_id": args.repo_id,
        "filename": args.filename,
        "path": str(path),
        "expected_size": args.expected_size,
        "observed_size": observed_size,
        "expected_sha256": args.expected_sha256,
        "observed_sha256": observed_sha256,
        "verified": verified,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    if not verified:
        raise ValueError(f"Asset verification failed for {args.filename}")


if __name__ == "__main__":
    main()
