#!/usr/bin/env python3
"""Download the public ESM-IF1 checkpoint with resumable HTTP and an audit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import requests


URL = (
    "https://dl.fbaipublicfiles.com/fair-esm/models/"
    "esm_if1_gvp4_t16_142M_UR50.pt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    with requests.get(URL, headers=headers, stream=True, timeout=120) as response:
        if existing and response.status_code != 206:
            existing = 0
            partial.unlink(missing_ok=True)
        response.raise_for_status()
        mode = "ab" if existing else "wb"
        with partial.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    handle.write(chunk)
    partial.replace(args.output)
    report = {
        "url": URL,
        "path": str(args.output),
        "bytes": args.output.stat().st_size,
        "sha256": sha256(args.output),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report))


if __name__ == "__main__":
    main()
