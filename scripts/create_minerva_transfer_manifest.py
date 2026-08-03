#!/usr/bin/env python3
"""Create or verify SHA-256 provenance for a transferred MD directory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXCLUDED = {".queue.lock", ".DS_Store", "transfer_manifest.json"}


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inventory(root: Path):
    rows = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if path.name in EXCLUDED:
            continue
        rows.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--md-dir", required=True, type=Path)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("transfer_manifest.json"),
    )
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    manifest = (
        args.manifest
        if args.manifest.is_absolute()
        else args.md_dir / args.manifest
    )
    if args.verify:
        expected = json.loads(manifest.read_text())
        observed = inventory(args.md_dir)
        if observed != expected["files"]:
            raise SystemExit("Transfer manifest verification failed")
        print(
            json.dumps(
                {
                    "status": "passed",
                    "files": len(observed),
                    "bytes": sum(row["bytes"] for row in observed),
                },
                indent=2,
            )
        )
        return

    rows = inventory(args.md_dir)
    report = {
        "status": "created",
        "files": rows,
        "file_count": len(rows),
        "total_bytes": sum(row["bytes"] for row in rows),
    }
    manifest.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": "created",
                "manifest": str(manifest),
                "files": len(rows),
                "bytes": report["total_bytes"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
