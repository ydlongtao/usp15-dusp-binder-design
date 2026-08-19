#!/usr/bin/env python3
"""Split a multi-record FASTA into deterministic per-candidate files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("output_dir", type=Path)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    name = None
    seq: list[str] = []
    count = 0
    for line in args.input.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                count += 1
                safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or f"candidate_{count}"
                (args.output_dir / f"{safe}.fasta").write_text(f">{name}\n{''.join(seq)}\n")
            name, seq = line[1:].strip(), []
        elif name is not None:
            seq.append(line)
    if name is not None:
        count += 1
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or f"candidate_{count}"
        (args.output_dir / f"{safe}.fasta").write_text(f">{name}\n{''.join(seq)}\n")
    if count == 0:
        raise SystemExit("No FASTA records found")
    print(count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
