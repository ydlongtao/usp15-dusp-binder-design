#!/usr/bin/env python3
"""Repair interrupted sparse Hugging Face downloads with audited HTTP ranges."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import time

import requests


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True, type=Path)
    parser.add_argument("--url", required=True)
    parser.add_argument("--expected-size", required=True, type=int)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def sparse_holes(file_descriptor: int, size: int) -> list[tuple[int, int]]:
    holes = []
    position = 0
    while position < size:
        try:
            data_start = os.lseek(file_descriptor, position, os.SEEK_DATA)
        except OSError:
            data_start = size
        if data_start > position:
            holes.append((position, data_start))
        if data_start >= size:
            break
        try:
            hole_start = os.lseek(file_descriptor, data_start, os.SEEK_HOLE)
        except OSError:
            hole_start = size
        position = hole_start
    return holes


def write_all(file_descriptor: int, data: bytes, offset: int) -> None:
    written = 0
    while written < len(data):
        count = os.pwrite(file_descriptor, data[written:], offset + written)
        if count <= 0:
            raise OSError(f"Short pwrite at offset {offset + written}")
        written += count


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def split_ranges(
    holes: list[tuple[int, int]], chunk_size: int
) -> list[tuple[int, int]]:
    ranges = []
    for start, end_exclusive in holes:
        position = start
        while position < end_exclusive:
            chunk_end = min(position + chunk_size, end_exclusive)
            ranges.append((position, chunk_end))
            position = chunk_end
    return ranges


def repair_range(
    *,
    file_descriptor: int,
    url: str,
    start: int,
    end_exclusive: int,
    expected_size: int,
    retries: int,
) -> dict[str, object]:
    end_inclusive = end_exclusive - 1
    expected_content_range = f"bytes {start}-{end_inclusive}/{expected_size}"
    expected_bytes = end_exclusive - start
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(
                url,
                headers={"Range": f"bytes={start}-{end_inclusive}"},
                timeout=(30, 300),
            )
            response.raise_for_status()
            observed_content_range = response.headers.get("Content-Range")
            if response.status_code != 206:
                raise ValueError(
                    f"Range {start}-{end_inclusive}: HTTP {response.status_code}"
                )
            if observed_content_range != expected_content_range:
                raise ValueError(
                    f"Range {start}-{end_inclusive}: expected Content-Range "
                    f"{expected_content_range!r}, found {observed_content_range!r}"
                )
            if len(response.content) != expected_bytes:
                raise ValueError(
                    f"Range {start}-{end_inclusive}: expected {expected_bytes} bytes, "
                    f"found {len(response.content)}"
                )
            write_all(file_descriptor, response.content, start)
            return {
                "start": start,
                "end_inclusive": end_inclusive,
                "bytes": expected_bytes,
                "content_range": observed_content_range,
                "attempt": attempt,
            }
        except (requests.RequestException, ValueError, OSError) as error:
            last_error = error
            if attempt < retries:
                time.sleep(2**attempt)
    raise RuntimeError(
        f"Range {start}-{end_inclusive} failed after {retries} attempts"
    ) from last_error


def main() -> None:
    args = parse_args()
    if args.file.stat().st_size != args.expected_size:
        raise ValueError(
            f"Expected logical size {args.expected_size}, found {args.file.stat().st_size}"
        )

    file_descriptor = os.open(args.file, os.O_RDWR)
    try:
        holes = sparse_holes(file_descriptor, args.expected_size)
        ranges = split_ranges(holes, args.chunk_size)
        repair_records = []
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    repair_range,
                    file_descriptor=file_descriptor,
                    url=args.url,
                    start=start,
                    end_exclusive=end_exclusive,
                    expected_size=args.expected_size,
                    retries=args.retries,
                )
                for start, end_exclusive in ranges
            ]
            for future in as_completed(futures):
                repair_records.append(future.result())
        os.fsync(file_descriptor)
    finally:
        os.close(file_descriptor)

    observed_sha256 = sha256(args.file)
    verified = observed_sha256 == args.expected_sha256
    report = {
        "file": str(args.file),
        "url": args.url,
        "expected_size": args.expected_size,
        "expected_sha256": args.expected_sha256,
        "observed_sha256": observed_sha256,
        "holes_repaired": len(repair_records),
        "bytes_repaired": sum(record["bytes"] for record in repair_records),
        "ranges": sorted(repair_records, key=lambda record: record["start"]),
        "chunk_size": args.chunk_size,
        "workers": args.workers,
        "verified": verified,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    if not verified:
        raise ValueError(
            f"SHA-256 mismatch: expected {args.expected_sha256}, found {observed_sha256}"
        )


if __name__ == "__main__":
    main()
