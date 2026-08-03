#!/usr/bin/env python3
"""Download and audit official MIT-licensed Boltz-2 inference assets."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from huggingface_hub import hf_hub_download


REPO_ID = "boltz-community/boltz-2"
EXPECTED_SIZES = {
    "mols.tar": 1_855_662_080,
    "boltz2_conf.ckpt": 2_286_561_469,
    "boltz2_aff.ckpt": 2_062_139_170,
}
EXPECTED_SHA256 = {
    "mols.tar": "39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7",
    "boltz2_conf.ckpt": (
        "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
    ),
    "boltz2_aff.ckpt": (
        "dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
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
    records = []
    for filename, expected_size in EXPECTED_SIZES.items():
        downloaded = Path(
            hf_hub_download(
                repo_id=REPO_ID,
                filename=filename,
                local_dir=args.output_dir,
            )
        )
        observed_size = downloaded.stat().st_size
        if observed_size != expected_size:
            raise ValueError(
                f"{filename}: expected {expected_size} bytes, found {observed_size}"
            )
        observed_sha256 = sha256(downloaded)
        if observed_sha256 != EXPECTED_SHA256[filename]:
            raise ValueError(
                f"{filename}: expected SHA-256 {EXPECTED_SHA256[filename]}, "
                f"found {observed_sha256}"
            )
        records.append(
            {
                "repo_id": REPO_ID,
                "filename": filename,
                "bytes": observed_size,
                "sha256": observed_sha256,
                "license": "MIT",
            }
        )

    manifest = {
        "phase": "R8 Boltz-2 official asset download",
        "source": f"https://huggingface.co/{REPO_ID}",
        "records": records,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
