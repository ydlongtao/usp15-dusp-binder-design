#!/usr/bin/env python3
"""Run a bounded, standalone ESM3 Phase-A smoke test.

This script deliberately does not import OVO or launch a design workflow. It
records runtime/model failures as structured JSON rather than claiming success.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path


def first_fasta(path: Path) -> tuple[str, str]:
    name = None
    seq_parts: list[str] = []
    records: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq_parts)))
            name = line[1:].strip() or "candidate_1"
            seq_parts = []
        elif name is not None:
            seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts)))
    if not records or not records[0][1]:
        raise ValueError(f"no FASTA sequence found in {path}")
    return records[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--model", default="esm3-sm-open-v1")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    ap.add_argument("--num-steps", type=int, default=1)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {
        "phase": "A",
        "status": "blocked",
        "model": args.model,
        "seed": args.seed,
        "python": sys.version,
        "platform": platform.platform(),
        "input_fasta": str(args.fasta),
        "started_at_epoch": time.time(),
        "checks": {},
    }
    try:
        name, sequence = first_fasta(args.fasta)
        summary["input_id"] = name
        summary["sequence_length"] = len(sequence)
        summary["sequence_sha256"] = hashlib.sha256(sequence.encode()).hexdigest()
        checks = summary["checks"]
        assert isinstance(checks, dict)
        checks["fasta"] = "passed"
    except Exception as exc:  # noqa: BLE001
        summary["error"] = f"FASTA: {type(exc).__name__}: {exc}"
        (args.output_dir / "phase_a_summary.json").write_text(json.dumps(summary, indent=2))
        return 2

    if importlib.util.find_spec("torch") is None:
        summary["error"] = "PyTorch is not installed; install a compatible torch/ESM3 runtime first."
        (args.output_dir / "phase_a_summary.json").write_text(json.dumps(summary, indent=2))
        return 3
    if importlib.util.find_spec("esm") is None:
        summary["error"] = "The ESM package is not installed; install the official ESM3 package first."
        (args.output_dir / "phase_a_summary.json").write_text(json.dumps(summary, indent=2))
        return 3

    try:
        import torch
        import esm
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESMProtein, GenerationConfig

        summary["torch_version"] = torch.__version__
        summary["esm_version"] = getattr(esm, "__version__", "unknown")
        if args.device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is not available")
        if args.device == "cuda" or (args.device == "auto" and torch.cuda.is_available()):
            device = "cuda"
        else:
            device = "cpu"
        summary["device"] = device
        if device == "cuda":
            summary["cuda_device"] = torch.cuda.get_device_name(0)
            summary["cuda_total_memory_bytes"] = torch.cuda.get_device_properties(0).total_memory

        torch.manual_seed(args.seed)
        model_name = args.model.replace("-", "_")
        model = ESM3.from_pretrained(model_name, device=torch.device(device))
        model.eval()
        checks = summary["checks"]
        assert isinstance(checks, dict)
        checks["model_load"] = "passed"

        protein = ESMProtein(sequence=sequence)
        generation = GenerationConfig(track="sequence", num_steps=max(1, min(args.num_steps, len(sequence))))
        start = time.time()
        with torch.no_grad():
            generated = model.generate(protein, generation)
        (args.output_dir / "sequence_generation.json").write_text(json.dumps({"input_id": name, "result": str(generated)}, indent=2))
        summary["sequence_generation_seconds"] = time.time() - start
        checks["sequence_track"] = "passed"
        summary["status"] = "passed"
    except (Exception, KeyboardInterrupt) as exc:  # noqa: BLE001
        summary["error"] = f"runtime: {type(exc).__name__}: {exc}"
    finally:
        summary["finished_at_epoch"] = time.time()
        (args.output_dir / "phase_a_summary.json").write_text(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "passed" else 4


if __name__ == "__main__":
    raise SystemExit(main())
