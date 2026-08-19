#!/usr/bin/env python3
"""Run the bounded ESM3 structure-track smoke test."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

from run_esm3_phase_a_smoke import first_fasta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--model", default="esm3-sm-open-v1")
    ap.add_argument("--num-steps", type=int, default=1)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    name, sequence = first_fasta(args.fasta)
    summary: dict[str, object] = {
        "phase": "A",
        "track": "structure",
        "status": "blocked",
        "model": args.model,
        "input_id": name,
        "sequence_length": len(sequence),
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
        "checks": {"fasta": "passed"},
    }
    try:
        import torch
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESMProtein, GenerationConfig

        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is not available")
        summary["torch_version"] = torch.__version__
        summary["esm_version"] = __import__("esm").__version__
        summary["device"] = "cuda"
        summary["cuda_device"] = torch.cuda.get_device_name(0)
        summary["cuda_total_memory_bytes"] = torch.cuda.get_device_properties(0).total_memory
        model = ESM3.from_pretrained(args.model.replace("-", "_"), device=torch.device("cuda"))
        model.eval()
        summary["checks"]["model_load"] = "passed"  # type: ignore[index]
        protein = ESMProtein(sequence=sequence)
        config = GenerationConfig(track="structure", num_steps=max(1, min(args.num_steps, len(sequence))))
        started = time.time()
        with torch.no_grad():
            generated = model.generate(protein, config)
        generated.to_pdb(args.output_dir / "esm3_structure_smoke.pdb")
        summary["generation_seconds"] = time.time() - started
        summary["checks"]["structure_track"] = "passed"  # type: ignore[index]
        summary["status"] = "passed"
    except (Exception, KeyboardInterrupt) as exc:  # noqa: BLE001
        summary["error"] = f"runtime: {type(exc).__name__}: {exc}"
    summary["finished_at_epoch"] = time.time()
    (args.output_dir / "structure_smoke_summary.json").write_text(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "passed" else 4


if __name__ == "__main__":
    raise SystemExit(main())
