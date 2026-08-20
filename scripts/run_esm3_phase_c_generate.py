#!/usr/bin/env python3
"""Generate ESM3 sequence variants from a fixed single-chain binder backbone.

This is a sequence-generation pilot only. It does not model the USP15 complex
and its outputs must go through the existing AF2 and selectivity gates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone-pdb", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--model", default="esm3-sm-open-v1")
    ap.add_argument("--num-samples", type=int, default=3)
    ap.add_argument("--num-steps", type=int, default=20)
    ap.add_argument("--temperature", type=float, default=0.8)
    ap.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    result = {
        "schema_version": "esm3_phase_c_generation_v1",
        "backbone_pdb": args.backbone_pdb.name,
        "model": args.model,
        "num_samples": args.num_samples,
        "num_steps": args.num_steps,
        "temperature": args.temperature,
        "status": "blocked",
        "checks": {"backbone_pdb": "passed"},
        "samples": [],
        "scope": "fixed single-chain binder backbone; no USP15 complex conditioning",
    }
    try:
        import torch
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESMProtein, GenerationConfig

        if args.device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but unavailable")
        device = torch.device(args.device)
        protein = ESMProtein.from_pdb(args.backbone_pdb)
        if protein.coordinates is None:
            raise ValueError("backbone PDB did not yield coordinates")
        sequence_length = int(protein.coordinates.shape[0])
        # Sequence generation conditioned on coordinates: mask the input sequence.
        protein.sequence = None
        torch.manual_seed(0)
        model = ESM3.from_pretrained(args.model.replace("-", "_"), device=device)
        model.eval()
        result["device"] = str(device)
        result["sequence_length"] = sequence_length
        result["esm_version"] = __import__("esm").__version__
        result["torch_version"] = torch.__version__
        for sample_id in range(max(1, args.num_samples)):
            torch.manual_seed(sample_id)
            cfg = GenerationConfig(
                track="sequence",
                num_steps=max(1, min(args.num_steps, sequence_length)),
                temperature=args.temperature,
                condition_on_coordinates_only=True,
            )
            with torch.no_grad():
                generated = model.generate(protein, cfg)
            sequence = generated.sequence
            if not sequence or len(sequence) != sequence_length:
                raise ValueError(f"invalid generated sequence length: {len(sequence or '')}")
            out_fasta = args.output_dir / f"sample_{sample_id:02d}.fasta"
            out_fasta.write_text(f">esm3_phase_c_sample_{sample_id:02d}\n{sequence}\n")
            result["samples"].append({
                "sample_id": sample_id,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "fasta": out_fasta.name,
                "has_cysteine": "C" in sequence,
            })
        result["checks"]["model_load"] = "passed"
        result["checks"]["sequence_generation"] = "passed"
        result["status"] = "passed"
    except Exception as exc:  # noqa: BLE001
        result["error"] = f"{type(exc).__name__}: {exc}"
    result["elapsed_seconds"] = time.time() - started
    (args.output_dir / "phase_c_generation.json").write_text(json.dumps(result, indent=2))
    if result["samples"]:
        with (args.output_dir / "phase_c_sequences.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample_id", "sequence_sha256", "fasta", "has_cysteine"])
            writer.writeheader()
            writer.writerows(result["samples"])
    return 0 if result["status"] == "passed" else 4


if __name__ == "__main__":
    raise SystemExit(main())
