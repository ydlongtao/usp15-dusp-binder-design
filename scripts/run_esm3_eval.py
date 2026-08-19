#!/usr/bin/env python3
"""OVO-facing ESM3 descriptor/refolding adapter.

The adapter is intentionally orthogonal to AF2 and selectivity gates. It
accepts one FASTA record, runs ESM3 sequence and structure tracks, and writes
portable JSON/CSV/PDB outputs suitable for an OVO/Nextflow result processor.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
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
    ap.add_argument("--candidate-id", default=None)
    ap.add_argument("--model", default="esm3-sm-open-v1")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--num-steps", type=int, default=1)
    ap.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    ap.add_argument("--backbone-pdb", type=Path, default=None)
    ap.add_argument("--af2-pdb", type=Path, default=None)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    name, sequence = first_fasta(args.fasta)
    candidate_id = args.candidate_id or name
    result: dict[str, object] = {
        "schema_version": "esm3_ovo_descriptor_v1",
        "candidate_id": candidate_id,
        "input_id": name,
        "model": args.model,
        "seed": args.seed,
        "sequence_length": len(sequence),
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "python_platform": platform.platform(),
        "checks": {"fasta": "passed"},
        "status": "blocked",
    }
    started = time.time()
    try:
        import torch
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESMProtein, GenerationConfig

        if args.device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but unavailable")
        device = torch.device(args.device)
        result["device"] = str(device)
        result["torch_version"] = torch.__version__
        result["esm_version"] = __import__("esm").__version__
        if device.type == "cuda":
            result["cuda_device"] = torch.cuda.get_device_name(0)
            result["cuda_total_memory_bytes"] = torch.cuda.get_device_properties(0).total_memory
        torch.manual_seed(args.seed)
        model = ESM3.from_pretrained(args.model.replace("-", "_"), device=device)
        model.eval()
        result["checks"]["model_load"] = "passed"  # type: ignore[index]
        protein = ESMProtein(sequence=sequence)
        steps = max(1, min(args.num_steps, len(sequence)))
        with torch.no_grad():
            seq_start = time.time()
            seq_generated = model.generate(protein, GenerationConfig(track="sequence", num_steps=steps))
            result["sequence_generation_seconds"] = time.time() - seq_start
            result["sequence_track"] = getattr(seq_generated, "sequence", None)
            result["checks"]["sequence_track"] = "passed"  # type: ignore[index]
            structure_start = time.time()
            struct_generated = model.generate(protein, GenerationConfig(track="structure", num_steps=steps))
            struct_generated.to_pdb(args.output_dir / "esm3_structure.pdb")
            result["structure_generation_seconds"] = time.time() - structure_start
            result["checks"]["structure_track"] = "passed"  # type: ignore[index]
        result["status"] = "passed"
    except (Exception, KeyboardInterrupt) as exc:  # noqa: BLE001
        result["error"] = f"{type(exc).__name__}: {exc}"
    result["elapsed_seconds"] = time.time() - started
    result["backbone_pdb_supplied"] = bool(args.backbone_pdb)
    result["af2_pdb_supplied"] = bool(args.af2_pdb)
    (args.output_dir / "esm3_eval.json").write_text(json.dumps(result, indent=2))
    with (args.output_dir / "esm3_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "status", "model", "device", "sequence_length", "sequence_track", "structure_track", "elapsed_seconds"])
        writer.writeheader()
        checks = result.get("checks", {})
        writer.writerow({
            "candidate_id": candidate_id,
            "status": result.get("status"),
            "model": args.model,
            "device": result.get("device"),
            "sequence_length": len(sequence),
            "sequence_track": checks.get("sequence_track", "blocked") if isinstance(checks, dict) else "blocked",
            "structure_track": checks.get("structure_track", "blocked") if isinstance(checks, dict) else "blocked",
            "elapsed_seconds": result.get("elapsed_seconds"),
        })
    return 0 if result["status"] == "passed" else 4


if __name__ == "__main__":
    raise SystemExit(main())
