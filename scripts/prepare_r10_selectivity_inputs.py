#!/usr/bin/env python3
"""Build same-pose USP4/USP11 complexes for R10 positive-screen passers."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-summary", required=True, type=Path)
    parser.add_argument("--interface-summary", required=True, type=Path)
    parser.add_argument("--panel-dir", required=True, type=Path)
    parser.add_argument("--usp4-target", required=True, type=Path)
    parser.add_argument("--usp11-target", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atom_lines(path: Path, chain: str) -> list[str]:
    lines = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if (
            line.startswith("ATOM")
            and line[21] == chain
            and line[16] in (" ", "A")
        ):
            lines.append(line[:16] + " " + line[17:])
    if not lines:
        raise ValueError(f"No chain {chain} ATOM records in {path}")
    return lines


def sequence_from_binder(lines: list[str]) -> str:
    sequence = []
    seen = set()
    for line in lines:
        if line[12:16].strip() != "CA":
            continue
        key = (line[22:26], line[26])
        if key in seen:
            continue
        seen.add(key)
        residue = line[17:20].strip()
        if residue not in AA3_TO_1:
            raise ValueError(f"Unsupported binder residue {residue}")
        sequence.append(AA3_TO_1[residue])
    if not sequence:
        raise ValueError("Binder sequence is empty")
    return "".join(sequence)


def write_complex(binder: list[str], target: list[str], output: Path) -> None:
    serial = 1
    lines = []
    for source_line in binder + target:
        lines.append(f"{source_line[:6]}{serial:5d}{source_line[11:]}")
        serial += 1
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\nTER\nEND\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    summary = json.loads(args.positive_summary.read_text(encoding="utf-8"))
    interface_summary = json.loads(
        args.interface_summary.read_text(encoding="utf-8")
    )
    passing_ids = list(interface_summary["passing_ids"])
    if summary.get("model") != "model_2_ptm" or summary.get("template_mode") != "ct":
        raise ValueError("Positive summary is not the R10 model-2 pTM ct screen")
    if not passing_ids:
        raise ValueError("No R10 positive-screen passer cleared the interface audit")
    if not set(passing_ids).issubset(set(summary["passing_ids"])):
        raise ValueError("Interface audit contains a non-positive candidate")

    targets = {
        "USP4": args.usp4_target,
        "USP11": args.usp11_target,
    }
    target_lines = {name: atom_lines(path, "B") for name, path in targets.items()}
    records = []
    for candidate_id in passing_ids:
        source = args.panel_dir / f"{candidate_id}.pdb"
        if not source.is_file():
            raise FileNotFoundError(source)
        binder = atom_lines(source, "A")
        sequence = sequence_from_binder(binder)
        if "C" in sequence:
            raise ValueError(f"{candidate_id} contains cysteine")
        outputs = {}
        for target_name in targets:
            output = args.output_dir / target_name / f"{candidate_id}.pdb"
            write_complex(binder, target_lines[target_name], output)
            outputs[target_name] = str(output)
        records.append(
            {
                "id": candidate_id,
                "sequence": sequence,
                "length": len(sequence),
                "source_pdb": str(source),
                "source_sha256": sha256(source),
                "outputs": outputs,
            }
        )

    report = {
        "phase": "R10 same-pose USP4/USP11 selectivity input preparation",
        "positive_summary": str(args.positive_summary),
        "interface_summary": str(args.interface_summary),
        "positive_count": len(records),
        "targets": {
            name: {
                "pdb": str(path),
                "sha256": sha256(path),
                "chain": "B",
            }
            for name, path in targets.items()
        },
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"positive_count": len(records)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
