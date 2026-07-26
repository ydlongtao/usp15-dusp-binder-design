#!/usr/bin/env python3
"""Extract and align USP4/USP11 DUSP targets into the USP15 reference frame."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.PDB import PDBParser


AA3_TO_1 = {
    name.upper(): letter.upper()
    for name, letter in protein_letters_3to1_extended.items()
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--usp15-reference", required=True, type=Path)
    parser.add_argument("--usp4-source", required=True, type=Path)
    parser.add_argument("--usp11-source", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def chain_residues(path: Path, chain_id: str) -> list[dict]:
    structure = PDBParser(QUIET=True).get_structure(path.stem, path)
    chain = structure[0][chain_id]
    residues = []
    for residue in chain:
        if residue.id[0] != " " or "CA" not in residue:
            continue
        aa = AA3_TO_1.get(residue.resname.upper())
        if aa is None:
            continue
        residues.append(
            {
                "key": (int(residue.id[1]), str(residue.id[2]).strip()),
                "aa": aa,
                "ca": np.asarray(residue["CA"].coord, dtype=float),
            }
        )
    if not residues:
        raise ValueError(f"No standard C-alpha residues in {path} chain {chain_id}")
    return residues


def available_chains(path: Path) -> list[str]:
    structure = PDBParser(QUIET=True).get_structure(path.stem, path)
    return [
        chain.id
        for chain in structure[0]
        if any(residue.id[0] == " " and "CA" in residue for residue in chain)
    ]


def kabsch(
    mobile: np.ndarray, reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u_matrix @ vt_matrix))
    rotation = u_matrix @ correction @ vt_matrix
    translation = reference_center - mobile_center @ rotation
    aligned = mobile @ rotation + translation
    rmsd = float(
        np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1)))
    )
    return rotation, translation, rmsd


def align_chain(reference: list[dict], homolog: list[dict]) -> dict:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(
        "".join(item["aa"] for item in reference),
        "".join(item["aa"] for item in homolog),
    )[0]
    pairs: list[tuple[int, int]] = []
    for (ref_start, ref_end), (hom_start, hom_end) in zip(*alignment.aligned):
        pairs.extend(zip(range(ref_start, ref_end), range(hom_start, hom_end)))
    if len(pairs) < 60:
        raise ValueError(f"Insufficient DUSP alignment: only {len(pairs)} pairs")

    reference_ca = np.asarray([reference[i]["ca"] for i, _ in pairs])
    homolog_ca = np.asarray([homolog[j]["ca"] for _, j in pairs])
    rotation, translation, rmsd = kabsch(homolog_ca, reference_ca)
    identity = sum(
        reference[i]["aa"] == homolog[j]["aa"] for i, j in pairs
    ) / len(pairs)

    ref_min = min(i for i, _ in pairs)
    ref_max = max(i for i, _ in pairs)
    hom_min = min(j for _, j in pairs)
    hom_max = max(j for _, j in pairs)
    start_index = max(0, hom_min - ref_min)
    end_index = min(
        len(homolog) - 1,
        hom_max + (len(reference) - 1 - ref_max),
    )
    return {
        "score": float(alignment.score),
        "pairs": pairs,
        "pair_count": len(pairs),
        "identity": float(identity),
        "rmsd": rmsd,
        "rotation": rotation,
        "translation": translation,
        "start_index": start_index,
        "end_index": end_index,
    }


def transform_target(
    source: Path,
    source_chain: str,
    selected_keys: set[tuple[int, str]],
    rotation: np.ndarray,
    translation: np.ndarray,
    output: Path,
) -> int:
    lines = []
    residues = set()
    atom_serial = 1
    for line in source.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[21] != source_chain:
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        key = (int(line[22:26]), line[26].strip())
        if key not in selected_keys:
            continue
        xyz = np.asarray(
            [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]
        )
        transformed = xyz @ rotation + translation
        updated = (
            f"{line[:6]}{atom_serial:5d}{line[11:16]} {line[17:21]}B"
            f"{line[22:30]}{transformed[0]:8.3f}{transformed[1]:8.3f}"
            f"{transformed[2]:8.3f}{line[54:]}"
        )
        lines.append(updated)
        residues.add(key)
        atom_serial += 1
    if not lines:
        raise ValueError(f"No selected atoms written from {source} chain {source_chain}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\nTER\nEND\n", encoding="utf-8")
    return len(residues)


def prepare_one(
    name: str,
    source: Path,
    reference: list[dict],
    output: Path,
) -> dict:
    candidates = []
    for chain_id in available_chains(source):
        homolog = chain_residues(source, chain_id)
        try:
            result = align_chain(reference, homolog)
        except ValueError:
            continue
        result["chain_id"] = chain_id
        result["homolog"] = homolog
        candidates.append(result)
    if not candidates:
        raise ValueError(f"No alignable chain found in {source}")
    best = sorted(
        candidates,
        key=lambda item: (-item["score"], item["rmsd"], item["chain_id"]),
    )[0]
    homolog = best["homolog"]
    selected = homolog[best["start_index"] : best["end_index"] + 1]
    selected_keys = {item["key"] for item in selected}
    residue_count = transform_target(
        source,
        best["chain_id"],
        selected_keys,
        best["rotation"],
        best["translation"],
        output,
    )
    return {
        "target": name,
        "source_pdb": str(source),
        "source_sha256": sha256(source),
        "source_chain": best["chain_id"],
        "source_observed_residues": len(homolog),
        "alignment_score": best["score"],
        "alignment_pair_count": best["pair_count"],
        "alignment_identity": best["identity"],
        "alignment_ca_rmsd_A": best["rmsd"],
        "selected_source_residue_start": selected[0]["key"][0],
        "selected_source_residue_end": selected[-1]["key"][0],
        "selected_residue_count": residue_count,
        "output_pdb": str(output),
        "output_chain": "B",
    }


def main() -> None:
    args = parse_args()
    reference = chain_residues(args.usp15_reference, "B")
    if len(reference) != 129:
        raise ValueError(
            f"Expected 129 USP15 target residues, found {len(reference)}"
        )
    records = [
        prepare_one(
            "USP4",
            args.usp4_source,
            reference,
            args.output_dir / "USP4_5CTR_DUSP_aligned_chainB.pdb",
        ),
        prepare_one(
            "USP11",
            args.usp11_source,
            reference,
            args.output_dir / "USP11_4MEL_DUSP_aligned_chainB.pdb",
        ),
    ]
    report = {
        "phase": "R10 USP4/USP11 DUSP structural preparation",
        "reference_pdb": str(args.usp15_reference),
        "reference_sha256": sha256(args.usp15_reference),
        "reference_chain": "B",
        "reference_residues": len(reference),
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    main()
