#!/usr/bin/env python3
"""Transplant four independent 6DJ9 UbV poses onto the complete USP15 target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


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
POSES = {
    "AK": ("A", "K"),
    "BL": ("B", "L"),
    "CJ": ("C", "J"),
    "DH": ("D", "H"),
}
HOTSPOTS = {50, 52, 53, 55, 57, 61}
INTERFACE_CENTER_RESIDUES = (50, 52, 53, 55, 57)
FIXED_INTERFACE = {4, 6, 7, 8, 9, 44, 46, 48, 49, 50, 51, 72, 73, 74, 75}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--complex", required=True, type=Path)
    parser.add_argument("--target-reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def atom_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.startswith("ATOM") and line[16] in {" ", "A"}
    ]


def coordinates(line: str) -> np.ndarray:
    return np.array(
        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
        dtype=float,
    )


def ca_coordinates(lines: list[str], chain: str) -> dict[int, np.ndarray]:
    return {
        int(line[22:26]): coordinates(line)
        for line in lines
        if line[21] == chain and line[12:16].strip() == "CA"
    }


def fit_transform(
    mobile: np.ndarray, fixed: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    mobile_center = mobile.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    u_matrix, _, vt_matrix = np.linalg.svd(
        (mobile - mobile_center).T @ (fixed - fixed_center)
    )
    row_rotation = u_matrix @ vt_matrix
    if np.linalg.det(row_rotation) < 0:
        vt_matrix[-1, :] *= -1
        row_rotation = u_matrix @ vt_matrix
    rotation = row_rotation.T
    translation = fixed_center - rotation @ mobile_center
    fitted = (rotation @ mobile.T).T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - fixed) ** 2, axis=1))))
    return rotation, translation, rmsd


def rewrite_atom(
    line: str,
    output_chain: str,
    rotation: np.ndarray | None = None,
    translation: np.ndarray | None = None,
    shift: np.ndarray | None = None,
) -> str:
    chars = list(line.ljust(80))
    chars[21] = output_chain
    if rotation is not None and translation is not None:
        transformed = rotation @ coordinates(line) + translation
        if shift is not None:
            transformed += shift
        chars[30:54] = (
            f"{transformed[0]:8.3f}{transformed[1]:8.3f}{transformed[2]:8.3f}"
        )
    return "".join(chars)


def compact_ranges(residues: list[int]) -> str:
    ranges: list[tuple[int, int]] = []
    start = previous = residues[0]
    for residue in residues[1:]:
        if residue == previous + 1:
            previous = residue
        else:
            ranges.append((start, previous))
            start = previous = residue
    ranges.append((start, previous))
    return "/".join(
        f"A{start}" if start == end else f"A{start}-{end}"
        for start, end in ranges
    )


def residue_sequence(lines: list[str], chain: str) -> str:
    residues: dict[int, str] = {}
    for line in lines:
        if line[21] != chain:
            continue
        residue = int(line[22:26])
        residues.setdefault(residue, AA3_TO_1[line[17:20].strip()])
    return "".join(residues[index] for index in sorted(residues))


def interface_diagnostics(
    binder_lines: list[str], target_lines: list[str]
) -> dict[str, object]:
    binder = [(int(line[22:26]), coordinates(line)) for line in binder_lines]
    target = [(int(line[22:26]), coordinates(line)) for line in target_lines]
    binder_coordinates = np.stack([coordinate for _, coordinate in binder])
    target_coordinates = np.stack([coordinate for _, coordinate in target])
    distances = np.linalg.norm(
        binder_coordinates[:, np.newaxis, :]
        - target_coordinates[np.newaxis, :, :],
        axis=2,
    )
    diagnostics: dict[str, object] = {}
    for cutoff in (5.0, 6.0):
        binder_indices, target_indices = np.where(distances < cutoff)
        diagnostics[f"binder_contact_residues_{cutoff:g}A"] = sorted(
            {binder[index][0] for index in binder_indices}
        )
        diagnostics[f"target_contact_residues_{cutoff:g}A"] = sorted(
            {target[index][0] for index in target_indices}
        )
    target_contacts = set(diagnostics["target_contact_residues_5A"])
    diagnostics.update(
        {
            "minimum_interchain_atom_distance_A": float(distances.min()),
            "interchain_atom_pairs_below_2A": int((distances < 2.0).sum()),
            "hotspots_contacted_5A": sorted(HOTSPOTS & target_contacts),
        }
    )
    return diagnostics


def find_declash_translation(
    binder_lines: list[str],
    target_lines: list[str],
    outward: np.ndarray,
) -> tuple[np.ndarray, dict[str, object], dict[str, object]]:
    """Find the smallest deterministic rigid translation satisfying hard gates."""
    binder_coordinates = np.stack([coordinates(line) for line in binder_lines])
    target_coordinates = np.stack([coordinates(line) for line in target_lines])
    target_residues = np.array([int(line[22:26]) for line in target_lines])
    successful: list[
        tuple[float, np.ndarray, dict[str, object], dict[str, object]]
    ] = []
    searches: list[dict[str, object]] = []
    seed_distances = np.arange(0.0, 2.5001, 0.25)
    for seed_distance in seed_distances:
        shift = float(seed_distance) * outward
        minimum_distance = 0.0
        hotspot_contacts: list[int] = []
        clash_count = 0
        for iteration in range(300):
            displaced = binder_coordinates + shift
            deltas = (
                displaced[:, np.newaxis, :]
                - target_coordinates[np.newaxis, :, :]
            )
            distances = np.linalg.norm(deltas, axis=2)
            clash_count = int((distances < 2.0).sum())
            minimum_distance = float(distances.min())
            _, target_indices = np.where(distances < 5.0)
            contacted_target = set(target_residues[target_indices].tolist())
            hotspot_contacts = sorted(HOTSPOTS & contacted_target)
            if clash_count == 0 and len(hotspot_contacts) >= 4:
                shifted_lines = []
                for line in binder_lines:
                    chars = list(line.ljust(80))
                    shifted = coordinates(line) + shift
                    chars[30:54] = (
                        f"{shifted[0]:8.3f}{shifted[1]:8.3f}{shifted[2]:8.3f}"
                    )
                    shifted_lines.append("".join(chars))
                diagnostics = interface_diagnostics(shifted_lines, target_lines)
                search_record = {
                    "method": "iterative_interchain_repulsion",
                    "seed_outward_shift_A": float(seed_distance),
                    "iterations": iteration,
                    "translation_vector_A": shift.tolist(),
                    "translation_norm_A": float(np.linalg.norm(shift)),
                }
                successful.append(
                    (
                        float(np.linalg.norm(shift)),
                        shift.copy(),
                        diagnostics,
                        search_record,
                    )
                )
                break

            repulsion_mask = distances < 2.1
            binder_indices, target_indices = np.where(repulsion_mask)
            if len(binder_indices) == 0:
                break
            pair_deltas = deltas[binder_indices, target_indices]
            pair_distances = distances[binder_indices, target_indices]
            pair_directions = pair_deltas / np.maximum(
                pair_distances[:, np.newaxis], 1e-6
            )
            weights = 2.1 - pair_distances
            repulsion = (weights[:, np.newaxis] * pair_directions).sum(axis=0)
            if np.linalg.norm(repulsion) < 1e-8:
                closest = int(np.argmin(pair_distances))
                repulsion = pair_directions[closest]
            direction = repulsion / np.linalg.norm(repulsion)
            step = max(0.01, min(0.08, 2.01 - minimum_distance))
            shift = shift + step * direction
            if np.linalg.norm(shift) > 4.0:
                break
        searches.append(
            {
                "seed_outward_shift_A": float(seed_distance),
                "final_translation_norm_A": float(np.linalg.norm(shift)),
                "final_minimum_interchain_atom_distance_A": minimum_distance,
                "final_interchain_atom_pairs_below_2A": clash_count,
                "final_hotspots_contacted_5A": hotspot_contacts,
            }
        )
    if not successful:
        raise ValueError(
            "No rigid translation satisfies zero clashes and four hotspots: "
            + json.dumps(searches)
        )
    _, translation, diagnostics, search_record = min(
        successful, key=lambda item: item[0]
    )
    return translation, diagnostics, search_record


def renumber_and_write(
    path: Path,
    binder_lines: list[str],
    target_lines: list[str],
    inpaint_seq: str,
) -> None:
    output = [
        'REMARK   1 Input contig: "A1-76/0 B6-134/0"',
        'REMARK   1 Standardized contig: "A1-76/0 B6-134/0"',
        f'REMARK   1 Inpaint seq: "{inpaint_seq}"',
        'REMARK   1 Chains: "A B"',
        'REMARK   1 Input hotspots: "B50,B52,B53,B55,B57,B61"',
        'REMARK   1 Standardized hotspots: "B50,B52,B53,B55,B57,B61"',
    ]
    serial = 0
    for line in binder_lines + target_lines:
        serial += 1
        chars = list(line.ljust(80))
        chars[6:11] = f"{serial:5d}"
        output.append("".join(chars))
    output.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    complex_lines = atom_lines(args.complex)
    reference_lines = atom_lines(args.target_reference)
    reference_ca = ca_coordinates(reference_lines, "A")
    target_lines = [
        rewrite_atom(line, "B")
        for line in reference_lines
        if line[21] == "A" and 6 <= int(line[22:26]) <= 134
    ]
    if sorted(reference_ca) != list(range(6, 135)):
        raise ValueError("Expected complete target-reference CA residues A6-134")

    design_positions = sorted(set(range(1, 77)) - FIXED_INTERFACE)
    inpaint_seq = compact_ranges(design_positions)
    records: list[dict[str, object]] = []
    output_paths: list[str] = []
    for pose_name, (target_chain, binder_chain) in POSES.items():
        source_target_ca = ca_coordinates(complex_lines, target_chain)
        fitting_residues = list(range(25, 71))
        if not all(
            residue in source_target_ca and residue in reference_ca
            for residue in fitting_residues
        ):
            raise ValueError(f"Pose {pose_name} lacks target CA residues 25-70")
        rotation, translation, alignment_rmsd = fit_transform(
            np.stack([source_target_ca[index] for index in fitting_residues]),
            np.stack([reference_ca[index] for index in fitting_residues]),
        )
        source_binder = [
            line
            for line in complex_lines
            if line[21] == binder_chain and 1 <= int(line[22:26]) <= 76
        ]
        binder_residues = sorted({int(line[22:26]) for line in source_binder})
        if binder_residues != list(range(1, 77)):
            raise ValueError(f"Pose {pose_name} lacks complete binder residues 1-76")
        transformed_zero = [
            rewrite_atom(line, "A", rotation, translation, np.zeros(3))
            for line in source_binder
        ]
        binder_center = np.stack(
            [coordinates(line) for line in transformed_zero]
        ).mean(axis=0)
        interface_center = np.stack(
            [reference_ca[index] for index in INTERFACE_CENTER_RESIDUES]
        ).mean(axis=0)
        outward = binder_center - interface_center
        outward /= np.linalg.norm(outward)

        selected_translation, selected_diagnostics, declash_search = (
            find_declash_translation(transformed_zero, target_lines, outward)
        )
        selected_binder = []
        for line in transformed_zero:
            chars = list(line.ljust(80))
            shifted = coordinates(line) + selected_translation
            chars[30:54] = (
                f"{shifted[0]:8.3f}{shifted[1]:8.3f}{shifted[2]:8.3f}"
            )
            selected_binder.append("".join(chars))

        sequence = residue_sequence(selected_binder, "A")
        fixed_sequence = "".join(sequence[index - 1] for index in sorted(FIXED_INTERFACE))
        if len(sequence) != 76:
            raise ValueError(f"Pose {pose_name} binder sequence is not 76 aa")
        if "C" in fixed_sequence:
            raise ValueError(f"Pose {pose_name} fixes a Cys at the interface")
        output_path = args.output_dir / f"pose_{pose_name}_fixed_interface.pdb"
        renumber_and_write(
            output_path,
            selected_binder,
            target_lines,
            inpaint_seq,
        )
        output_paths.append(output_path.name)
        records.append(
            {
                "pose": pose_name,
                "source_target_chain": target_chain,
                "source_binder_chain": binder_chain,
                "target_alignment": "6DJ9 residues 25-70 to 3T9L residues 25-70",
                "target_ca_alignment_rmsd_A": alignment_rmsd,
                "declash_translation": declash_search,
                "sequence": sequence,
                "length": len(sequence),
                "binder_cys_count_before_design": sequence.count("C"),
                "fixed_interface_positions": sorted(FIXED_INTERFACE),
                "fixed_interface_sequence": fixed_sequence,
                "interface_diagnostics": selected_diagnostics,
                "output_pdb": output_path.name,
            }
        )

    report = {
        "phase": "R4 6DJ9 crystal-pose ensemble",
        "source_complex": "PDB 6DJ9 asymmetric unit",
        "target_reference": args.target_reference.name,
        "selected_backbones": len(records),
        "pose_pairs": {
            pose: {"target_chain": chains[0], "binder_chain": chains[1]}
            for pose, chains in POSES.items()
        },
        "fixed_crystallographic_interface_positions": sorted(FIXED_INTERFACE),
        "redesigned_core_positions": design_positions,
        "fixed_position_count": len(FIXED_INTERFACE),
        "redesigned_position_count": len(design_positions),
        "inpaint_seq": inpaint_seq,
        "sequence_model": "protein_mpnn",
        "temperature": 0.1,
        "sequences_per_backbone": 3,
        "omit_amino_acids": "C",
        "output_pdbs": output_paths,
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
