#!/usr/bin/env python3
"""Place a stable 1UBQ scaffold at the 6DJ9 USP15 interface."""

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
HOTSPOTS = {50, 52, 53, 55, 57, 61}
UBIQUITIN_CORE_MAPPING = (
    [(index, index) for index in range(1, 8)]
    + [(index, index + 2) for index in range(8, 75)]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ubiquitin", required=True, type=Path)
    parser.add_argument("--complex", required=True, type=Path)
    parser.add_argument("--target-reference", required=True, type=Path)
    parser.add_argument("--wt-output", required=True, type=Path)
    parser.add_argument("--design-output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def atom_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.startswith("ATOM") and line[16] in {" ", "A"}
    ]


def ca_coordinates(lines: list[str], chain: str) -> dict[int, np.ndarray]:
    return {
        int(line[22:26]): np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=float,
        )
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


def transform_atom(
    line: str,
    output_chain: str,
    rotation_1: np.ndarray,
    translation_1: np.ndarray,
    rotation_2: np.ndarray,
    translation_2: np.ndarray,
    shift: np.ndarray,
) -> str:
    coordinate = np.array(
        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
        dtype=float,
    )
    coordinate = (
        rotation_2 @ (rotation_1 @ coordinate + translation_1)
        + translation_2
        + shift
    )
    chars = list(line.ljust(80))
    chars[21] = output_chain
    chars[30:54] = (
        f"{coordinate[0]:8.3f}{coordinate[1]:8.3f}{coordinate[2]:8.3f}"
    )
    return "".join(chars)


def rewrite_target(line: str) -> str:
    chars = list(line.ljust(80))
    chars[21] = "B"
    return "".join(chars)


def residue_sequence(lines: list[str], chain: str) -> str:
    residues = {}
    for line in lines:
        if line.startswith("ATOM") and line[21] == chain:
            residues.setdefault(int(line[22:26]), AA3_TO_1[line[17:20].strip()])
    return "".join(residues[index] for index in sorted(residues))


def compact_ranges(residues: list[int]) -> str:
    ranges = []
    start = previous = residues[0]
    for residue in residues[1:]:
        if residue == previous + 1:
            previous = residue
            continue
        ranges.append((start, previous))
        start = previous = residue
    ranges.append((start, previous))
    return "/".join(
        f"A{start}" if start == end else f"A{start}-{end}"
        for start, end in ranges
    )


def interface_diagnostics(
    binder_lines: list[str], target_lines: list[str]
) -> dict[str, object]:
    binder = [
        (
            int(line[22:26]),
            np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            ),
        )
        for line in binder_lines
        if line.startswith("ATOM")
    ]
    target = [
        (
            int(line[22:26]),
            np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            ),
        )
        for line in target_lines
        if line.startswith("ATOM")
    ]
    binder_coordinates = np.stack([coordinate for _, coordinate in binder])
    target_coordinates = np.stack([coordinate for _, coordinate in target])
    distances = np.sqrt(
        (
            (
                binder_coordinates[:, np.newaxis, :]
                - target_coordinates[np.newaxis, :, :]
            )
            ** 2
        ).sum(axis=2)
    )
    diagnostics = {}
    for cutoff in (5.0, 6.0):
        indices = np.where(distances < cutoff)
        diagnostics[f"binder_contact_residues_{cutoff:g}A"] = sorted(
            {binder[index][0] for index in indices[0]}
        )
        diagnostics[f"target_contact_residues_{cutoff:g}A"] = sorted(
            {target[index][0] for index in indices[1]}
        )
    target_contacts = diagnostics["target_contact_residues_5A"]
    diagnostics.update(
        {
            "minimum_interchain_atom_distance_A": float(distances.min()),
            "interchain_atom_pairs_below_2A": int((distances < 2.0).sum()),
            "hotspots_contacted_5A": sorted(HOTSPOTS & set(target_contacts)),
        }
    )
    return diagnostics


def main() -> int:
    args = parse_args()
    ubiquitin_lines = atom_lines(args.ubiquitin)
    complex_lines = atom_lines(args.complex)
    target_lines = atom_lines(args.target_reference)
    ubiquitin_ca = ca_coordinates(ubiquitin_lines, "A")
    ubv_ca = ca_coordinates(complex_lines, "K")
    source_target_ca = ca_coordinates(complex_lines, "A")
    reference_target_ca = ca_coordinates(target_lines, "A")

    # The UbV has a two-residue insertion after structurally equivalent position
    # 7. The well-conserved C-terminal core gives a stable pose transform.
    core_mapping = [
        pair for pair in UBIQUITIN_CORE_MAPPING if 20 <= pair[0] <= 74
    ]
    rotation_1, translation_1, scaffold_rmsd = fit_transform(
        np.stack([ubiquitin_ca[wt] for wt, _ in core_mapping]),
        np.stack([ubv_ca[ubv] for _, ubv in core_mapping]),
    )
    target_mapping = list(range(25, 71))
    rotation_2, translation_2, target_rmsd = fit_transform(
        np.stack([source_target_ca[index] for index in target_mapping]),
        np.stack([reference_target_ca[index] for index in target_mapping]),
    )

    preliminary_binder = []
    for line in ubiquitin_lines:
        if line[21] == "A" and 1 <= int(line[22:26]) <= 76:
            preliminary_binder.append(
                transform_atom(
                    line,
                    "A",
                    rotation_1,
                    translation_1,
                    rotation_2,
                    translation_2,
                    np.zeros(3),
                )
            )
    output_target = [
        rewrite_target(line)
        for line in target_lines
        if line[21] == "A" and 6 <= int(line[22:26]) <= 134
    ]

    binder_coordinates = np.stack(
        [
            np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
            for line in preliminary_binder
        ]
    )
    interface_center = np.stack(
        [reference_target_ca[index] for index in (50, 52, 53, 55, 57)]
    ).mean(axis=0)
    outward = binder_coordinates.mean(axis=0) - interface_center
    outward /= np.linalg.norm(outward)
    shift = 0.75 * outward

    output_binder = []
    for line in ubiquitin_lines:
        if line[21] == "A" and 1 <= int(line[22:26]) <= 76:
            output_binder.append(
                transform_atom(
                    line,
                    "A",
                    rotation_1,
                    translation_1,
                    rotation_2,
                    translation_2,
                    shift,
                )
            )
    diagnostics = interface_diagnostics(output_binder, output_target)
    design_residues = diagnostics["binder_contact_residues_6A"]
    if len(diagnostics["hotspots_contacted_5A"]) < 4:
        raise ValueError("Placed ubiquitin scaffold contacts fewer than four hotspots")
    if diagnostics["interchain_atom_pairs_below_2A"] != 0:
        raise ValueError("Placed ubiquitin scaffold retains interchain clashes")

    def write_output(path: Path, inpaint: str | None) -> None:
        header = [
            'REMARK   1 Input contig: "A1-76/0 B6-134/0"',
            'REMARK   1 Standardized contig: "A1-76/0 B6-134/0"',
            'REMARK   1 Chains: "A B"',
            'REMARK   1 Input hotspots: "B50,B52,B53,B55,B57,B61"',
            'REMARK   1 Standardized hotspots: "B50,B52,B53,B55,B57,B61"',
        ]
        if inpaint:
            header.append(f'REMARK   1 Inpaint seq: "{inpaint}"')
        output = header + output_binder + output_target
        serial = 0
        for index, line in enumerate(output):
            if line.startswith("ATOM"):
                serial += 1
                chars = list(line.ljust(80))
                chars[6:11] = f"{serial:5d}"
                output[index] = "".join(chars)
        output.append("END")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(output) + "\n", encoding="utf-8")

    inpaint_seq = compact_ranges(design_residues)
    write_output(args.wt_output, None)
    write_output(args.design_output, inpaint_seq)
    sequence = residue_sequence(output_binder, "A")
    report = {
        "source_scaffold": "PDB 1UBQ chain A residues 1-76",
        "source_interface": "PDB 6DJ9 chains K/A",
        "target_reference": args.target_reference.name,
        "scaffold_core_mapping": "1UBQ A20-74 to 6DJ9 K22-76",
        "scaffold_core_ca_rmsd_A": scaffold_rmsd,
        "target_interface_mapping": "6DJ9 A25-70 to target A25-70",
        "target_interface_ca_rmsd_A": target_rmsd,
        "outward_shift_A": 0.75,
        "sequence": sequence,
        "length": len(sequence),
        "binder_cys_count": sequence.count("C"),
        "designed_interface_residues": design_residues,
        "inpaint_seq": inpaint_seq,
        "interface_diagnostics": diagnostics,
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
