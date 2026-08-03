#!/usr/bin/env python3
"""Extract 6DJ9 A/K and create WT plus Cys11-free UbV controls."""

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
VARIANTS = {
    "ubv15d_wt_control": "CYS",
    "ubv15d_c11a": "ALA",
    "ubv15d_c11s": "SER",
    "ubv15d_c11v": "VAL",
}
HOTSPOTS = {50, 52, 53, 55, 57, 61}
ALIGNMENT_RESIDUES = set(range(25, 71))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--target-reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--native-diagnostic-output", type=Path)
    return parser.parse_args()


def atom_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.startswith("ATOM") and line[16] in {" ", "A"}
    ]


def residue_sequence(lines: list[str], chain: str) -> tuple[str, list[int]]:
    seen: set[tuple[int, str]] = set()
    sequence: list[str] = []
    residues: list[int] = []
    for line in lines:
        if line[21] != chain:
            continue
        residue = int(line[22:26])
        insertion = line[26]
        key = (residue, insertion)
        if key in seen:
            continue
        seen.add(key)
        if insertion != " ":
            raise ValueError(f"Insertion code in chain {chain} residue {residue}")
        name = line[17:20].strip()
        if name not in AA3_TO_1:
            raise ValueError(f"Unsupported residue {name}")
        residues.append(residue)
        sequence.append(AA3_TO_1[name])
    return "".join(sequence), residues


def rewrite_atom(
    line: str,
    output_chain: str,
    output_residue: int,
    residue_name: str | None = None,
    atom_name: str | None = None,
    coordinates: np.ndarray | None = None,
) -> str:
    chars = list(line.ljust(80))
    chars[21] = output_chain
    chars[22:26] = f"{output_residue:4d}"
    chars[26] = " "
    if residue_name is not None:
        chars[17:20] = f"{residue_name:>3}"
    if atom_name is not None:
        chars[12:16] = f"{atom_name:>4}"
        chars[76:78] = f"{atom_name[0]:>2}"
    if coordinates is not None:
        chars[30:54] = f"{coordinates[0]:8.3f}{coordinates[1]:8.3f}{coordinates[2]:8.3f}"
    return "".join(chars)


def ca_coordinates(lines: list[str], chain: str) -> dict[int, np.ndarray]:
    coordinates = {}
    for line in lines:
        if line[21] == chain and line[12:16].strip() == "CA":
            coordinates[int(line[22:26])] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
    return coordinates


def alignment_transform(
    source_lines: list[str], reference_lines: list[str]
) -> tuple[np.ndarray, np.ndarray, float, int]:
    source_ca = ca_coordinates(source_lines, "A")
    reference_ca = ca_coordinates(reference_lines, "A")
    common = sorted(set(source_ca) & set(reference_ca) & ALIGNMENT_RESIDUES)
    if len(common) < 3:
        raise ValueError("Too few common target C-alpha atoms for alignment")
    mobile = np.stack([source_ca[index] for index in common])
    fixed = np.stack([reference_ca[index] for index in common])
    mobile_center = mobile.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (fixed - fixed_center)
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = vt_matrix.T @ u_matrix.T
    if np.linalg.det(rotation) < 0:
        vt_matrix[-1, :] *= -1
        rotation = vt_matrix.T @ u_matrix.T
    translation = fixed_center - rotation @ mobile_center
    fitted = (rotation @ mobile.T).T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - fixed) ** 2, axis=1))))
    return rotation, translation, rmsd, len(common)


def interface_diagnostics(lines: list[str]) -> dict[str, object]:
    binder_atoms = []
    target_atoms = []
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        coordinate = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=float,
        )
        residue = int(line[22:26])
        if line[21] == "A":
            binder_atoms.append((residue, coordinate))
        elif line[21] == "B":
            target_atoms.append((residue, coordinate))
    binder_coordinates = np.stack([coordinate for _, coordinate in binder_atoms])
    target_coordinates = np.stack([coordinate for _, coordinate in target_atoms])
    distances = np.sqrt(
        (
            (
                binder_coordinates[:, np.newaxis, :]
                - target_coordinates[np.newaxis, :, :]
            )
            ** 2
        ).sum(axis=2)
    )
    contact_indices = np.where(distances < 5.0)
    target_contacts = sorted({target_atoms[index][0] for index in contact_indices[1]})
    binder_contacts = sorted({binder_atoms[index][0] for index in contact_indices[0]})
    return {
        "minimum_interchain_atom_distance_A": float(distances.min()),
        "interchain_atom_pairs_below_2A": int((distances < 2.0).sum()),
        "binder_contact_residues_5A": binder_contacts,
        "target_contact_residues_5A": target_contacts,
        "hotspots_contacted_5A": sorted(HOTSPOTS & set(target_contacts)),
    }


def build_variant(
    source_lines: list[str],
    reference_lines: list[str],
    variant_name: str,
    position_11_name: str,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> tuple[list[str], str]:
    output: list[str] = [
        'REMARK   1 Input contig: "A6-134/0 76-76"',
        'REMARK   1 Standardized contig: "76-76 A6-134/0"',
        'REMARK   1 Chains: "A B"',
        'REMARK   1 Input hotspots: "A50,A52,A53,A55,A57,A61"',
        'REMARK   1 Standardized hotspots: "B50,B52,B53,B55,B57,B61"',
    ]
    binder_names: dict[int, str] = {}
    for line in source_lines:
        source_chain = line[21]
        source_residue = int(line[22:26])
        atom = line[12:16].strip()
        if source_chain == "K" and 1 <= source_residue <= 76:
            residue_name = line[17:20].strip()
            atom_name = None
            if source_residue == 11:
                residue_name = position_11_name
                if position_11_name == "ALA" and atom == "SG":
                    continue
                if position_11_name == "SER" and atom == "SG":
                    atom_name = "OG"
                if position_11_name == "VAL" and atom == "SG":
                    continue
            binder_names[source_residue] = residue_name
            source_coordinates = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
            output.append(
                rewrite_atom(
                    line,
                    output_chain="A",
                    output_residue=source_residue,
                    residue_name=residue_name,
                    atom_name=atom_name,
                    coordinates=rotation @ source_coordinates + translation,
                )
            )
    for line in reference_lines:
        source_chain = line[21]
        source_residue = int(line[22:26])
        if source_chain == "A" and 6 <= source_residue <= 134:
            output.append(
                rewrite_atom(
                    line,
                    output_chain="B",
                    output_residue=source_residue,
                )
            )

    if sorted(binder_names) != list(range(1, 77)):
        raise ValueError(f"{variant_name}: incomplete binder residues")
    sequence = "".join(AA3_TO_1[binder_names[index]] for index in range(1, 77))
    atom_serial = 0
    for index, line in enumerate(output):
        if line.startswith("ATOM"):
            atom_serial += 1
            chars = list(line.ljust(80))
            chars[6:11] = f"{atom_serial:5d}"
            output[index] = "".join(chars)
    output.append("END")
    return output, sequence


def build_native_diagnostic(source_lines: list[str]) -> list[str]:
    output = [
        'REMARK   1 Diagnostic-only native 6DJ9 A/K complex',
        'REMARK   1 Binder chain A is source K1-76',
        'REMARK   1 Target chain B is resolved source A3-134 (A76 unresolved)',
    ]
    for line in source_lines:
        source_chain = line[21]
        source_residue = int(line[22:26])
        if source_chain == "K" and 1 <= source_residue <= 76:
            output.append(rewrite_atom(line, "A", source_residue))
    for line in source_lines:
        source_chain = line[21]
        source_residue = int(line[22:26])
        if source_chain == "A" and 3 <= source_residue <= 134:
            output.append(rewrite_atom(line, "B", source_residue))
    atom_serial = 0
    for index, line in enumerate(output):
        if line.startswith("ATOM"):
            atom_serial += 1
            chars = list(line.ljust(80))
            chars[6:11] = f"{atom_serial:5d}"
            output[index] = "".join(chars)
    output.append("END")
    return output


def main() -> int:
    args = parse_args()
    source_lines = atom_lines(args.source)
    reference_lines = atom_lines(args.target_reference)
    _, binder_residues = residue_sequence(source_lines, "K")
    _, source_target_residues = residue_sequence(source_lines, "A")
    _, target_residues = residue_sequence(reference_lines, "A")
    if not set(range(1, 77)).issubset(binder_residues):
        raise ValueError("6DJ9 chain K does not contain residues 1–76")
    if not set(range(6, 135)).issubset(target_residues):
        raise ValueError("Target reference does not contain residues 6–134")
    if not HOTSPOTS.issubset(target_residues):
        raise ValueError("Target hotspots are incomplete")
    if not HOTSPOTS.issubset(source_target_residues):
        raise ValueError("6DJ9 target hotspots are incomplete")
    rotation, translation, alignment_rmsd, aligned_ca_count = alignment_transform(
        source_lines, reference_lines
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    fasta_lines: list[str] = []
    for variant_name, position_11_name in VARIANTS.items():
        output, sequence = build_variant(
            source_lines,
            reference_lines,
            variant_name,
            position_11_name,
            rotation,
            translation,
        )
        output_path = args.output_dir / f"{variant_name}.pdb"
        output_path.write_text("\n".join(output) + "\n", encoding="utf-8")
        fasta_lines.extend([f">{variant_name}", sequence])
        records.append(
            {
                "id": variant_name,
                "length": len(sequence),
                "sequence": sequence,
                "position_11": position_11_name,
                "binder_cys_count": sequence.count("C"),
                "candidate_eligible": position_11_name != "CYS",
                "interface_diagnostics": interface_diagnostics(output),
            }
        )

    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    args.fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    native_diagnostic = None
    if args.native_diagnostic_output is not None:
        native_output = build_native_diagnostic(source_lines)
        args.native_diagnostic_output.parent.mkdir(parents=True, exist_ok=True)
        args.native_diagnostic_output.write_text(
            "\n".join(native_output) + "\n", encoding="utf-8"
        )
        native_diagnostic = {
            "id": "ubv15d_native_wt_diagnostic",
            "path": args.native_diagnostic_output.name,
            "candidate_eligible": False,
            "target_residues_resolved": "3-75,77-134",
            "purpose": "diagnose target-conformation transfer; not a candidate",
            "interface_diagnostics": interface_diagnostics(native_output),
        }
    report = {
        "source": "PDB 6DJ9",
        "target_reference": args.target_reference.name,
        "source_target_chain": "A",
        "source_binder_chain": "K",
        "target_residues": "6-134",
        "binder_residues": "1-76",
        "standardized_binder_chain": "A",
        "standardized_target_chain": "B",
        "target_alignment_residues": "25-70",
        "target_alignment_ca_count": aligned_ca_count,
        "target_alignment_ca_rmsd_A": alignment_rmsd,
        "variants": records,
        "native_diagnostic": native_diagnostic,
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
