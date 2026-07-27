#!/usr/bin/env python3
"""Audit all prepared R10 systems before any production MD starts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import openmm
from openmm import CustomExternalForce, MonteCarloBarostat, XmlSerializer, unit
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, PDBFile
import parmed as pmd


EXPECTED_FILES = {
    "fixed_heavy.pdb",
    "protein_protonated.pdb",
    "solvated.pdb",
    "complex_tleap.prmtop",
    "complex.prmtop",
    "complex.inpcrd",
    "solvated.prmtop",
    "solvated.inpcrd",
    "binder.prmtop",
    "target.prmtop",
    "system_nvt.xml",
    "system_npt.xml",
    "system_production.xml",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def version_tuple(value: str):
    parts = re.findall(r"\d+", value)
    return tuple(int(part) for part in parts[:3])


def vector_angle_degrees(a, b):
    av = a.value_in_unit(unit.nanometer)
    bv = b.value_in_unit(unit.nanometer)
    dot = sum(float(x) * float(y) for x, y in zip(av, bv))
    anorm = math.sqrt(sum(float(x) ** 2 for x in av))
    bnorm = math.sqrt(sum(float(x) ** 2 for x in bv))
    return math.degrees(math.acos(dot / (anorm * bnorm)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-root", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--ranks",
        type=int,
        nargs="+",
        default=list(range(1, 11)),
    )
    args = parser.parse_args()
    rows = []
    for rank in args.ranks:
        path = args.prepared_root / f"rank{rank:02d}"
        metadata = json.loads((path / "preparation.json").read_text())
        if set(metadata["files"]) != EXPECTED_FILES:
            raise RuntimeError(f"rank {rank:02d} file manifest differs")
        for filename, expected in metadata["files"].items():
            artifact = path / filename
            if artifact.stat().st_size != expected["bytes"]:
                raise RuntimeError(
                    f"rank {rank:02d} byte count differs for {filename}"
                )
            if sha256(artifact) != expected["sha256"]:
                raise RuntimeError(
                    f"rank {rank:02d} SHA-256 differs for {filename}"
                )
        input_pdb = args.input_root / f"USP15_R10_rank{rank:02d}.pdb"
        if sha256(input_pdb) != metadata["input_sha256"]:
            raise RuntimeError(f"rank {rank:02d} input SHA-256 differs")

        system = XmlSerializer.deserialize(
            (path / "system_production.xml").read_text()
        )
        restraints = [
            force.getName()
            for force in system.getForces()
            if isinstance(force, CustomExternalForce)
        ]
        if restraints or metadata["production_restraints"] is not False:
            raise RuntimeError(f"rank {rank:02d} has production restraints")
        if version_tuple(metadata["openmm_version"]) < (8, 5):
            raise RuntimeError(f"rank {rank:02d} uses old OpenMM")
        if metadata["tleap_force_fields"] != [
            "leaprc.protein.ff19SB",
            "leaprc.water.opc",
        ]:
            raise RuntimeError(f"rank {rank:02d} force fields differ")
        if metadata["binder_residues"] != 76:
            raise RuntimeError(f"rank {rank:02d} binder length differs")
        if metadata["target_residues"] != 129:
            raise RuntimeError(f"rank {rank:02d} target length differs")
        if metadata["water_model"] != "OPC":
            raise RuntimeError(f"rank {rank:02d} water model differs")
        if metadata["box_shape"] != "truncated_octahedron":
            raise RuntimeError(f"rank {rank:02d} box shape differs")

        dry = pmd.load_file(
            str(path / "complex.prmtop"), str(path / "complex.inpcrd")
        )
        binder = pmd.load_file(str(path / "binder.prmtop"))
        target = pmd.load_file(str(path / "target.prmtop"))
        if len(dry.residues) != 205:
            raise RuntimeError(f"rank {rank:02d} dry residue count differs")
        if len(binder.residues) != 76 or len(target.residues) != 129:
            raise RuntimeError(
                f"rank {rank:02d} split topology residue counts differ"
            )
        if len(binder.atoms) + len(target.atoms) != len(dry.atoms):
            raise RuntimeError(
                f"rank {rank:02d} split topology atom counts differ"
            )

        protein = PDBFile(str(path / "protein_protonated.pdb"))
        chain_lengths = [
            sum(1 for _ in chain.residues())
            for chain in protein.topology.chains()
        ]
        if chain_lengths != [76, 129]:
            raise RuntimeError(
                f"rank {rank:02d} protein PDB chain lengths differ: "
                f"{chain_lengths}"
            )

        prmtop = AmberPrmtopFile(str(path / "solvated.prmtop"))
        inpcrd = AmberInpcrdFile(str(path / "solvated.inpcrd"))
        if prmtop.topology.getNumAtoms() != metadata["solvated_atoms"]:
            raise RuntimeError(f"rank {rank:02d} solvated atom count differs")
        waters = [
            residue
            for residue in prmtop.topology.residues()
            if residue.name in {"WAT", "HOH"}
        ]
        if not waters or any(len(list(residue.atoms())) != 4 for residue in waters):
            raise RuntimeError(f"rank {rank:02d} is not four-site OPC water")
        vectors = inpcrd.boxVectors
        if vectors is None:
            raise RuntimeError(f"rank {rank:02d} has no periodic box")
        lengths_nm = [
            float(
                math.sqrt(
                    sum(
                        float(component) ** 2
                        for component in vector.value_in_unit(unit.nanometer)
                    )
                )
            )
            for vector in vectors
        ]
        angles_deg = [
            vector_angle_degrees(vectors[1], vectors[2]),
            vector_angle_degrees(vectors[0], vectors[2]),
            vector_angle_degrees(vectors[0], vectors[1]),
        ]
        if max(lengths_nm) - min(lengths_nm) > 1e-3:
            raise RuntimeError(f"rank {rank:02d} box edge lengths differ")
        if any(abs(angle - 109.471) > 0.05 for angle in angles_deg):
            raise RuntimeError(
                f"rank {rank:02d} is not a truncated-octahedron box: "
                f"{angles_deg}"
            )
        barostats = sum(
            isinstance(force, MonteCarloBarostat)
            for force in system.getForces()
        )
        if barostats != 1:
            raise RuntimeError(
                f"rank {rank:02d} production barostat count is {barostats}"
            )
        rows.append(
            {
                "rank": rank,
                "atoms": metadata["solvated_atoms"],
                "waters": len(waters),
                "salt_pairs": metadata["salt_pairs"],
                "box_edge_lengths_nm": lengths_nm,
                "box_angles_deg": angles_deg,
                "protein_chain_lengths": chain_lengths,
                "production_barostats": barostats,
                "production_restraint_forces": restraints,
            }
        )
    report = {
        "status": "passed",
        "prepared_systems": len(rows),
        "runtime_openmm": openmm.__version__,
        "openmm_minimum": "8.5",
        "force_field": "AMBER ff19SB",
        "water": "OPC",
        "production_restraints": False,
        "systems": rows,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
