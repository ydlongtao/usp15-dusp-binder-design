#!/usr/bin/env python3
"""Prepare one USP15 R10 complex with ff19SB and explicit OPC solvent."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

import openmm
from openmm import MonteCarloBarostat, XmlSerializer, unit
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME, PDBFile
import parmed as pmd
from pdbfixer import PDBFixer


HOTSPOTS = [45, 47, 48, 50, 52, 56]
EXPECTED_BINDER_RESIDUES = 76
EXPECTED_TARGET_RESIDUES = 129
EXPECTED_PROTEIN_RESIDUES = EXPECTED_BINDER_RESIDUES + EXPECTED_TARGET_RESIDUES


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_tleap(script_path: Path, log_path: Path):
    with log_path.open("w") as log:
        subprocess.run(
            ["tleap", "-f", str(script_path)],
            cwd=str(script_path.parent),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run_ante_mmpbsa(
    solvated_prmtop: Path,
    complex_prmtop: Path,
    target_prmtop: Path,
    binder_prmtop: Path,
    log_path: Path,
):
    with log_path.open("w") as log:
        subprocess.run(
            [
                "ante-MMPBSA.py",
                "-p",
                str(solvated_prmtop),
                "-c",
                str(complex_prmtop),
                "-r",
                str(target_prmtop),
                "-l",
                str(binder_prmtop),
                "-s",
                ":WAT,Na+,Cl-",
                "-n",
                ":1-76",
                "--radii",
                "mbondi3",
            ],
            cwd=str(complex_prmtop.parent),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def write_tleap_script(path: Path, lines):
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def count_waters(pdb_path: Path) -> int:
    residues = set()
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() in {"WAT", "HOH"}:
            residues.add((line[21:22], line[22:26], line[26:27]))
    return len(residues)


def add_restraints(system, topology, positions, protein_atom_count, k_default):
    force = openmm.CustomExternalForce(
        "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2"
    )
    force.setName("protein_heavy_atom_positional_restraints")
    force.addGlobalParameter("k", float(k_default))
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    restrained = 0
    for atom, position in zip(topology.atoms(), positions):
        if atom.index >= protein_atom_count:
            break
        if atom.element is None or atom.element.symbol == "H":
            continue
        force.addParticle(
            atom.index, list(position.value_in_unit(unit.nanometer))
        )
        restrained += 1
    system.addForce(force)
    return restrained


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--padding-nm", type=float, default=1.2)
    parser.add_argument("--ionic-strength-m", type=float, default=0.15)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fixer = PDBFixer(filename=str(args.input_pdb))
    fixer.findMissingResidues()
    if fixer.missingResidues:
        raise ValueError(
            f"Missing residues are not rebuilt: {fixer.missingResidues}"
        )
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise ValueError(
            f"Nonstandard residues found: {fixer.nonstandardResidues}"
        )
    fixer.findMissingAtoms()
    missing_atoms = {
        f"{res.chain.id}:{res.id}:{res.name}": [atom.name for atom in atoms]
        for res, atoms in fixer.missingAtoms.items()
    }
    fixer.addMissingAtoms()

    chain_residues = {
        chain.id: list(chain.residues()) for chain in fixer.topology.chains()
    }
    if set(chain_residues) != {"A", "B"}:
        raise ValueError(f"Expected chains A/B, found {sorted(chain_residues)}")
    if len(chain_residues["A"]) != EXPECTED_BINDER_RESIDUES:
        raise ValueError("Chain A is not the expected 76-residue binder")
    if len(chain_residues["B"]) != EXPECTED_TARGET_RESIDUES:
        raise ValueError("Chain B is not the expected 129-residue USP15 DUSP")

    fixed_heavy = args.output_dir / "fixed_heavy.pdb"
    with fixed_heavy.open("w") as handle:
        PDBFile.writeFile(
            fixer.topology, fixer.positions, handle, keepIds=True
        )

    # OpenMM's Modeller cannot generate an OPC solvent box directly.  Use
    # AmberTools/tleap with the canonical ff19SB and OPC leaprc files, then
    # import the resulting AMBER topology into OpenMM.
    preliminary_pdb = args.output_dir / "preliminary_solvated.pdb"
    preliminary_leap = args.output_dir / "tleap_preliminary.in"
    write_tleap_script(
        preliminary_leap,
        [
            "source leaprc.protein.ff19SB",
            "source leaprc.water.opc",
            f"complex = loadPdb {fixed_heavy}",
            f"solvateOct complex OPCBOX {args.padding_nm * 10.0:.3f}",
            f"savePdb complex {preliminary_pdb}",
            "quit",
        ],
    )
    run_tleap(
        preliminary_leap, args.output_dir / "tleap_preliminary.log"
    )
    preliminary_waters = count_waters(preliminary_pdb)
    if preliminary_waters < 1000:
        raise RuntimeError(
            f"Unexpectedly small OPC box: {preliminary_waters} waters"
        )

    # 55.5 mol/L pure water provides a stable molecule-count estimate for
    # salt pairs.  Neutralizing counterions are added separately.
    salt_pairs = int(
        round(preliminary_waters * args.ionic_strength_m / 55.5)
    )

    tleap_dry_prmtop = args.output_dir / "complex_tleap.prmtop"
    dry_prmtop = args.output_dir / "complex.prmtop"
    dry_inpcrd = args.output_dir / "complex.inpcrd"
    solvated_prmtop = args.output_dir / "solvated.prmtop"
    solvated_inpcrd = args.output_dir / "solvated.inpcrd"
    final_leap = args.output_dir / "tleap_final.in"
    write_tleap_script(
        final_leap,
        [
            "source leaprc.protein.ff19SB",
            "source leaprc.water.opc",
            f"complex = loadPdb {fixed_heavy}",
            f"saveAmberParm complex {tleap_dry_prmtop} {dry_inpcrd}",
            f"solvateOct complex OPCBOX {args.padding_nm * 10.0:.3f}",
            "addIonsRand complex Na+ 0",
            "addIonsRand complex Cl- 0",
            f"addIonsRand complex Na+ {salt_pairs}",
            f"addIonsRand complex Cl- {salt_pairs}",
            f"saveAmberParm complex {solvated_prmtop} {solvated_inpcrd}",
            "quit",
        ],
    )
    run_tleap(final_leap, args.output_dir / "tleap_final.log")

    dry_structure = pmd.load_file(
        str(tleap_dry_prmtop), str(dry_inpcrd)
    )
    if len(dry_structure.residues) != EXPECTED_PROTEIN_RESIDUES:
        raise ValueError(
            f"Dry AMBER topology has {len(dry_structure.residues)} residues"
        )
    for index, residue in enumerate(dry_structure.residues):
        residue.chain = "A" if index < EXPECTED_BINDER_RESIDUES else "B"
    protein_pdb = args.output_dir / "protein_protonated.pdb"
    dry_structure.save(str(protein_pdb), overwrite=True)
    prmtop = AmberPrmtopFile(str(solvated_prmtop))
    inpcrd = AmberInpcrdFile(str(solvated_inpcrd))
    if inpcrd.boxVectors is None:
        raise ValueError("Solvated AMBER coordinates have no periodic box")
    solvated_pdb = args.output_dir / "solvated.pdb"
    with solvated_pdb.open("w") as handle:
        PDBFile.writeFile(
            prmtop.topology, inpcrd.positions, handle, keepIds=False
        )

    # Generate component topologies using AmberTools' supported MM/GBSA
    # topology workflow.  Direct ParmEd slicing can leave inconsistent
    # Lennard-Jones type tables in the component prmtops.
    binder_prmtop = args.output_dir / "binder.prmtop"
    target_prmtop = args.output_dir / "target.prmtop"
    run_ante_mmpbsa(
        solvated_prmtop,
        dry_prmtop,
        target_prmtop,
        binder_prmtop,
        args.output_dir / "ante_mmpbsa.log",
    )
    # Fail preparation immediately if any generated topology is unreadable.
    mmpbsa_complex = pmd.load_file(str(dry_prmtop))
    mmpbsa_binder = pmd.load_file(str(binder_prmtop))
    mmpbsa_target = pmd.load_file(str(target_prmtop))
    if (
        len(mmpbsa_complex.residues) != EXPECTED_PROTEIN_RESIDUES
        or len(mmpbsa_binder.residues) != EXPECTED_BINDER_RESIDUES
        or len(mmpbsa_target.residues) != EXPECTED_TARGET_RESIDUES
        or len(mmpbsa_binder.atoms) + len(mmpbsa_target.atoms)
        != len(mmpbsa_complex.atoms)
    ):
        raise RuntimeError("AmberTools MM/GBSA topology split is inconsistent")

    protein_atom_count = len(dry_structure.atoms)
    common = dict(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=5e-4,
    )
    system_nvt = prmtop.createSystem(**common)
    restrained_atoms = add_restraints(
        system_nvt,
        prmtop.topology,
        inpcrd.positions,
        protein_atom_count,
        1000.0,
    )
    system_npt = prmtop.createSystem(**common)
    add_restraints(
        system_npt,
        prmtop.topology,
        inpcrd.positions,
        protein_atom_count,
        1000.0,
    )
    system_npt.addForce(
        MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25)
    )
    system_production = prmtop.createSystem(**common)
    system_production.addForce(
        MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25)
    )

    for name, system in {
        "system_nvt.xml": system_nvt,
        "system_npt.xml": system_npt,
        "system_production.xml": system_production,
    }.items():
        (args.output_dir / name).write_text(
            XmlSerializer.serialize(system), encoding="utf-8"
        )

    production_restraints = [
        force.getName()
        for force in system_production.getForces()
        if isinstance(force, openmm.CustomExternalForce)
    ]
    if production_restraints:
        raise RuntimeError(
            f"Production system contains restraints: {production_restraints}"
        )

    neutralizing_ions = sum(
        1
        for residue in prmtop.topology.residues()
        if residue.name in {"Na+", "Cl-"}
    ) - 2 * salt_pairs
    metadata = {
        "status": "prepared",
        "openmm_version": openmm.__version__,
        "input_pdb": str(args.input_pdb),
        "input_sha256": sha256(args.input_pdb),
        "pdbfixer_added_missing_atoms": missing_atoms,
        "pdbfixer_added_missing_residues": False,
        "parameterization": "AmberTools tleap imported by OpenMM",
        "mmpbsa_topology_generation": (
            "ante-MMPBSA.py with mbondi3 radii; receptor residues 77-205 "
            "and ligand residues 1-76"
        ),
        "tleap_force_fields": [
            "leaprc.protein.ff19SB",
            "leaprc.water.opc",
        ],
        "water_model": "OPC",
        "ph_assignment_note": (
            "Standard tleap protonation; histidine states recorded in "
            "protein_protonated.pdb and must be reviewed before interpretation."
        ),
        "padding_nm": args.padding_nm,
        "box_shape": "truncated_octahedron",
        "ionic_strength_target_m": args.ionic_strength_m,
        "preliminary_waters": preliminary_waters,
        "salt_pairs": salt_pairs,
        "neutralizing_ions": neutralizing_ions,
        "nonbonded_method": "PME",
        "nonbonded_cutoff_nm": 1.0,
        "constraints": "HBonds",
        "timestep_fs": 2.0,
        "temperature_k": 300.0,
        "pressure_bar": 1.0,
        "binder_chain": "A",
        "target_chain": "B",
        "binder_residues": EXPECTED_BINDER_RESIDUES,
        "target_residues": EXPECTED_TARGET_RESIDUES,
        "target_hotspots_prediction_numbering": HOTSPOTS,
        "protein_atoms": protein_atom_count,
        "solvated_atoms": prmtop.topology.getNumAtoms(),
        "restrained_protein_heavy_atoms_during_equilibration": restrained_atoms,
        "production_custom_external_forces": production_restraints,
        "production_restraints": False,
        "files": {
            path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in (
                fixed_heavy,
                protein_pdb,
                solvated_pdb,
                tleap_dry_prmtop,
                dry_prmtop,
                dry_inpcrd,
                solvated_prmtop,
                solvated_inpcrd,
                binder_prmtop,
                target_prmtop,
                args.output_dir / "system_nvt.xml",
                args.output_dir / "system_npt.xml",
                args.output_dir / "system_production.xml",
            )
        },
    }
    (args.output_dir / "preparation.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
