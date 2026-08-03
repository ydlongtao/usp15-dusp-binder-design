#!/usr/bin/env python3
"""Run one independently seeded, restartable USP15 OpenMM MD replica."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
import sys
import time

import openmm
from openmm import XmlSerializer, unit
from openmm.app import (
    CheckpointReporter,
    PDBFile,
    Simulation,
    StateDataReporter,
    XTCReporter,
)


def deserialize(path: Path):
    return XmlSerializer.deserialize(path.read_text(encoding="utf-8"))


def serialize(path: Path, value):
    path.write_text(XmlSerializer.serialize(value), encoding="utf-8")


def make_integrator(seed):
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(seed)
    return integrator


def set_barostat_seed(system, seed):
    for force in system.getForces():
        if isinstance(force, openmm.MonteCarloBarostat):
            force.setRandomNumberSeed(seed)


def assert_no_production_restraints(system):
    restraint_names = [
        force.getName()
        for force in system.getForces()
        if isinstance(force, openmm.CustomExternalForce)
    ]
    if restraint_names:
        raise RuntimeError(
            f"Production system contains CustomExternalForce: {restraint_names}"
        )


def new_simulation(topology, system, seed, platform, properties):
    integrator = make_integrator(seed)
    return Simulation(topology, system, integrator, platform, properties)


def save_state_and_pdb(simulation, state_path, pdb_path):
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    serialize(state_path, state)
    with pdb_path.open("w") as handle:
        PDBFile.writeFile(
            simulation.topology,
            state.getPositions(),
            handle,
            keepIds=True,
        )
    return state


def run_stage(
    name,
    topology,
    system,
    seed,
    platform,
    properties,
    input_state,
    output_dir,
    steps,
    restraint_k=None,
):
    state_path = output_dir / f"{name}.state.xml"
    pdb_path = output_dir / f"{name}.pdb"
    if state_path.exists() and pdb_path.exists():
        return state_path

    simulation = new_simulation(topology, system, seed, platform, properties)
    simulation.context.setState(deserialize(input_state))
    if restraint_k is not None:
        simulation.context.setParameter("k", float(restraint_k))
    simulation.reporters.append(
        StateDataReporter(
            str(output_dir / f"{name}.csv"),
            max(1, min(steps, 5000)),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            separator=",",
        )
    )
    simulation.step(steps)
    save_state_and_pdb(simulation, state_path, pdb_path)
    # Minerva H100 nodes reject creation of a second OpenMM CUDA Context while
    # the preceding Context is still alive in the same Python process.
    # Explicitly release each completed stage before constructing the next one.
    del simulation
    gc.collect()
    return state_path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--production-ns", type=float, default=100.0)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--device-index", default="0")
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    status_path = args.output_dir / "status.json"
    if status_path.exists():
        previous = json.loads(status_path.read_text())
        if previous.get("status") == "completed":
            print(json.dumps(previous, indent=2))
            return

    pdb = PDBFile(str(args.prepared_dir / "solvated.pdb"))
    preparation = json.loads(
        (args.prepared_dir / "preparation.json").read_text()
    )
    protein_atom_count = int(preparation["protein_atoms"])
    if protein_atom_count <= 0 or protein_atom_count >= pdb.topology.getNumAtoms():
        raise RuntimeError(
            f"Invalid audited protein atom count {protein_atom_count} for "
            f"{pdb.topology.getNumAtoms()} total atoms"
        )
    system_nvt = deserialize(args.prepared_dir / "system_nvt.xml")
    system_npt = deserialize(args.prepared_dir / "system_npt.xml")
    system_production = deserialize(
        args.prepared_dir / "system_production.xml"
    )
    assert_no_production_restraints(system_production)
    set_barostat_seed(system_npt, args.seed + 10000)
    set_barostat_seed(system_production, args.seed + 20000)

    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {
        "DeviceIndex": args.device_index,
        "Precision": "mixed",
    }

    if args.smoke:
        nvt_steps, npt1_steps, npt2_steps = 5000, 5000, 10000
        production_steps = 25000
        production_ns = 0.05
    else:
        nvt_steps = 250000       # 0.5 ns
        npt1_steps = 250000      # 0.5 ns
        npt2_steps = 500000      # 1.0 ns
        production_steps = int(round(args.production_ns * 500000))
        production_ns = args.production_ns

    started = time.time()
    status = {
        "status": "running",
        "seed": args.seed,
        "openmm_version": openmm.__version__,
        "platform": platform.getName(),
        "platform_properties": properties,
        "smoke": args.smoke,
        "production_ns_requested": production_ns,
        "production_restraints": False,
        "started_unix": started,
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n")

    try:
        minimized_state = args.output_dir / "minimized.state.xml"
        minimized_pdb = args.output_dir / "minimized.pdb"
        if not minimized_state.exists():
            simulation = new_simulation(
                pdb.topology,
                system_nvt,
                args.seed,
                platform,
                properties,
            )
            simulation.context.setPositions(pdb.positions)
            simulation.context.setParameter("k", 1000.0)
            simulation.minimizeEnergy(
                tolerance=10.0
                * unit.kilojoule_per_mole
                / unit.nanometer,
                maxIterations=20000 if not args.smoke else 500,
            )
            simulation.context.setVelocitiesToTemperature(
                300.0 * unit.kelvin, args.seed
            )
            save_state_and_pdb(simulation, minimized_state, minimized_pdb)
            del simulation
            gc.collect()

        nvt = run_stage(
            "nvt_0p5ns" if not args.smoke else "nvt_smoke",
            pdb.topology,
            system_nvt,
            args.seed + 1,
            platform,
            properties,
            minimized_state,
            args.output_dir,
            nvt_steps,
            restraint_k=1000.0,
        )
        npt1 = run_stage(
            "npt_0p5ns_k100" if not args.smoke else "npt1_smoke_k100",
            pdb.topology,
            system_npt,
            args.seed + 2,
            platform,
            properties,
            nvt,
            args.output_dir,
            npt1_steps,
            restraint_k=100.0,
        )
        npt2 = run_stage(
            "npt_1ns_k10" if not args.smoke else "npt2_smoke_k10",
            pdb.topology,
            system_npt,
            args.seed + 3,
            platform,
            properties,
            npt1,
            args.output_dir,
            npt2_steps,
            restraint_k=10.0,
        )

        production_checkpoint = args.output_dir / "production.chk"
        production_xtc = args.output_dir / "production_protein.xtc"
        production_log = args.output_dir / "production.csv"
        simulation = new_simulation(
            pdb.topology,
            system_production,
            args.seed + 4,
            platform,
            properties,
        )
        if production_checkpoint.exists():
            simulation.loadCheckpoint(production_checkpoint.read_bytes())
        else:
            simulation.context.setState(deserialize(npt2))
            # A serialized equilibration State carries its accumulated step
            # count and time.  Production is a new stage, so reset both here.
            # A true production checkpoint restart intentionally retains them.
            simulation.context.setStepCount(0)
            simulation.context.setTime(0.0 * unit.picosecond)

        # tleap/OpenMM may assign solvent to reused PDB chain identifiers.
        # The audited AMBER topology always stores the binder and target first,
        # followed by water and ions, so chain-ID selection is unsafe here.
        protein_atom_indices = list(range(protein_atom_count))
        trajectory_interval = 5000     # 10 ps
        log_interval = 50000           # 100 ps
        checkpoint_interval = 500000   # 1 ns
        simulation.reporters.append(
            XTCReporter(
                str(production_xtc),
                trajectory_interval,
                append=production_xtc.exists(),
                # Keep the two protein chains in one continuous coordinate
                # frame for target-aligned interface analysis.
                enforcePeriodicBox=False,
                atomSubset=protein_atom_indices,
            )
        )
        simulation.reporters.append(
            StateDataReporter(
                str(production_log),
                log_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
                remainingTime=True,
                totalSteps=production_steps,
                separator=",",
                append=production_log.exists(),
            )
        )
        simulation.reporters.append(
            CheckpointReporter(
                str(production_checkpoint), checkpoint_interval
            )
        )
        remaining = production_steps - simulation.currentStep
        if remaining > 0:
            simulation.step(remaining)
        if simulation.currentStep != production_steps:
            raise RuntimeError(
                f"Production ended at step {simulation.currentStep}; "
                f"expected {production_steps}"
            )

        final_state = args.output_dir / "production_final.state.xml"
        final_pdb = args.output_dir / "production_final.pdb"
        save_state_and_pdb(simulation, final_state, final_pdb)
        production_checkpoint.write_bytes(
            simulation.context.createCheckpoint()
        )

        status.update(
            {
                "status": "completed",
                "production_steps": production_steps,
                "trajectory_interval_steps": trajectory_interval,
                "trajectory_interval_ps": 10.0,
                "protein_trajectory_atoms": len(protein_atom_indices),
                "completed_unix": time.time(),
                "wall_seconds": time.time() - started,
                "production_restraints": False,
            }
        )
    except Exception as exc:
        status.update(
            {
                "status": "failed",
                "failed_unix": time.time(),
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
        )
        status_path.write_text(
            json.dumps(status, indent=2, sort_keys=True) + "\n"
        )
        raise

    status_path.write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(status, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
