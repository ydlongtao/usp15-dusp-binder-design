#!/usr/bin/env python3
"""Analyze one protein-only USP15 OpenMM trajectory."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.lib.distances import distance_array
import mdtraj as md
import numpy as np


HOTSPOTS = [45, 47, 48, 50, 52, 56]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--frame-interval-ps", type=float, default=10.0)
    parser.add_argument("--burn-in-ns", type=float, default=20.0)
    parser.add_argument("--contact-cutoff-a", type=float, default=4.5)
    parser.add_argument("--hotspot-cutoff-a", type=float, default=5.0)
    parser.add_argument("--sasa-stride", type=int, default=100)
    return parser.parse_args()


def chain_selection(chain_id, suffix):
    return f"(chainID {chain_id} or segid {chain_id}) and ({suffix})"


def quantiles(values):
    data = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(data)),
        "median": float(np.median(data)),
        "q05": float(np.quantile(data, 0.05)),
        "q95": float(np.quantile(data, 0.95)),
        "max": float(np.max(data)),
    }


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    universe = mda.Universe(str(args.topology), str(args.trajectory))
    reference = mda.Universe(str(args.topology))
    binder = universe.select_atoms(chain_selection("A", "protein"))
    target = universe.select_atoms(chain_selection("B", "protein"))
    binder_ca = universe.select_atoms(chain_selection("A", "name CA"))
    target_ca = universe.select_atoms(chain_selection("B", "name CA"))
    binder_heavy = universe.select_atoms(chain_selection("A", "not name H*"))
    target_heavy = universe.select_atoms(chain_selection("B", "not name H*"))
    ref_binder_ca = reference.select_atoms(chain_selection("A", "name CA"))
    ref_target_ca = reference.select_atoms(chain_selection("B", "name CA"))

    if len(binder.residues) != 76 or len(target.residues) != 129:
        raise ValueError(
            f"Unexpected chain sizes A={len(binder.residues)} "
            f"B={len(target.residues)}"
        )
    target_hotspot_atoms = {
        hotspot: target.residues[hotspot - 1].atoms.select_atoms(
            "not name H*"
        )
        for hotspot in HOTSPOTS
    }

    # Native residue contacts are defined once from the protonated starting
    # structure, then their occupancy is measured across the trajectory.
    ref_binder_heavy = reference.select_atoms(
        chain_selection("A", "not name H*")
    )
    ref_target_heavy = reference.select_atoms(
        chain_selection("B", "not name H*")
    )
    native_pairs = []
    for bres in ref_binder_heavy.residues:
        batoms = bres.atoms.select_atoms("not name H*")
        for tres in ref_target_heavy.residues:
            tatoms = tres.atoms.select_atoms("not name H*")
            minimum = float(
                distance_array(batoms.positions, tatoms.positions).min()
            )
            if minimum <= args.contact_cutoff_a:
                native_pairs.append((bres.resid, tres.resid, minimum))

    pair_hits = {(b, t): 0 for b, t, _ in native_pairs}
    hotspot_hits = {hotspot: 0 for hotspot in HOTSPOTS}
    rows = []
    ca_positions = []

    for frame_index, timestep in enumerate(universe.trajectory):
        time_ns = (
            (frame_index + 1) * args.frame_interval_ps / 1000.0
        )
        if time_ns <= args.burn_in_ns:
            continue
        align.alignto(
            universe,
            reference,
            select=chain_selection("B", "name CA"),
            weights="mass",
        )
        binder_rmsd = float(
            rmsd(
                binder_ca.positions,
                ref_binder_ca.positions,
                center=False,
                superposition=False,
            )
        )
        target_rmsd = float(
            rmsd(
                target_ca.positions,
                ref_target_ca.positions,
                center=False,
                superposition=False,
            )
        )
        com_distance = float(
            np.linalg.norm(
                binder.center_of_mass() - target.center_of_mass()
            )
        )

        native_contacts_present = 0
        for binder_resid, target_resid, _ in native_pairs:
            batoms = binder_heavy.select_atoms(f"resid {binder_resid}")
            tatoms = target_heavy.select_atoms(f"resid {target_resid}")
            present = (
                distance_array(batoms.positions, tatoms.positions).min()
                <= args.contact_cutoff_a
            )
            if present:
                pair_hits[(binder_resid, target_resid)] += 1
                native_contacts_present += 1

        hotspot_contacts = 0
        for hotspot in HOTSPOTS:
            hatoms = target_hotspot_atoms[hotspot]
            present = (
                distance_array(binder_heavy.positions, hatoms.positions).min()
                <= args.hotspot_cutoff_a
            )
            if present:
                hotspot_hits[hotspot] += 1
                hotspot_contacts += 1

        rows.append(
            {
                "frame": frame_index,
                "time_ns": time_ns,
                "binder_rmsd_a": binder_rmsd,
                "target_rmsd_a": target_rmsd,
                "com_distance_a": com_distance,
                "native_contacts_present": native_contacts_present,
                "native_contact_fraction": (
                    native_contacts_present / len(native_pairs)
                    if native_pairs
                    else 0.0
                ),
                "hotspots_contacting": hotspot_contacts,
            }
        )
        ca_positions.append(
            binder.select_atoms("name CA").positions.copy()
        )

    frame_count = len(rows)
    if frame_count == 0:
        raise RuntimeError(
            "Trajectory contains no post-burn-in frames; lower --burn-in-ns "
            "only for a technical smoke test"
        )

    with (args.output_dir / "per_frame.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with (args.output_dir / "native_contact_occupancy.csv").open(
        "w", newline=""
    ) as handle:
        fieldnames = [
            "binder_resid",
            "target_resid",
            "initial_min_distance_a",
            "occupancy",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for binder_resid, target_resid, initial in native_pairs:
            writer.writerow(
                {
                    "binder_resid": binder_resid,
                    "target_resid": target_resid,
                    "initial_min_distance_a": initial,
                    "occupancy": pair_hits[
                        (binder_resid, target_resid)
                    ]
                    / frame_count,
                }
            )

    ca_array = np.asarray(ca_positions)
    mean_ca = ca_array.mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum((ca_array - mean_ca) ** 2, axis=2), axis=0))
    with (args.output_dir / "binder_ca_rmsf.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["binder_resid", "ca_rmsf_a"])
        for residue, value in zip(binder.residues, rmsf):
            writer.writerow([residue.resid, float(value)])

    # SASA is evaluated every 1 ns by default to keep the analysis bounded.
    trajectory = md.load_xtc(str(args.trajectory), top=str(args.topology))
    burn_in_frames = int(
        args.burn_in_ns * 1000.0 / args.frame_interval_ps
    )
    sampled = trajectory[burn_in_frames :: args.sasa_stride]
    if sampled.n_frames == 0:
        raise RuntimeError("No post-burn-in frames available for SASA")
    complex_sasa = md.shrake_rupley(sampled, mode="atom").sum(axis=1)
    binder_indices = sampled.topology.select("chainid 0")
    target_indices = sampled.topology.select("chainid 1")
    binder_only = sampled.atom_slice(binder_indices)
    target_only = sampled.atom_slice(target_indices)
    binder_sasa = md.shrake_rupley(binder_only, mode="atom").sum(axis=1)
    target_sasa = md.shrake_rupley(target_only, mode="atom").sum(axis=1)
    buried_sasa_a2 = (binder_sasa + target_sasa - complex_sasa) * 100.0

    hbonds = md.baker_hubbard(
        sampled, freq=0.10, periodic=False, exclude_water=True
    )
    interface_hbonds = []
    for donor, hydrogen, acceptor in hbonds:
        donor_chain = sampled.topology.atom(int(donor)).residue.chain.index
        acceptor_chain = sampled.topology.atom(int(acceptor)).residue.chain.index
        if {donor_chain, acceptor_chain} == {0, 1}:
            interface_hbonds.append(
                [int(donor), int(hydrogen), int(acceptor)]
            )

    hotspot_occupancy = {
        str(hotspot): hotspot_hits[hotspot] / frame_count
        for hotspot in HOTSPOTS
    }
    summary = {
        "status": "completed",
        "frames": frame_count,
        "frame_interval_ps": args.frame_interval_ps,
        "burn_in_ns": args.burn_in_ns,
        "duration_ns_observed": rows[-1]["time_ns"],
        "alignment": "USP15 target chain B C-alpha atoms",
        "binder_rmsd_a": quantiles(
            [row["binder_rmsd_a"] for row in rows]
        ),
        "target_rmsd_a": quantiles(
            [row["target_rmsd_a"] for row in rows]
        ),
        "native_contact_fraction": quantiles(
            [row["native_contact_fraction"] for row in rows]
        ),
        "hotspots_contacting": quantiles(
            [row["hotspots_contacting"] for row in rows]
        ),
        "hotspot_contact_occupancy": hotspot_occupancy,
        "hotspots_with_occupancy_ge_0p5": sum(
            value >= 0.5 for value in hotspot_occupancy.values()
        ),
        "buried_sasa_a2": quantiles(buried_sasa_a2),
        "persistent_interface_hbonds_freq_ge_0p10": len(interface_hbonds),
        "persistent_interface_hbond_atom_triplets": interface_hbonds,
        "com_distance_a": quantiles(
            [row["com_distance_a"] for row in rows]
        ),
        "binder_ca_rmsf_a": quantiles(rmsf),
        "interpretation": (
            "Structural stability diagnostics only; not evidence of "
            "experimental binding or inhibition."
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
