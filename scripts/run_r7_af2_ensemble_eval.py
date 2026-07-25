#!/usr/bin/env python3
"""Run one OVO/ColabDesign AF2 model across a bounded dropout-seed ensemble."""

from __future__ import annotations

import argparse
import glob
import json
import os
import time
from pathlib import Path

from colabdesign import mk_af_model


METRICS = {
    "rmsd": "target_aligned_binder_rmsd",
    "plddt": "binder_plddt",
    "pae": "binder_pae",
    "ptm": "ptm",
    "con": "con_loss",
    "i_pae": "ipae",
    "i_ptm": "iptm",
    "i_con": "i_con_loss",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--params", required=True)
    parser.add_argument("--architecture", choices=("ptm", "multimer"), required=True)
    parser.add_argument("--model-number", type=int, required=True)
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--num-recycles", type=int, default=3)
    parser.add_argument("--dropout", action="store_true")
    parser.add_argument("--use-binder-template", action="store_true")
    parser.add_argument("--use-interface-template", action="store_true")
    return parser.parse_args()


def total_length(path: str) -> int:
    residues = set()
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                residues.add((line[21], line[22:27]))
    return len(residues)


def flip_chains(pdb_string: str, output_path: Path) -> None:
    binder_lines = []
    target_lines = []
    for line in pdb_string.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if line[21] == "A":
                target_lines.append(line[:21] + "B" + line[22:])
            elif line[21] == "B":
                binder_lines.append(line[:21] + "A" + line[22:])
            else:
                raise ValueError(f"Expected only chains A and B, found {line[21]!r}")
        elif line.strip() and line.startswith(("MODEL", "END")):
            continue
        elif line.strip():
            raise ValueError(f"Unsupported ColabDesign PDB line: {line}")
    output_path.write_text(
        "\n".join(binder_lines + target_lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    seeds = tuple(int(item) for item in args.seeds.split(","))
    if seeds != (0, 1, 2):
        raise ValueError(f"R7 requires exactly seeds 0,1,2, got {seeds}")
    if args.num_recycles != 3:
        raise ValueError("R7 requires exactly three recycles")
    if args.use_interface_template and not args.use_binder_template:
        raise ValueError("Interface template requires binder template")
    if args.architecture == "ptm" and args.model_number not in (1, 2):
        raise ValueError("Template-enabled pTM calibration supports models 1 and 2")
    if args.architecture == "multimer" and args.model_number not in range(1, 6):
        raise ValueError("Multimer calibration supports models 1 through 5")

    model_name = (
        f"model_{args.model_number}_ptm"
        if args.architecture == "ptm"
        else f"model_{args.model_number}_multimer_v3"
    )
    model = mk_af_model(
        protocol="binder",
        data_dir=args.params,
        use_multimer=args.architecture == "multimer",
        model_names=[model_name],
        use_initial_guess=True,
    )

    paths = sorted(glob.glob(os.path.join(args.input_dir, "*.pdb")), key=total_length)
    if not paths:
        raise ValueError(f"No PDB inputs found in {args.input_dir}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = output_dir.with_suffix(".jsonl")

    with jsonl_path.open("w", encoding="utf-8") as handle:
        for path in paths:
            basename = Path(path).stem
            model.prep_inputs(
                path,
                binder_chain="A",
                target_chain="B",
                rm_target=False,
                rm_binder=not args.use_binder_template,
                rm_template_ic=not args.use_interface_template,
            )
            model.set_seq(mode="wildtype")
            model.set_opt(num_recycles=args.num_recycles)

            for seed in seeds:
                started = time.time()
                model.predict(
                    num_models=1,
                    verbose=False,
                    seed=seed,
                    dropout=args.dropout,
                )
                record = {
                    "id": f"{basename}__seed{seed}",
                    "input_id": basename,
                    "architecture": args.architecture,
                    "model_name": model_name,
                    "seed": seed,
                    "dropout": args.dropout,
                    "num_recycles": args.num_recycles,
                    "template_mode": (
                        "ct"
                        if args.use_interface_template
                        else "tbt"
                        if args.use_binder_template
                        else "tt"
                    ),
                }
                record.update(
                    {
                        new_key: model.aux["log"].get(old_key)
                        for old_key, new_key in METRICS.items()
                    }
                )
                record["binder_plddt"] *= 100.0
                record["binder_pae"] *= 31.0
                record["ipae"] *= 31.0
                record["time"] = time.time() - started
                json.dump(record, handle)
                handle.write("\n")
                handle.flush()
                flip_chains(
                    model.save_pdb(),
                    output_dir / f"{record['id']}__{model_name}_ct.pdb",
                )
                print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
