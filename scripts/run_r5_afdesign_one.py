#!/usr/bin/env python3
"""Optimize one fixed RFD1 complex sequence with AFDesign and save an audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from colabdesign import clear_mem, mk_afdesign_model


AA1_TO_3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}
HOTSPOTS = "B50,B52,B53,B55,B57,B61"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument(
        "--design-models",
        default="model_3_multimer_v3,model_4_multimer_v3",
    )
    parser.add_argument("--num-recycles", default=1, type=int)
    parser.add_argument("--soft-iters", default=80, type=int)
    parser.add_argument("--hard-iters", default=16, type=int)
    parser.add_argument("--tries", default=10, type=int)
    parser.add_argument(
        "--binder-template",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--allow-model-1-design", action="store_true")
    return parser.parse_args()


def serializable_log(log: dict[str, object]) -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in log.items():
        if hasattr(value, "item"):
            try:
                output[key] = value.item()
                continue
            except (TypeError, ValueError):
                pass
        if isinstance(value, (str, int, float, bool)) or value is None:
            output[key] = value
        else:
            output[key] = str(value)
    return output


def write_sequence_on_input(
    input_path: Path, output_path: Path, sequence: str
) -> None:
    source = input_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    binder_residues = sorted(
        {
            int(line[22:26])
            for line in source
            if line.startswith("ATOM") and line[21] == "A"
        }
    )
    if binder_residues != list(range(1, 77)):
        raise ValueError("Expected source binder A1-76")
    if len(sequence) != 76 or "C" in sequence:
        raise ValueError("AFDesign sequence must be 76 aa and Cys-free")

    output: list[str] = []
    serial = 0
    for line in source:
        if line.startswith("REMARK"):
            output.append(line)
            continue
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        if chain == "A":
            atom_name = line[12:16].strip()
            if atom_name not in {"N", "CA", "C", "O", "OXT"}:
                continue
            residue = int(line[22:26])
            chars = list(line.ljust(80))
            chars[17:20] = AA1_TO_3[sequence[residue - 1]]
            line = "".join(chars)
        elif chain != "B":
            continue
        serial += 1
        chars = list(line.ljust(80))
        chars[6:11] = f"{serial:5d}"
        output.append("".join(chars))
    output.append("END")
    output_path.write_text("\n".join(output) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    model_names = [
        model.strip() for model in args.design_models.split(",") if model.strip()
    ]
    if "model_1_multimer_v3" in model_names and not args.allow_model_1_design:
        raise ValueError(
            "Reserve model_1_multimer_v3 unless model-in-the-loop design is "
            "explicitly enabled"
        )
    clear_mem()
    model = mk_afdesign_model(
        protocol="binder",
        use_multimer=True,
        num_recycles=args.num_recycles,
        recycle_mode="sample",
        data_dir=args.data_dir,
    )
    unknown_models = sorted(set(model_names) - set(model._model_names))
    if unknown_models:
        raise ValueError(f"Unknown design models: {unknown_models}")
    model.prep_inputs(
        pdb_filename=str(args.input_pdb),
        chain="B",
        binder_chain="A",
        hotspot=HOTSPOTS,
        ignore_missing=False,
        use_binder_template=args.binder_template,
    )
    if model._target_len != 129 or model._binder_len != 76:
        raise ValueError(
            f"Unexpected target/binder lengths {model._target_len}/{model._binder_len}"
        )
    model.restart(seed=args.seed, mode="wt", rm_aa="C")
    model.set_optimizer(
        optimizer="sgd",
        learning_rate=0.1,
        norm_seq_grad=True,
    )
    model.set_weights(
        dgram_cce=0.05,
        plddt=1.0,
        pae=0.5,
        i_pae=5.0,
        con=0.1,
        i_con=5.0,
    )
    model.design_pssm_semigreedy(
        args.soft_iters,
        args.hard_iters,
        tries=args.tries,
        num_recycles=args.num_recycles,
        num_models=len(model_names),
        sample_models=False,
        models=model_names,
        dropout=True,
        verbose=1,
    )
    sequence = model.get_seqs()[0]
    if len(sequence) != 76 or "C" in sequence:
        raise ValueError("AFDesign emitted an invalid binder sequence")

    predicted_path = args.output_dir / f"seed_{args.seed}_afdesign_prediction.pdb"
    input_path = args.output_dir / f"seed_{args.seed}_validation_input.pdb"
    report_path = args.output_dir / f"seed_{args.seed}_afdesign_report.json"
    model.save_pdb(str(predicted_path))
    write_sequence_on_input(args.input_pdb, input_path, sequence)
    report = {
        "phase": "R5 AFDesign fixed-RFD1 sequence optimization",
        "seed": args.seed,
        "source_pdb": args.input_pdb.name,
        "target_chain": "B",
        "binder_chain": "A",
        "target_hotspots": HOTSPOTS.split(","),
        "binder_structure_template_during_design": args.binder_template,
        "binder_sequence_template_during_design": False,
        "sequence": sequence,
        "length": len(sequence),
        "binder_cys_count": sequence.count("C"),
        "design_models": model_names,
        "validation_model_reserved": (
            "model_1_multimer_v3"
            if "model_1_multimer_v3" not in model_names
            else None
        ),
        "model_in_the_loop": "model_1_multimer_v3" in model_names,
        "num_recycles": args.num_recycles,
        "soft_iters": args.soft_iters,
        "hard_iters": args.hard_iters,
        "semigreedy_tries": args.tries,
        "loss_weights": {
            "dgram_cce": 0.05,
            "plddt": 1.0,
            "pae": 0.5,
            "i_pae": 5.0,
            "con": 0.1,
            "i_con": 5.0,
        },
        "final_log": serializable_log(model.aux["log"]),
        "predicted_pdb": predicted_path.name,
        "validation_input_pdb": input_path.name,
    }
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
