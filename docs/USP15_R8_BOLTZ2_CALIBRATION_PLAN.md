# USP15 R8 independent Boltz-2 calibration

## Purpose

R7 found that AF2 pTM model 2 reproducibly recognizes the 6DJ9 interface, while
all five AF2 multimer-v3 models fail it. The user authorized an independent
predictor path. R8 uses the existing OVO Boltz image and official Boltz-2 model
as an independent sequence-only fold-and-dock validator.

[Boltz code](https://github.com/jwohlwend/boltz) and
[weights](https://huggingface.co/boltz-community/boltz-2) are MIT licensed.
The server image is Boltz 2.2.1. Official weights are stored only in a
campaign-isolated runtime cache and are not committed.

Asset integrity is checked against the official Hugging Face LFS metadata:

- `mols.tar`: 1,855,662,080 bytes, SHA-256
  `39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7`;
- `boltz2_conf.ckpt`: 2,286,561,469 bytes, SHA-256
  `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`.
- `boltz2_aff.ckpt`: 2,062,139,170 bytes, SHA-256
  `dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e`.

Boltz 2.2.1 initializes the affinity checkpoint even when this campaign does
not request or use affinity output. The file is therefore integrity-checked as
a runtime dependency, but affinity remains excluded from all protein-binder
acceptance and ranking decisions.

The Hugging Face `xetHash` is retained only as transport metadata and is not
treated as the file SHA-256. Interrupted or sparse downloads are never promoted
unless both the exact byte count and LFS SHA-256 match. Failed attempts and their
audit reports are preserved under the campaign R8 directory.

## Invariants

- USP15 target and hotspots do not change.
- No PyRosetta or RFD3.
- Three recycles.
- iPAE <= 10.
- target-aligned binder RMSD <= 2 Å.
- binder pLDDT >= 80.
- GPU-heavy jobs strictly serial on one V100.
- USP4/USP11 selectivity requirements remain unchanged.

## R8A smoke and calibration

Prepare sequence-only YAML inputs for:

1. exact-native 6DJ9 binder/target sequences;
2. 6DJ9 UbV and complete 3T9L A6–134 sequences.

Do not provide a structural template, forced pocket/contact constraint, or
inference-time potential. Generate MSAs through the public ColabFold MMseqs2
service. These sequences and structures are public or generated within this
campaign.

Run Boltz-2 with:

- three recycles;
- one diffusion sample per run;
- default 200 sampling steps;
- seeds 0, 1, and 2;
- PDB output;
- full PAE output.

Run seed 0 as a smoke first. Verify chain mapping and independently calculate:

- iPAE as the mean bidirectional cross-chain PAE;
- mean chain-A binder pLDDT;
- chain-B-target-aligned chain-A binder C-alpha RMSD.

Then run seeds 1 and 2. A seed passes only if all three unchanged gates pass on
both controls. Boltz-2 calibrates only if at least two of three seeds pass.

## R8B candidate screen and selectivity

If calibrated, screen the exact-sequence-unique R7 panel with:

- AF2 pTM model 2 ct-mode seeds 0–2, requiring at least 2/3 passes;
- independent Boltz-2 seeds 0–2, requiring at least 2/3 passes.

Boltz predictions are sequence-only and must independently recover the designed
USP15 interface. Positive candidates then enter the original USP4/USP11
target-template pTM/multimer off-target screen. Export only candidates passing
both homolog screens and the existing de-redundancy rules.

If the Boltz positive controls fail, stop without candidate promotion or
threshold relaxation.
