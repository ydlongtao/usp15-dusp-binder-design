# USP15 R7 AF2 model/seed ensemble

## Authorization and rationale

The user authorized continued execution after OVO model 1 failed to recover the
known 6DJ9 positive interface in target-, binder-, and interface-template modes.
R7 tests whether that failure is specific to model 1 or to a single
deterministic prediction.

R7 changes only the AF2 model/seed ensemble. It preserves:

- complete 3T9L A6–134 as the design target;
- the six USP15 hotspots;
- three AF2 recycles;
- iPAE <= 10;
- target-aligned binder RMSD <= 2 Å;
- binder pLDDT >= 80;
- no PyRosetta and no RFD3;
- one V100 with GPU-heavy jobs strictly serial;
- all original USP4/USP11 selectivity rules.

## R7A positive-control calibration

Use OVO's existing `ovo-colabdesign` image and installed weights. Test
interface-template (`ct`) predictions because it was the only R6 mode that
approached the iPAE gate.

Models:

- pTM: `model_1_ptm`, `model_2_ptm`;
- multimer: `model_1_multimer_v3` through `model_5_multimer_v3`.

For every model run dropout-enabled seeds 0, 1, and 2, with three recycles, on:

1. the exact-native 6DJ9 A–K complex;
2. the 6DJ9 UbV pose transplanted onto complete 3T9L A6–134.

A seed passes only when all three unchanged numerical gates pass for both
controls. A model calibrates only when at least two of its three seeds pass.
R7A succeeds only if at least one pTM model and at least one multimer model
calibrate. This prevents a single stochastic result from defining the validator.
Neither positive control is candidate-eligible.

## R7B candidate re-screen

If R7A succeeds, assemble all exact-sequence-unique, complete-target, no-Cys,
45–80-aa candidate-eligible inputs from R3–R5. Recheck the original backbone
hotspot gates before AF2.

Run the selected pTM and multimer models with the same three dropout seeds and
`ct` mode. A design passes the positive screen only when at least two of three
seeds pass every numerical gate for both selected architectures. Label all
positive results as geometry-conditioned.

Cluster passing binders by sequence identity and prioritize different
<=80%-identity clusters, source backbones, and topology families. If fewer than
three pass, a bounded no-Cys LigandMPNN extension may be run on the best eligible
backbones, using the same calibrated validator.

## R7C selectivity and export

For positive designs, prepare USP4 5CTR and USP11 4MEL DUSP-aligned complexes.
Run the original non-interface-conditioned target-template pTM and multimer
off-target tests. For both models and both homologs require:

- delta iPAE >= 5 relative to USP15; and
- off-target iPAE >= 15 or binder RMSD > 4 Å.

Export at least three de-redundant candidates only after positive and off-target
gates pass. Deliver FASTA, input/predicted PDBs, parameters, metrics, provenance,
and elimination reasons. These remain computational candidates, not evidence of
experimental binding or inhibition.

If R7A does not calibrate, stop without screening prior designs or relaxing a
gate.

## Completed R7A result

All 42 expected records completed. `model_2_ptm` passed both controls for all
three seeds:

| control | seed | iPAE | binder RMSD (Å) | binder pLDDT |
| --- | ---: | ---: | ---: | ---: |
| complete 3T9L target + 6DJ9 UbV | 0 | 9.914 | 1.164 | 83.107 |
| complete 3T9L target + 6DJ9 UbV | 1 | 9.938 | 1.196 | 83.012 |
| complete 3T9L target + 6DJ9 UbV | 2 | 9.010 | 1.175 | 85.301 |
| exact-native 6DJ9 | 0 | 9.297 | 1.233 | 87.853 |
| exact-native 6DJ9 | 1 | 9.943 | 1.260 | 86.920 |
| exact-native 6DJ9 | 2 | 9.280 | 1.215 | 88.377 |

`model_1_ptm` passed zero seeds. Every multimer-v3 model from 1 through 5
also passed zero seeds and predicted the binder in a wrong pose, with
target-aligned binder RMSDs of roughly 46–60 Å. R7 therefore failed its
predeclared two-architecture rule and did not run R7B. The calibrated pTM
model is retained only as one branch of the separately authorized R8
independent-predictor workflow.
