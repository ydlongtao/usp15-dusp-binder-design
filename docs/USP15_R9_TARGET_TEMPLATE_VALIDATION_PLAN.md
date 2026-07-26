# USP15 R9 calibrated target-template validation

## Rationale and authorization

R8 Boltz-2 sequence-only calibration technically completed but did not recover
either known 6DJ9 binder pose. The user authorized continued optimization until
at least three computational candidates are obtained. R9 does not reinterpret
the R8 failure or relax any numerical threshold.

R7 established that AF2 `model_2_ptm` with interface templates (`ct`) recovered
both positive controls for 3/3 dropout seeds. R9 tests the less conditioned
target-template-only mode (`tt`) before using it as a common positive and
off-target validator. This removes binder coordinates and template interface
contacts from inference while retaining the target structure template.

## Invariants

- Target: complete USP15 3T9L chain A residues 6–134.
- Hotspots: A50, A52, A53, A55, A57, and A61.
- No PyRosetta or RFD3.
- Three recycles and dropout seeds 0, 1, and 2.
- iPAE `<=10`.
- target-aligned binder RMSD `<=2 Å`.
- binder pLDDT `>=80`.
- Backbone hotspot gates remain unchanged.
- GPU-heavy stages run serially on one V100.

## R9A target-template calibration

Run AF2 `model_2_ptm` in `tt` mode on:

1. exact-native 6DJ9;
2. the 6DJ9 UbV pose transplanted onto complete 3T9L A6–134.

A seed passes only if both controls meet all three fixed gates. The protocol
calibrates only when at least two of three seeds pass. If it does not calibrate,
do not screen candidates with this mode.

## R9B fixed-panel positive screen

If R9A calibrates, run the exact 52-member R7 panel with the same model, mode,
dropout seeds, and recycles. A design passes only if at least two of three seeds
meet every fixed positive gate.

No sequence is called selective or exported at this stage.

## R9C homolog selectivity

Prepare human USP4 5CTR and USP11 4MEL DUSP targets by structural alignment to
the complete USP15 target. For every R9B passer, run the same calibrated
`model_2_ptm` target-template protocol on USP15, USP4, and USP11.

For each homolog, at least two of three seeds must satisfy:

- `delta iPAE >= 5` relative to the paired USP15 seed; and
- off-target `iPAE >=15` or target-aligned binder RMSD `>4 Å`.

The numerical selectivity thresholds are unchanged. Using one calibrated
model/mode across the on-target and off-target comparisons replaces the
uncalibrated model-1/multimer pair, which failed known positive controls and
therefore cannot provide an interpretable negative result.

## R9D de-redundancy and export

Cluster successful binders at 80% sequence identity and retain structural,
length, and source-backbone diversity. Export at least three representatives
only if they pass the complete positive and homolog screens. Deliver FASTA,
input and predicted PDBs, parameters, metrics, provenance, and elimination
reasons. These remain computational candidates requiring experimental
validation.
