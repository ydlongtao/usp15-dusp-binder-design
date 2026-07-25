# USP15 R6 AF2 template calibration and candidate recovery

## Authorization and scope

R3–R5 showed that OVO's original target-template test does not recover the exact
native 6DJ9 USP15 DUSP–UbV complex. The user has authorized a validation-level
calibration, without changing the target or any numerical acceptance threshold.
R6 remains a computational-design workflow and does not establish biochemical
binding, inhibition, or cellular activity.

The following invariants remain unchanged:

- target: complete 3T9L chain A residues 6–134;
- hotspots: A50, A52, A53, A55, A57, and A61 in the source numbering;
- no PyRosetta and no RFD3;
- one V100, with all GPU-heavy jobs strictly serial;
- iPAE <= 10;
- target-aligned binder RMSD <= 2 Å;
- binder pLDDT >= 80;
- candidate binder length 45–80 aa and no Cys;
- USP4/USP11 selectivity gates from the original plan.

## R6A: exact-native positive-control calibration

Run the exact 6DJ9 native complex through four OVO tests:

1. `af2_model_1_ptm_tbt_3rec`
2. `af2_model_1_multimer_tbt_3rec`
3. `af2_model_1_ptm_ct_3rec`
4. `af2_model_1_multimer_ct_3rec`

`tbt` provides the binder structure as a template but removes template
inter-chain contacts. `ct` retains binder and interface template geometry. The
least geometry-conditioned mode for which both pTM and multimer pass all three
unchanged gates is selected, with `tbt` preferred over `ct`.

The exact-native input has an incomplete crystallographic target and is strictly
diagnostic. It can never count toward the requested candidates.

If neither mode passes both model tests, R6 stops because the calibrated
validator still fails its positive control. No numerical gate is relaxed.

If the exact-native control fails and an independent coordinate audit shows that
the incomplete crystallographic target itself is not preserved, one technical
follow-up is permitted before stopping: repeat the same four tests on the
diagnostic 6DJ9 UbV pose transplanted onto the complete 3T9L A6–134 target. This
follow-up changes neither model, recycle count, nor threshold. Its wild-type UbV
contains Cys and remains candidate-ineligible. A mode is calibrated only if both
pTM and multimer pass on this complete-target control.

## R6B: prior-design re-screen

Build an exact-sequence-unique panel from candidate-eligible R3–R5 structures.
Each structure must retain the complete target, satisfy the original backbone
hotspot gates, contain a 45–80-aa chain-A binder, and contain no Cys.

For each panel member:

- run model 1 multimer with the calibrated template mode as the primary positive
  gate;
- run model 1 pTM with the same mode as corroborating evidence;
- retain the original target-template scores as independent reference evidence;
- reject any design that fails a primary numerical gate.

Passing sequences are clustered by pairwise sequence identity. Final candidates
must be exact-sequence unique, and the export preferentially selects members
from different <=80%-identity clusters, backbones, or topology/source families.
If fewer than three independent passers exist, a bounded sequence-design
extension may be performed on the best eligible backbones using the already
authorized no-Cys ProteinMPNN settings. It is validated with exactly the same
calibrated protocol.

## R6C: USP4/USP11 selectivity

For each de-redundant positive candidate, align the human USP4 5CTR and USP11
4MEL DUSP cores to 3T9L and preserve the candidate binder's starting pose.
Standardize binder as chain A and target as chain B.

Run the original, non-interface-conditioned off-target tests:

- `af2_model_1_ptm_tt_3rec`
- `af2_model_1_multimer_tt_3rec`

For each homolog and both tests require:

- delta iPAE relative to the calibrated USP15 result >= 5; and
- off-target iPAE >= 15 or binder RMSD > 4 Å.

This asymmetric test is intentional: `ct` would inject the transplanted
off-target interface and would therefore be unsuitable as a negative fold-and-
dock screen. Positive candidates will be labeled geometry-conditioned whenever
`ct` is the selected R6A mode.

## Deliverables and stopping rule

Export at least three de-redundant candidates only if they pass the calibrated
USP15 gate and both USP4/USP11 tests. Deliver FASTA, input and predicted complex
PDBs, parameters, positive/off-target metrics, provenance, and explicit
elimination reasons.

If R6 cannot produce three such candidates, report non-convergence without
loosening any threshold or substituting another target or predictor.
