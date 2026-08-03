# USP15 DUSP R3 UbV-scaffold rescue plan

## Objective

Obtain at least three non-cysteine computational USP15 DUSP binder candidates
without changing the target, using PyRosetta, or relaxing the established AF2
and selectivity gates.

## Rationale

R2 completed 120 target-template AF2 predictions without a combined gate pass.
The best iPAE and target-aligned binder RMSD remained far from their thresholds,
so scaling the same de novo distribution is not justified.

PDB 6DJ9 contains the experimentally crystallized USP15 DUSP–UbV 15.D complex.
The resolved UbV fold is 76 residues, within the campaign's 45–80 residue range.
R3 therefore uses the bound UbV as a validated fold/interface scaffold rather
than generating another unconstrained de novo topology.

## Stage 1: positive-control smoke

- Locally align 6DJ9 target chain A residues 25–70 onto the corresponding,
  interface-bearing region of the complete 3T9L A6–134 target reference, then
  transplant 6DJ9 binder chain K residues 1–76 into the 3T9L frame. A local
  alignment is required because the two crystal structures have a large
  C-terminal inter-lobe conformational difference despite matching sequence.
- Standardize the transplanted binder to chain A and the complete 3T9L target
  to chain B. This also avoids the unresolved 6DJ9 target residue 76.
- Test the original UbV sequence as a positive control.
- Test three cysteine-free variants: C11A, C11S, and C11V.
- Run only `af2_model_1_multimer_tt_3rec`.
- Preserve the existing positive gates:
  - iPAE `<= 10`
  - target-aligned binder RMSD `<= 2 Å`
  - binder pLDDT `>= 80`

At least three cysteine-free variants must pass before focused diversification.
The original Cys-containing structure is a control and cannot be counted toward
the requested candidates.

If the transplanted wild-type control fails, run one diagnostic prediction on
the exact resolved 6DJ9 A/K complex with the same `tt_3rec` model and unchanged
gates. This structure lacks resolved target residue A76 and is therefore
diagnostic-only; it cannot be counted as a candidate. Its purpose is to
distinguish failure of the 3T9L pose transfer from failure to recover the native
UbV sequence under the target-template protocol.

If even the exact native complex fails because the binder sequence does not
self-refold, run one minimal LigandMPNN smoke on the transplanted UbV backbone:
exactly three sequences, temperature `0.1`, omit `C`, and no amino-acid bias.
Evaluate those sequences with the same target-template AF2 model and unchanged
gates before any larger sequence sampling.

If full-backbone LigandMPNN still does not recover a stable fold, use PDB 1UBQ
as a no-Cys, experimentally resolved ubiquitin scaffold. Align its conserved
core to the bound UbV core, transfer that pose to 3T9L, remove clashes with a
documented 0.75 Å outward rigid-body translation, and redesign only residues
within 6 Å of USP15. The remaining 1UBQ scaffold sequence stays fixed. Run one
wild-type scaffold control plus exactly three interface designs through the
unchanged AF2 gate before replicating sampling.

If sequence-only interface grafting fails, run a bounded RFD1 partial-diffusion
matrix from the stable 3T9L–ubiquitin pose. Test `partial_T` 5, 10, and 15 at
10 backbones each with contig `76-76/0 B6-134`. Keep the target chain fixed,
retain all six hotspots, and apply `N_contact_hotspots >= 8` plus
`N_hotspots_on_interface >= 4`. Select at most three backbones per condition,
generate three LigandMPNN sequences per backbone at temperature `0.1` with
omit-C, and apply the unchanged target-template AF2 gate.

If the partial-diffusion geometry passes but LigandMPNN sequences fail AF2,
test the originally planned ProteinMPNN weights without FastRelax/PyRosetta.
Use only the top backbone from each partial-diffusion condition, exactly three
sequences per backbone, temperature `0.1`, and omit `C`. The AF2 model and all
acceptance thresholds remain unchanged.

If full-sequence ProteinMPNN fails, preserve the experimentally observed 6DJ9
UbV contact residues (binder positions 4, 6–9, 44, 46, 48–51, and 72–75) and
redesign only the remaining fold-core positions. Run a single three-sequence
smoke at temperature `0.1`, omit `C`, followed by the unchanged AF2 gate before
any expansion.

## Stage 2: focused diversification

If Stage 1 passes, generate a near-native library around the passing Cys-free
scaffold. Keep the ubiquitin fold core and the crystallographic interface
geometry fixed. Diversify only explicitly documented solvent/interface
positions and retain omit-C behavior. Each sequence must independently pass the
same AF2 target-template gate.

## Stage 3: selectivity and finalization

Cluster passing sequences and require at least three non-identical
representatives. Transplant each binder pose onto the aligned USP4 and USP11
DUSP domains and run both target-template AF2 tests from the original plan.
Candidates must retain the USP15 positive gate and satisfy the established
off-target delta criteria. Export FASTA, complex PDB, metrics, provenance, and
rejection reasons.

## Stop rules

- Do not count the Cys-containing control as a candidate.
- Do not relax AF2 or selectivity thresholds.
- Do not claim biochemical inhibition or cellular activity.
- If the crystallographic positive control itself fails the target-template
  protocol, stop and diagnose the evaluation method before generating a larger
  library.
