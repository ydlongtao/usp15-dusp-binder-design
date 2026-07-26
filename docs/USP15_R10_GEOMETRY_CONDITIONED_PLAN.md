# USP15 R10 geometry-conditioned candidate campaign

## Why R10 is a separate campaign

R9 tested AF2 `model_2_ptm` with only the target template. All three seeds
failed both known positive controls, with iPAE near 29 and binder RMSD between
52 and 79 Å. R9 therefore stopped before candidate screening.

The only protocol that has reproducibly recovered both positive controls is
AF2 `model_2_ptm` with binder and interface templates (`ct`): all three R7
dropout seeds passed the unchanged iPAE, binder RMSD, and binder pLDDT gates.
R10 uses this calibrated protocol as a geometry-conditioned compatibility test.
It does not claim sequence-only pose recovery and does not reinterpret the
failed R8 or R9 validators.

## Invariants

- Complete USP15 3T9L A6–134 and the six predefined hotspots.
- No PyRosetta or RFD3.
- Original backbone hotspot gates.
- AF2 three recycles, dropout seeds 0, 1, and 2.
- iPAE `<=10`.
- target-aligned binder RMSD `<=2 Å`.
- binder pLDDT `>=80`.
- GPU-heavy work serial on one V100.

## R10A fixed-panel positive screen

Run the exact 52-member R7 panel with calibrated `model_2_ptm` `ct`. A design
passes only when at least two of three seeds meet all fixed positive gates.
Every passer remains labelled `geometry-conditioned`.

## R10B same-pose homolog challenge

Structurally align human USP4 5CTR and USP11 4MEL DUSP cores to the USP15
target. Transfer each binder pose without changing its coordinates, then run
the same model, template mode, seeds, recycles, and dropout setting.

For each homolog, at least two paired seeds must satisfy:

- on-target minus off-target confidence separation expressed as
  `off-target iPAE - USP15 iPAE >=5`; and
- off-target `iPAE >=15` or binder RMSD `>4 Å`.

This is deliberately a conservative same-pose compatibility challenge: the
interface template is supplied to both the on-target and homolog models. A
homolog is rejected only if it remains compatible despite that identical
conditioning.

## R10C de-redundancy and export

Cluster binders at 80% sequence identity and retain representatives from
different sequence clusters, source backbones, lengths, and topology families.
Export at least three only after both homolog challenges pass. Deliver FASTA,
input/predicted complexes, metrics, provenance, and elimination reasons.

R10 outputs are computational, geometry-conditioned candidates. They are not
evidence of experimental binding, selectivity, inhibition, or cellular
activity.
