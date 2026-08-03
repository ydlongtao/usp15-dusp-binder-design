# USP15 R8 Boltz-2 calibration results

## Outcome

R8 stopped at the predeclared seed-0 positive-control gate. Boltz-2 completed
both sequence-only predictions with the official model and verified assets, but
neither control recovered the experimental binder pose within the unchanged
target-aligned binder RMSD threshold. Seeds 1–2, the fixed 52-member candidate
panel, and USP4/USP11 off-target screening were therefore not started.

No R8 prediction is promoted as a candidate.

## Runtime qualification

- Boltz version: 2.2.1.
- Model assets: all three official files passed exact byte-count and LFS
  SHA-256 verification.
- V100 runtime: the derived image retained Boltz 2.2.1 and the official
  weights, while using PyTorch 2.5.1/cu121 with `sm_70` support.
- A real V100 tensor smoke passed before inference.
- Docker shared memory was raised from 64 MB to 8 GB after the first prediction
  attempt exhausted `/dev/shm`. The failed log was preserved before retrying.
- The successful retry produced two PDB, PAE, pLDDT, and confidence records
  with zero failed Boltz examples.

The Boltz output is nested below
`seed_0/boltz_results_input/predictions`; the R8 drivers and independent
summarizers were corrected to audit this actual output location.

## Fixed-gate metrics

| Positive control | iPAE | Target-aligned binder RMSD (Å) | Binder pLDDT | Pass |
|---|---:|---:|---:|---|
| exact-native 6DJ9 | 7.033 | 48.444 | 90.938 | no |
| 6DJ9 UbV pose with complete 3T9L target | 9.897 | 48.481 | 90.910 | no |

The fixed thresholds remained iPAE `<=10`, target-aligned binder RMSD `<=2 Å`,
and binder pLDDT `>=80`. Both controls passed iPAE and binder pLDDT but failed
the pose-recovery requirement by a large margin.

An additional geometry audit confirmed that the complete-target prediction
itself aligned to its reference at 1.378 Å C-alpha RMSD, while the binder
remained 48.481 Å from the reference pose after that target alignment. This
supports a genuine pose-recovery failure rather than a simple chain-label or
target-alignment error.

## Decision

The planned calibration rule required seed 0 to pass both controls before
running seeds 1–2, and at least two of three seeds to pass before candidate
screening. That rule was not met. In accordance with the R8 plan:

- thresholds were not relaxed;
- no structural template, forced contact, or inference potential was added;
- no alternate model was substituted;
- no candidate or selectivity screen was launched.

Further work requires an explicitly authorized new validation strategy. It
cannot be represented as continuation of the calibrated R8 path.
