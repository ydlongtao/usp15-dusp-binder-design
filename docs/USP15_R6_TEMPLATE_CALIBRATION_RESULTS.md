# USP15 R6 AF2 template-calibration results

## Outcome

R6 completed both authorized positive-control calibrations with no technical
failures. Neither OVO binder-template mode (`tbt`) nor complex/interface-template
mode (`ct`) recovered the known 6DJ9 USP15 DUSP–UbV complex under all three
unchanged gates. Therefore no R3–R5 design was re-screened and no candidate was
promoted.

The fixed gates were:

- iPAE <= 10;
- target-aligned binder RMSD <= 2 Å;
- binder pLDDT >= 80.

All GPU-heavy tests ran serially on one V100. No PyRosetta or RFD3 step was used.

## Exact-native 6DJ9 control

The diagnostic input used the directly resolved 6DJ9 A–K complex, standardized
to binder chain A and target chain B. It is candidate-ineligible because the
crystallographic target is incomplete.

| Test | iPAE | Binder RMSD (Å) | Binder pLDDT | Pass |
|---|---:|---:|---:|---|
| model 1 pTM + tbt, 3 recycles | 26.58 | 45.71 | 85.81 | no |
| model 1 multimer + tbt, 3 recycles | 21.85 | 50.03 | 77.01 | no |
| model 1 pTM + ct, 3 recycles | 12.55 | 9.56 | 85.92 | no |
| model 1 multimer + ct, 3 recycles | 20.64 | 50.21 | 77.84 | no |

An independent C-alpha audit confirmed identical binder and target sequences,
correct A/B chain mapping, and agreement with OVO's reported binder RMSD. For
pTM+ct the predicted target aligned to the input at 0.96 Å but the binder moved
9.56 Å. For multimer+ct the target itself moved 10.25 Å and the binder moved
50.21 Å. The failure is therefore not a chain-swap or metric-reporting artifact.

## Complete-3T9L-target technical control

Because the exact crystal target has one unresolved residue, the same four tests
were repeated on the 6DJ9 UbV pose transplanted onto complete 3T9L A6–134. This
wild-type UbV control contains Cys and was diagnostic only.

| Test | iPAE | Binder RMSD (Å) | Binder pLDDT | Pass |
|---|---:|---:|---:|---|
| model 1 pTM + tbt, 3 recycles | 27.05 | 43.36 | 85.31 | no |
| model 1 multimer + tbt, 3 recycles | 24.04 | 50.75 | 76.65 | no |
| model 1 pTM + ct, 3 recycles | 9.62 | 9.83 | 87.30 | no |
| model 1 multimer + ct, 3 recycles | 24.06 | 50.21 | 77.34 | no |

The complete-target result rules out the unresolved crystallographic residue as
the sole cause. pTM+ct crossed the iPAE and pLDDT gates but remained 7.83 Å
outside the fixed RMSD gate. Multimer+ct failed all three gates.

## Decision

The calibrated-mode selection rule required both pTM and multimer to pass on a
positive control. No mode qualified. Screening prior candidates with `ct` would
inject an interface geometry that the same validator cannot preserve for its
known positive control, so doing so would not support defensible candidate
promotion.

The campaign remains non-converged at the validation layer. Advancing requires
a newly authorized change such as an AF2 model/seed ensemble or an independently
calibrated structure predictor. The target, numerical gates, and USP4/USP11
selectivity requirements need not change, but the new validator must first pass
the same positive-control calibration before any design can count.
