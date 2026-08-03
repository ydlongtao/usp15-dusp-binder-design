# USP15 R9 target-template calibration results

## Outcome

AF2 `model_2_ptm` target-template-only (`tt`) completed all six expected
positive-control predictions but failed calibration. None of seeds 0, 1, or 2
passed both controls. The R9 fixed-panel screen was therefore not started.

## Metrics

| Control | Seed | iPAE | Binder RMSD (Å) | Binder pLDDT | Pass |
|---|---:|---:|---:|---:|---|
| complete 3T9L + 6DJ9 UbV | 0 | 29.123 | 65.677 | 49.232 | no |
| complete 3T9L + 6DJ9 UbV | 1 | 29.141 | 52.439 | 53.341 | no |
| complete 3T9L + 6DJ9 UbV | 2 | 29.293 | 53.938 | 39.269 | no |
| exact-native 6DJ9 | 0 | 28.806 | 74.340 | 59.133 | no |
| exact-native 6DJ9 | 1 | 29.195 | 79.023 | 53.002 | no |
| exact-native 6DJ9 | 2 | 28.691 | 75.394 | 36.590 | no |

The unchanged gates were iPAE `<=10`, target-aligned binder RMSD `<=2 Å`, and
binder pLDDT `>=80`. The result is not a near miss and does not support using
target-template-only predictions as either a positive or off-target validator.

## Decision

R9 stopped at its calibration gate without screening candidates or relaxing a
threshold. The separately documented R10 campaign uses only the R7-calibrated
interface-template protocol and labels every resulting output
`geometry-conditioned`.
