# USP15 R3 known-scaffold rescue results

## Outcome

R3 completed 51 technically successful
`af2_model_1_multimer_tt_3rec` predictions across seven bounded diagnostic and
design stages. Zero designs passed all three unchanged positive gates:

- iPAE `<= 10`
- target-aligned binder RMSD `<= 2 Å`
- binder pLDDT `>= 80`

No R3 sequence is a candidate, and no USP4/USP11 selectivity screen was started.

## Stage results

| Stage | AF2 records | Gate passes | Key observation |
|---|---:|---:|---|
| 6DJ9 UbV transplanted to complete 3T9L | 4 | 0 | WT control: iPAE 25.71, RMSD 50.87 Å, pLDDT 42.92 |
| Exact resolved 6DJ9 native-complex diagnostic | 1 | 0 | iPAE 22.56, RMSD 49.82 Å, pLDDT 54.17 |
| Full-backbone LigandMPNN on UbV | 3 | 0 | Binder pLDDT remained 51.30–59.99 |
| Stable 1UBQ scaffold/interface-only design | 4 | 0 | Best designed pLDDT 79.57; iPAE remained 26.10 |
| RFD1 partial diffusion + LigandMPNN | 27 | 0 | 30/30 backbones passed the unchanged hotspot filters |
| RFD1 partial diffusion + ProteinMPNN | 9 | 0 | Switching sequence weights did not recover the interface |
| Fixed 6DJ9 interface + ProteinMPNN core redesign | 3 | 0 | All 15 crystal contact positions were retained; pLDDT remained 58.86–61.56 |

The partial-diffusion matrix used `partial_T` 5, 10, and 15 with 10 backbones
per condition. All 30 backbones passed
`N_contact_hotspots >= 8` and `N_hotspots_on_interface >= 4`. The ranges were:

- P5: `N_contact_hotspots=23–32`, `N_hotspots_on_interface=4–5`
- P10: `N_contact_hotspots=25–34`, `N_hotspots_on_interface=4–5`
- P15: `N_contact_hotspots=28–39`, `N_hotspots_on_interface=4–5`

Nine backbones entered LigandMPNN, producing 27 expected 76-aa, no-Cys
sequences. The closest numerical result had iPAE 11.38 and binder pLDDT 85.30,
but target-aligned binder RMSD was 50.27 Å. Independent structure inspection
showed that its predicted binder moved to a different USP15 surface and no
longer contacted the required hotspot region, so it was rejected.

## Interpretation

The limiting factor is not backbone hotspot contact count. The fixed
target-template AF2 protocol fails to recover even the experimentally
crystallized 6DJ9 UbV complex. Consequently, adding more samples from the same
backbone or sequence distributions is not justified and cannot honestly be
expected to produce three candidates.

This result does not disprove binding of the experimental UbV. It shows that
the current no-binder-template AF2 test is a false negative for this positive
control and is therefore not calibrated for the known USP15 DUSP–UbV system.

## Required decision before R4

R4 requires explicit authorization for a validation-protocol change. The
smallest diagnostic change is to run the crystallographic positive control with
`af2_model_1_multimer_ct_3rec` (binder plus interface template) and
`af2_model_1_multimer_tbt_3rec` (binder template only), without counting either
as a final candidate. If a structure-aware diagnostic recovers 6DJ9, it can be
used for design guidance while retaining a separately documented
template-independent test.

Until that change is authorized, the campaign stops with zero candidates and
does not relax thresholds, start selectivity claims, or promote near-misses.
