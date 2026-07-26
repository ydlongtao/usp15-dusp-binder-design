# USP15 R10 geometry-conditioned candidate results

## Outcome

R10 completed the fixed 52-member screen, same-pose USP4/USP11 challenges,
80% sequence-identity clustering, and ProteinQC. The campaign produced 10
exported representatives from 24 non-redundant sequence clusters. All 10
exported candidates:

- are 76 aa and contain no cysteine;
- passed all three USP15 AF2 dropout seeds;
- passed all six paired USP4/USP11 selectivity seed tests;
- passed the fixed non-PyRosetta interface audit;
- passed sequence-complexity, hydrophobicity, and Protein-Sol gates.

These are geometry-conditioned computational candidates. R10 tests whether a
supplied binder pose remains compatible with USP15 and becomes incompatible
with aligned homolog surfaces. It is not an independent sequence-only
fold-and-dock validation and is not evidence of binding, biochemical
selectivity, inhibition, or cellular activity.

## Funnel

| Stage | Input | Passing |
|---|---:|---:|
| Fixed-panel USP15 screen | 52 | 41 |
| Non-PyRosetta interface audit | 41 | 36 |
| USP4 and USP11 same-pose selectivity | 36 | 28 |
| 80% sequence-identity clusters | 28 | 24 |
| Exported representatives | 24 | 10 |
| ProteinQC | 10 | 10 |

The USP15 screen used AF2 `model_2_ptm`, interface-template conditioning,
three recycles, dropout seeds 0–2, and the unchanged gates iPAE `<=10`,
target-aligned binder RMSD `<=2 Å`, and binder pLDDT `>=80`. A candidate
required at least two of three seeds; every exported candidate passed three of
three.

The interface audit required interface ΔSASA `>=600 Å²`, zero severe clashes
below 1.5 Å, `N_contact_hotspots>=8`, and
`N_hotspots_on_interface>=4`.

For each homolog, at least two paired seeds were required to satisfy
`off-target iPAE - USP15 iPAE >=5` and either off-target iPAE `>=15` or
binder RMSD `>4 Å`. Every exported candidate passed all three paired seeds
against each homolog.

## Exported candidates

`Min ΔiPAE` is the least passing off-target/on-target separation across the
USP4 and USP11 paired-seed tests. USP15 metrics are means over seeds 0–2.

| Rank | Candidate | Min ΔiPAE | ΔSASA (Å²) | USP15 iPAE | Binder RMSD (Å) | Binder pLDDT |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `r3_partial__P5_rank3_rfdiffusion_7_standardized_packed_2_1` | 13.475 | 712.22 | 8.923 | 1.282 | 90.85 |
| 2 | `r3_partial__P10_rank1_rfdiffusion_8_standardized_packed_3_1` | 13.347 | 878.93 | 7.747 | 1.256 | 90.98 |
| 3 | `r3_partial__P10_rank1_rfdiffusion_8_standardized_packed_2_1` | 12.936 | 923.87 | 6.662 | 0.927 | 92.45 |
| 4 | `r3_partial__P15_rank1_rfdiffusion_1_standardized_packed_1_1` | 12.404 | 906.03 | 7.926 | 1.433 | 87.72 |
| 5 | `r3_partial__P5_rank2_rfdiffusion_0_standardized_packed_1_1` | 11.823 | 811.70 | 8.826 | 1.112 | 90.53 |
| 6 | `r3_mpnn__P10_rank1_rfdiffusion_8_standardized_packed_1_1` | 11.649 | 887.68 | 7.993 | 1.168 | 91.19 |
| 7 | `r3_partial__P10_rank1_rfdiffusion_8_standardized_packed_1_1` | 11.486 | 840.61 | 7.886 | 1.662 | 90.64 |
| 8 | `r3_partial__P5_rank1_rfdiffusion_5_standardized_packed_1_1` | 11.266 | 862.28 | 7.654 | 1.412 | 89.25 |
| 9 | `r3_partial__P15_rank3_rfdiffusion_6_standardized_packed_2_1` | 11.071 | 960.24 | 6.466 | 1.009 | 92.91 |
| 10 | `r3_mpnn__P5_rank1_rfdiffusion_5_standardized_packed_3_1` | 10.988 | 944.61 | 7.647 | 0.832 | 90.46 |

The maximum pairwise sequence identity among the 10 exports is 73.68%.

## Sequences

```fasta
>USP15_R10_rank01
MKIKLVFSDGTEVEVEVDPSDTVLELKKKIEELTGYKPEQLLLFHKGKKLEDGKSLTYHGVKEGDTIHVNIVKEEE
>USP15_R10_rank02
MKITVELSDGTKVEVEIDESDTVQKLAEKIGEITGYKPEDLILLYKGKILPRDKKLSEIGIKEGDTIFVNVNKEKE
>USP15_R10_rank03
MKIKVQLSDGTVVDVEADESDTVRRLAEKIEEITGIKAEDMILLYKGRHLPPDKTLEEIGIKEGDVIYVVVREPEK
>USP15_R10_rank04
MKVKIEFNDGRVVELELKPSDTIAQVLAKLEERYGYKGESLTVLHKGKILNKDDTLEDVGVKEGDTLLINVNEPGS
>USP15_R10_rank05
MKITVELSDGTKVTLELDPSDTMAQVKAKIGEVTGYKPESLMLMYKGRVLKDDETLEDVGIKEGDTILVRVIKYPE
>USP15_R10_rank06
MNIKVKLSDGTVITVEAKESDKVLDLKKKIEEKTGIKPEDLILLYKGKILEDDKTLKEFGIKEGDTIHVNEIKKNE
>USP15_R10_rank07
MKIKVELSDGTVVDVEVDESDTVQKLAEKIGEKTGYPPESLTLLYKGRILEPDKTLAEHGIKEGDVVKVVVNEPEA
>USP15_R10_rank08
MTITVKFSDGTEVDVEIDESEKVSELKAKIEEKVGYKPKDLRLLYKGRVLKDDETLEEVGVKDGDTLLATIVKENE
>USP15_R10_rank09
MTIRVRFPSGETVDLELAPSDTARQIKERIEERLGYRAEELVLLYRGRVLADDDTLADVGIEAGAEILARRNLPGE
>USP15_R10_rank10
EKITVLFPDGTKVTVEVPLDATIKELKAKIQEKTGYNPKDLKLIYKGKVLKDDETLKEAGIKPGDEIHATIELEGE
```

## Homolog preparation and ProteinQC

Human USP4 was prepared from 5CTR chain C residues 10–138 and aligned to the
USP15 DUSP core with 1.416 Å CA RMSD. Human USP11 was prepared from 4MEL chain
B residues 35–157 with 1.444 Å CA RMSD.

ProteinQC used fixed gates of average sequence entropy `>=2.5`, GRAVY `<=1.0`,
and Protein-Sol scaled solubility `>=0.30`; 10/10 passed. ESM-IF1 native
sequence probability and PyDSSP three-state composition were retained as
ranking diagnostics only. The audited ESM-IF1 checkpoint was
1,700,450,121 bytes, and the official Protein-Sol archive matched its fixed
byte count and SHA-256.

## Deliverables

The final archive contains:

- FASTA, candidate metrics, elimination reasons, and manifest;
- one design complex plus three USP15, three USP4, and three USP11 predictions
  for each exported candidate;
- positive-screen, interface-audit, selectivity, homolog-preparation, and
  ProteinQC reports;
- the R10 plan and this results document.

Wet-lab validation remains required. A practical first pass should express all
10 representatives, then measure USP15 binding and USP4/USP11 counterspecificity
before testing any mechanistic or cellular effect.
