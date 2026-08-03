# USP15 DUSP Binder Design

Computational design and validation of a compact protein binder for the N-terminal DUSP domain of human USP15.

**Project status:** R10 geometry-conditioned design produced ten computational representatives. A separate rank01 association campaign completed three independent 100 ns OpenMM trajectories from an approximately 35 Å separated state. The trajectories did not show a global binder–target association event; therefore this repository does not claim experimental binding, affinity, inhibition, or selectivity.

## Start here

- [Latest English association report](https://ydlongtao.github.io/usp15-dusp-binder-design/USP15_rank01_association_results_en.html)
- [English rank01/seed0 MD visual report](https://ydlongtao.github.io/usp15-dusp-binder-design/USP15_rank01_seed0_MD_visual_report_en.html)
- [R10 complete report — English](https://ydlongtao.github.io/usp15-dusp-binder-design/USP15_R10_complete_report_en.html)
- [Association protocol](docs/USP15_R10_RANK01_ASSOCIATION_PROTOCOL.md)
- [Association summaries and metrics](docs/results/USP15_rank01_association/)

The self-contained association report includes a real trajectory-based MP4/GIF animation, the rank01 sequence, three-run summary statistics, center-of-mass distance traces, hotspot-contact diagnostics, and the limitations of interpreting an unbiased association attempt.

## Scientific scope

The target is human USP15 (UniProt `Q9Y4E8`), using the DUSP region represented by chain A residues 6–134 of PDB `3T9L`. The interface reference is PDB `6DJ9`; the design hotspot set is `A50,A52,A53,A55,A57,A61`. USP4 (`5CTR`) and USP11 (`4MEL`) were used for computational off-target screening in R10.

The project covers:

1. RFdiffusion RFD1 backbone generation.
2. LigandMPNN sequence design with omit-C constraints.
3. AF2/R10 geometry-conditioned ranking and USP4/USP11 computational counter-screening.
4. OpenMM 8.5.2 association and bound-state diagnostics.
5. Reproducible trajectory analysis and self-contained HTML reporting.

It does not include wet-lab validation and does not establish a mechanism of USP15 inhibition or disruption of a native interaction.

## Rank01 association campaign

The rank01 binder sequence is:

```text
MKIKLVFSDGTEVEVEVDPSDTVLELKKKIEELTGYKPEQLLLFHKGKKLEDGKSLTYHGVKEGDTIHVNIVKEEE
```

Simulation settings:

- OpenMM 8.5.2; AMBER ff19SB; OPC explicit water.
- 300 K, 1 bar, 2 fs timestep.
- Three independent 100 ns NPT runs (`seed0`–`seed2`).
- Initial binder–target separation approximately 35 Å.
- No binder–target distance or positional restraints during production.
- Dell V100 OpenCL backend was used because the local CUDA/NVRTC path was incompatible with this runtime; the molecular model and sampling protocol were unchanged.

The analysis reports binder–target center-of-mass distance, minimum heavy-atom distance, Cα contacts, and six hotspot contact occupancies. A global association event was not observed: the center-of-mass distance remained approximately 32–37 Å across all three runs. Local contacts are reported as diagnostics and must not be interpreted as a bound-state affinity measurement.

## R10 design protocol

- Target: USP15 DUSP, 3T9L chain A residues 6–134.
- Hotspot hard filters: `N_contact_hotspots >= 8` and `N_hotspots_on_interface >= 4`.
- Sequence design: LigandMPNN/ProteinMPNN weights, omit `C`, no amino-acid bias.
- R10 positive ranking: calibrated AF2 `model_2_ptm` interface-template protocol.
- Counter-screen: USP4 and USP11 using the fixed selectivity criteria.
- PyRosetta was not used in the current workflow.

R10 candidates are geometry-conditioned computational designs. Their AF2 structures are predictions, not experimental structures, and the association trajectories do not provide an experimental KD, kon, koff, or absolute binding free energy.

## Reproducibility

Relevant scripts are in `scripts/`, including:

- `scripts/minerva/run_rank01_association.py`
- `scripts/analyze_rank01_association.py`
- `scripts/generate_rank01_association_batch_report.py`
- `scripts/make_rank01_association_movie.pml`
- `scripts/render_rank01_association_movie.sh`

The corrected analysis maps DUSP-local hotspot numbers onto the prepared protein-only PDB numbering and preserves the original failed analysis directory for auditability.

The three raw XTC trajectories and checkpoints are approximately 600 MB in total and are retained in the local `local_results/rank01_association_dell/` directory. The repository contains the derived metrics, final PDB files, animation, and self-contained HTML report; raw trajectories are intentionally not committed because each XTC exceeds GitHub's standard 100 MB file limit.

## License

MIT. See [LICENSE](LICENSE).
