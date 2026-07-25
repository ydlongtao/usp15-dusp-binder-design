# Agent Instructions

## Mission

Maintain and execute the USP15 DUSP minibinder design workflow reproducibly. The repository covers computational design and in-silico screening only.

## Non-negotiable constraints

- Do not use PyRosetta in the current protocol.
- Use RFdiffusion RFD1 only; do not silently substitute RFD3.
- Keep the four pools and target hotspots defined in `config/campaign.json`.
- Never relax acceptance thresholds automatically.
- Do not start 1000-backbone scaling until all preview, pilot, and LigandMPNN + AF2 smoke gates pass.
- Run GPU-heavy stages serially on a single-GPU host.
- Never commit credentials, tokens, SSH endpoints, usernames, database files, absolute private paths, model weights, or generated Nextflow work directories.
- Do not make claims of biochemical inhibition, cellular activity, or selectivity without experimental evidence.

## Workflow invariants

- Target input must be 3T9L chain A residues 6–134.
- Required hotspots are `A50,A52,A53,A55,A57,A61`.
- Backbone hard filters are:
  - `N_contact_hotspots >= 8`
  - `N_hotspots_on_interface >= 4`
- LigandMPNN must generate 3 sequences per backbone with temperature `0.1`, omit `C`, and no amino-acid bias.
- AF2 smoke uses `af2_model_1_multimer_tt_3rec`.
- AF2 positive thresholds are iPAE ≤ 10, target-aligned binder RMSD ≤ 2 Å, and binder pLDDT ≥ 80.
- OpenMM energy and buried-polar diagnostics are ranking features, not Rosetta-equivalent metrics.

## Authorized R2 Phase A diagnostic

R1 remains immutable. The user has authorized the documented R2 Phase A diagnostic in
`docs/USP15_PARAMETER_OPTIMIZATION_PLAN.md`:

- Test exactly these existing Complex_beta backbones:
  - `USP15_R1_short_beta/rfdiffusion_55_standardized.pdb`
  - `USP15_R1_short_beta/rfdiffusion_75_standardized.pdb`
  - `USP15_R1_long_beta/rfdiffusion_76_standardized.pdb`
- For each backbone, generate exactly 3 LigandMPNN sequences at temperature `0.05`
  and exactly 3 at temperature `0.10`.
- Omit `C` and use no amino-acid bias.
- Run `af2_model_1_multimer_tt_3rec` for all 18 sequences.
- Preserve the R1 AF2 thresholds without modification.
- Do not start scaling automatically. Record a machine-readable summary and stop for
  review whether Phase A passes or fails.

This authorization applies only to R2 Phase A. It does not change the R1 temperature
or candidate-selection invariants.

## Authorized R2 Phase B/C optimization

Phase A chain/template alignment audit passed for all 18 predictions and Phase A
failed the unchanged AF2 gates. The user subsequently authorized execution of the
documented R2 Phase B/C matrix:

- Run exactly B1, B2, B3, B4, S1, and S2 from `config/r2_phase_b.tsv`, 50
  backbones per condition.
- Evaluate all conditions with the full six target hotspots and preserve
  `N_contact_hotspots >= 8` plus `N_hotspots_on_interface >= 4`.
- Apply the compactness, topology, loop, and contact-density thresholds recorded
  in `docs/USP15_PARAMETER_OPTIMIZATION_PLAN.md`.
- Use the official RFD1 `Complex_Fold_base_ckpt.pt` only for scaffold-guided
  S1/S2. Do not substitute `Complex_base` when the scaffold checkpoint is absent.
- OVO's binder pipeline cannot represent the official auto-contig scaffold PPI
  call. S1/S2 may therefore run via direct Docker CLI, with raw TRB files
  preserved and OVO-compatible metadata copies generated for standardization.
- Select at most five passing backbones per condition. For each selected
  backbone, generate three sequences at temperature `0.05` and three at `0.10`,
  omit `C`, and run only `af2_model_1_multimer_tt_3rec`.
- Keep all GPU-heavy stages serial and do not initiate 1000-backbone scaling
  unless at least one Phase C sequence passes every unchanged AF2 gate.

R2 Phase B/C has now completed. The 300-backbone matrix yielded 20 selected
backbones, and all 120 expected LigandMPNN/AF2 designs completed without binder
Cys or technical failures. Zero of 120 designs passed all three unchanged AF2
gates. Treat R2 as non-converged, do not start scaling, and require a newly
authorized round for any further topology or generation strategy.

## Repository layout

- `config/campaign.json`: authoritative campaign parameters.
- `config/pools.tsv`: shell-friendly four-pool matrix.
- `scripts/prepare_usp15_target.py`: target cleaning and validation.
- `scripts/run_preview.sh`: RFdiffusion preview driver.
- `scripts/run_backbone_pilot.sh`: full 100-backbone pilot and hard filtering.
- `scripts/run_pilot_queue.sh`: serial pilot queue.
- `scripts/run_sequence_af2_smoke.sh`: LigandMPNN and AF2 smoke.
- `scripts/run_smoke_after_pilots.sh`: post-pilot gate.
- `scripts/run_r2_phase_a.sh`: serial R2 Phase A LigandMPNN and AF2 driver.
- `scripts/summarize_r2_phase_a.py`: machine-readable R2 Phase A AF2 gate summary.
- `config/r2_phase_b.tsv`: authorized six-condition R2 backbone matrix.
- `scripts/audit_r2_phase_a_alignment.py`: independent Phase A chain/RMSD audit.
- `scripts/install_rfdiffusion_scaffold_checkpoint.sh`: verified official scaffold
  model installer.
- `scripts/prepare_r2_phase_b_resources.sh`: PyRosetta-free scaffold resource
  preparation.
- `scripts/run_r2_phase_b_backbones.sh`: one-condition R2 backbone and filter
  driver.
- `scripts/run_r2_phase_b_queue.sh`: serial six-condition Phase B queue.
- `scripts/normalize_scaffold_trb.py`: audited scaffold TRB metadata normalization.
- `scripts/build_r2_phase_c_matrix.py`: build Phase C matrix from passing top-5.
- `scripts/run_r2_phase_c.sh`: serial LigandMPNN/AF2 Phase C driver.
- `scripts/summarize_r2_phase_c.py`: unchanged AF2-gate Phase C summary.
- `scripts/validate_backbones.py`: independent contact validation.
- `scripts/summarize_backbone_metrics.py`: pilot summary generation.

## Editing rules

- Keep scripts idempotent and compatible with `set -euo pipefail`.
- Require deployment paths through environment variables; do not hard-code host paths.
- Preserve completed Nextflow stages and use `-resume`.
- Validate shell changes with:

```bash
bash -n scripts/*.sh
```

- Validate Python changes with:

```bash
python -m py_compile scripts/*.py
```

- Validate JSON before committing:

```bash
python -m json.tool config/campaign.json >/dev/null
```

- Run a secret scan before pushing:

```bash
rg -n -i 'token|password|secret|ssh-rsa|BEGIN .*PRIVATE KEY|[0-9]{1,3}(\\.[0-9]{1,3}){3}' .
```

Review every match; documentation terms are acceptable, real values are not.

## Runtime behavior

- Treat nonzero Nextflow exit codes, missing output counts, and missing trace files as hard failures.
- Preserve failed attempts for audit instead of overwriting them.
- A 15-step preview validates inputs and scheduler behavior; formal hotspot acceptance is decided by the full 50-step pilot hard filters.
- If no pilot backbone passes, stop before sequence design.
- If AF2 smoke fails, stop before scaling.
- Record exclusion reasons in machine-readable output.

## Scientific scope

Candidates are hypotheses for experimental testing. Maintain explicit provenance for PDB structures, target residue numbering, model versions, parameters, metrics, and rejection reasons.
