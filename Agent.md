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

## Repository layout

- `config/campaign.json`: authoritative campaign parameters.
- `config/pools.tsv`: shell-friendly four-pool matrix.
- `scripts/prepare_usp15_target.py`: target cleaning and validation.
- `scripts/run_preview.sh`: RFdiffusion preview driver.
- `scripts/run_backbone_pilot.sh`: full 100-backbone pilot and hard filtering.
- `scripts/run_pilot_queue.sh`: serial pilot queue.
- `scripts/run_sequence_af2_smoke.sh`: LigandMPNN and AF2 smoke.
- `scripts/run_smoke_after_pilots.sh`: post-pilot gate.
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
