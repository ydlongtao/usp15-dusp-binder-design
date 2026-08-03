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

## Authorized R3 known-scaffold rescue

The user has authorized continued optimization with the explicit objective of
obtaining at least three computational candidates, while retaining the USP15
DUSP target, no-PyRosetta rule, and unchanged gates. R3 is documented in
`docs/USP15_R3_UBV_RESCUE_PLAN.md` and may:

- use the crystallographic 6DJ9 UbV interface as a diagnostic/scaffold source;
- use the no-Cys 1UBQ ubiquitin fold as a stable 76-aa scaffold;
- run bounded RFD1 partial diffusion at `partial_T` 5, 10, and 15 with the
  target fixed, 10 backbones per condition;
- require `N_contact_hotspots >= 8` and
  `N_hotspots_on_interface >= 4`;
- select at most three backbones per partial-diffusion condition and generate
  exactly three LigandMPNN sequences per backbone at temperature `0.1`, omit
  `C`, with no amino-acid bias;
- use only `af2_model_1_multimer_tt_3rec` for the R3 positive gate.

Wild-type UbV/ubiquitin controls and any diagnostic structures with unresolved
target residues cannot count as candidates. Candidate status still requires all
unchanged positive AF2 gates and the original USP4/USP11 selectivity screen.

R3 has now completed 51 technically successful target-template AF2 predictions
with zero combined-gate passes. The experimental 6DJ9 positive control also
failed this protocol. Treat R3 as non-converged, do not scale any R3
distribution, and require explicit user authorization before running binder- or
interface-template AF2 diagnostics or changing the validation protocol. The
complete evidence is in `docs/USP15_R3_RESULTS.md`.

## Authorized R4 crystallographic-pose ensemble

The user has explicitly asked to continue optimization until at least three
computational peptide candidates are obtained while retaining all prior target,
no-PyRosetta, and AF2/selectivity constraints. R4 is documented in
`docs/USP15_R4_POSE_ENSEMBLE_PLAN.md` and may:

- transplant the four directly resolved 6DJ9 target/binder pairs A–K, B–L,
  C–J, and D–H onto complete 3T9L chain A residues 6–134;
- apply only the smallest deterministic three-dimensional de-clashing
  translation that leaves zero interchain atom pairs below 2 Å and at least
  four hotspots within 5 Å;
- preserve the 15 crystallographic interface positions listed in the plan;
- generate exactly three ProteinMPNN sequences per pose at temperature `0.1`,
  omit `C`, with no amino-acid bias;
- run only `af2_model_1_multimer_tt_3rec` and keep every positive gate
  unchanged.

If at least three sequences pass all positive gates, proceed immediately to the
original USP4/USP11 selectivity screen. If none passes, stop this bounded
ensemble branch rather than scaling it or relaxing a gate.

## Authorized R5 AFDesign sequence optimization

R4 completed with zero positive-gate passes. The user's continuing instruction
to obtain at least three candidates authorizes the bounded R5 sequence-only
optimization in `docs/USP15_R5_AFDESIGN_PLAN.md`:

- use the existing OVO `ovo-colabdesign` image; do not modify the OVO Python
  environment or download another AF2 weight set;
- optimize only the R3 P10 rank-1 RFD1 complex, which is the historical
  near-gate design;
- use multimer model 3 and model 4 during AFDesign and reserve multimer model 1
  for the unchanged OVO target-template validation;
- remove Cys and retain the complete target and all six hotspots;
- run one bounded low-iteration smoke before three full AFDesign seeds;
- keep all GPU work serial.

If model 3/4 cannot form the hotspot interface, R5 may run one 20-soft-step
model-1-in-the-loop diagnostic without a binder coordinate or sequence
template. It must be labeled as circular design/validation evidence, may not
count from internal loss, and may only advance through the unchanged OVO
positive and selectivity gates.

AFDesign internal losses are not acceptance gates. Only the unchanged
`af2_model_1_multimer_tt_3rec` metrics and the original USP4/USP11 screen may
promote a final candidate.

R4 and R5 are now complete. R4 produced 12/12 technically valid AF2 records and
R5 produced 4/4; both had zero combined-gate passes. Do not scale either
distribution. The cumulative failure, including the crystallographic 6DJ9
positive control, indicates a validation-protocol blocker rather than a
near-threshold candidate. See `docs/USP15_R4_R5_RESULTS.md`. Any binder-template
diagnostic, alternate predictor, or validation-decision change requires explicit
user authorization.

## Authorized R6 AF2 template calibration

The user has explicitly authorized validation-level calibration and execution
toward at least three computational candidates after the original target-template
test failed to recover the exact 6DJ9 native complex. R6 is documented in
`docs/USP15_R6_TEMPLATE_CALIBRATION_PLAN.md` and may:

- calibrate OVO binder-template (`tbt`) and complex/interface-template (`ct`)
  modes on the exact native 6DJ9 complex using model 1 pTM and model 1 multimer,
  each with three recycles;
- keep iPAE <= 10, target-aligned binder RMSD <= 2 Å, and binder pLDDT >= 80
  unchanged;
- select the least geometry-conditioned mode whose pTM and multimer tests both
  recover the positive control, preferring `tbt` over `ct`;
- use the calibrated mode to re-evaluate complete-target, no-Cys, 45–80-aa
  designs from prior authorized rounds;
- retain the original target-template pTM and multimer tests for the USP4/USP11
  off-target screen.

The native 6DJ9 positive control and any structure with an incomplete target
remain ineligible as candidates. If only `ct` calibrates, final results must be
reported as geometry-conditioned computational candidates rather than independent
fold-and-dock validation. GPU-heavy tests remain strictly serial.

R6 calibration has completed. Neither `tbt` nor `ct` passed all three unchanged
gates for both model 1 pTM and model 1 multimer on the exact-native 6DJ9
positive control. The permitted complete-3T9L-target follow-up also failed. The
best record was complete-target pTM+ct (iPAE 9.62, binder pLDDT 87.30, binder
RMSD 9.83 Å); model 1 multimer+ct remained a wrong-pose prediction. Do not
re-screen prior candidates with either uncalibrated mode. See
`docs/USP15_R6_TEMPLATE_CALIBRATION_RESULTS.md`. Any model ensemble, alternate
predictor, recycle/seed change, or replacement acceptance logic requires a new
explicit authorization.

## Authorized R7 AF2 model/seed ensemble

After R6 failed, the user explicitly authorized continued execution with an AF2
model/seed ensemble. R7 is documented in
`docs/USP15_R7_AF2_ENSEMBLE_PLAN.md` and may:

- use OVO's existing `ovo-colabdesign` image and installed AF2 weights without
  modifying the OVO environment;
- evaluate pTM models 1–2 and multimer-v3 models 1–5 in `ct` mode;
- run dropout-enabled seeds 0, 1, and 2 for every model, with exactly three
  recycles;
- calibrate first on both the exact-native 6DJ9 control and the 6DJ9 UbV pose
  on complete 3T9L A6–134;
- qualify a model only when at least two of three same-numbered seeds pass every
  unchanged gate on both controls;
- require independently qualified pTM and multimer models before candidate
  re-screening.

R7 does not change the target, positive thresholds, off-target thresholds,
candidate eligibility rules, no-PyRosetta rule, or serial single-V100 policy.
Single stochastic hits do not qualify a model or candidate.

R7 calibration completed 42/42 records. `model_2_ptm` calibrated for all three
seeds on both controls, but every multimer-v3 model 1–5 had zero passing seeds.
Therefore R7 did not meet its two-architecture rule and did not re-screen
candidates.

## Authorized R8 independent Boltz-2 calibration

The user's authorization covered either an AF2 model/seed ensemble or an
independently calibrated predictor. After R7's multimer branch failed, R8 may
use the already installed `ovo-boltz` 2.2.1 image and official MIT-licensed
Boltz-2 weights:

- use sequence-only protein-complex prediction without a complex template,
  forced contact, pocket constraint, or inference potential;
- use the public MSA service only for the already public USP15/UbV controls and
  in-scope designed sequences;
- retain three recycles and seeds 0, 1, and 2;
- compute cross-chain mean PAE, target-aligned binder RMSD, and mean binder
  pLDDT using the same numerical gates;
- calibrate first on both exact-native and complete-target 6DJ9 controls;
- require at least two of three seeds to pass both controls before candidate
  screening.

Downloading official Boltz-2 weights into a campaign-isolated cache is
permitted as a normal R8 implementation step. Require the exact official LFS
byte counts and SHA-256 values recorded in the R8 plan; do not confuse the
transport `xetHash` with the LFS file hash. Do not use Boltz affinity output as
a protein-protein acceptance metric.

## Authorized R10 Minerva migration

The user paused the original-server OpenMM run before any formal production
replica completed. Preserve the following migration invariants:

- do not reconnect to restart the original `usp15-r10-md` tmux queue unless the
  user explicitly reverses the migration decision;
- retain OpenMM 8.5.2, AMBER ff19SB, OPC water, 300 K, 1 bar, 2 fs, and no
  binder-target restraint in production;
- verify the transfer SHA-256 manifest before running on Minerva;
- convert the audited Docker archive to SIF with Apptainer/Singularity;
- require a new Minerva GPU-node CUDA smoke and a passed
  `smoke_minerva/audit.json` before production;
- submit the 30 replicas through Minerva LSF with maximum concurrency 1;
- keep all original failed attempts and paused-state records for provenance.

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
- `docs/USP15_R3_UBV_RESCUE_PLAN.md`: authorized known-scaffold rescue and
  bounded partial-diffusion plan.
- `scripts/prepare_r3_ubv_controls.py`: transplant and validate 6DJ9 UbV
  controls on the complete 3T9L target.
- `scripts/prepare_r3_ubiquitin_scaffold.py`: place a stable 1UBQ scaffold at
  the crystallographic interface.
- `config/r3_partial_diffusion.tsv`: bounded RFD1 partial-diffusion matrix.
- `scripts/run_r3_partial_diffusion.sh`: serial R3 backbone, sequence, and AF2
  pilot driver.
- `scripts/validate_backbones.py`: independent contact validation.
- `scripts/summarize_backbone_metrics.py`: pilot summary generation.
- `docs/USP15_R4_POSE_ENSEMBLE_PLAN.md`: authorized four-pose experimental
  interface ensemble and decision rules.
- `scripts/prepare_r4_6dj9_pose_ensemble.py`: transplant, de-clash, and audit
  four independent 6DJ9 poses.
- `scripts/validate_r4_pose_sequences.py`: enforce R4 sequence counts,
  no-Cys, uniqueness, and fixed-interface invariants.
- `scripts/run_r4_pose_ensemble.sh`: serial ProteinMPNN and AF2 R4 driver.
- `docs/USP15_R5_AFDESIGN_PLAN.md`: bounded AFDesign sequence-only rescue on
  the best RFD1 near-gate complex.
- `scripts/run_r5_afdesign_one.py`: one deterministic AFDesign seed with
  independent validation-model reservation and machine-readable audit.
- `scripts/summarize_r5_afdesign.py`: machine-readable R5 OVO gate summary.
- `docs/USP15_R4_R5_RESULTS.md`: complete R4/R5 metrics and validation blocker.
- `docs/USP15_R6_TEMPLATE_CALIBRATION_PLAN.md`: authorized template calibration
  and unchanged decision rules.
- `scripts/run_r6_template_calibration.sh`: strictly serial R6 calibration
  driver.
- `scripts/summarize_r6_template_calibration.py`: machine-readable R6 gate
  summary.
- `docs/USP15_R6_TEMPLATE_CALIBRATION_RESULTS.md`: exact-native and
  complete-target R6 calibration evidence.
- `docs/USP15_R7_AF2_ENSEMBLE_PLAN.md`: authorized R7 AF2 model/seed
  calibration and decision rule.
- `scripts/run_r7_af2_ensemble_eval.py`: one bounded AF2 architecture/model
  ensemble evaluator.
- `scripts/run_r7_calibration.sh`: strictly serial R7 positive-control driver.
- `scripts/summarize_r7_calibration.py`: R7 control metrics and calibrated-model
  decision.
- `scripts/prepare_r7_rescreen_panel.py`: complete-target, no-Cys, hotspot-gated
  and exact-sequence-unique prior-design panel.
- `docs/USP15_R8_BOLTZ2_CALIBRATION_PLAN.md`: independent sequence-only
  Boltz-2 calibration and candidate decision rules.
- `scripts/download_one_hf_asset.py`: exact-size and LFS-SHA-256 asset
  downloader.
- `scripts/repair_sparse_hf_asset.py`: range-audited sparse asset recovery and
  full LFS SHA-256 verification.
- `scripts/prepare_r8_boltz_controls.py`: sequence-only Boltz positive-control
  YAML generation.
- `scripts/run_r8_boltz2_calibration.sh`: strictly serial seeds 0–2 Boltz
  calibration driver.
- `scripts/summarize_r8_boltz2.py`: independent cross-chain PAE, binder pLDDT,
  and target-aligned binder RMSD audit.
- `scripts/run_r8_af2_candidate_screen.sh`: calibrated pTM-model-2 candidate
  ensemble.
- `scripts/run_r8_boltz2_candidate_screen.sh`: independent Boltz candidate
  ensemble for AF2-positive designs.
- `scripts/run_r8_pipeline_queue.sh`: persistent asset verification, seed-0
  smoke, seeds 1–2 calibration, and serial AF2-to-Boltz positive-screen queue.
- `docker/boltz-v100.Dockerfile`: same Boltz 2.2.1 runtime with an official
  PyTorch/cu121 wheel that retains V100 `sm_70` support.
- `docs/USP15_R10_MINERVA_MIGRATION.md`: paused-state provenance and Minerva
  LSF/Apptainer transfer procedure.
- `scripts/minerva/`: Minerva-specific container conversion, CUDA smoke, and
  strictly serial LSF array submission scripts.
- `scripts/create_minerva_transfer_manifest.py`: SHA-256 write/verify utility
  for the transferred MD directory.

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
