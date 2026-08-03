#!/usr/bin/env bash
set -euo pipefail

campaign_dir="${USP15_CAMPAIGN_DIR:?Set USP15_CAMPAIGN_DIR}"
ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR}"
ovo_env_dir="${OVO_ENV_DIR:?Set OVO_ENV_DIR}"
python_bin="${PYTHON_BIN:-${ovo_env_dir}/bin/python}"
pipeline_root="${ovo_env_dir}/lib/python3.13/site-packages/ovo/pipelines"
phase_dir="${R10_QC_DIR:-${campaign_dir}/r10/protein_qc}"
final_dir="${campaign_dir}/r10/final_candidates"
manifest="${final_dir}/manifest.json"
input_dir="${phase_dir}/input"
report_dir="${phase_dir}/reports"
esm_models_dir="${ovo_home_dir}/reference_files/esm_models"
esm_checkpoint="${esm_models_dir}/esm_if1_gvp4_t16_142M_UR50.pt"
proteinsol_zip="${phase_dir}/protein_sol.zip"

mkdir -p "${input_dir}" "${report_dir}" "${phase_dir}/vendor"
exec 9>"${phase_dir}/r10_protein_qc.lock"
if ! flock -n 9; then
    echo "Another R10 ProteinQC driver holds the lock"
    exit 1
fi

if [[ ! -f "${manifest}" ]]; then
    echo "Missing R10 final candidate manifest"
    exit 1
fi
if [[ ! -f "${esm_checkpoint}" ]] || [[ "$(stat -c %s "${esm_checkpoint}")" -lt 100000000 ]]; then
    echo "Missing complete ESM-IF1 checkpoint"
    exit 1
fi

rm -f "${input_dir}"/*.pdb
while IFS= read -r candidate_dir; do
    candidate_name="$(basename "${candidate_dir}")"
    candidate_id="${candidate_name#rank_??_}"
    cp -f "${candidate_dir}/design_complex.pdb" "${input_dir}/${candidate_id}.pdb"
done < <(find "${final_dir}" -mindepth 1 -maxdepth 1 -type d -name 'rank_*' | sort)

expected="$(
    "${python_bin}" - "${manifest}" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1]))["candidate_count"])
PY
)"
observed="$(find "${input_dir}" -maxdepth 1 -type f -name '*.pdb' | wc -l)"
if [[ "${observed}" -ne "${expected}" ]] || [[ "${observed}" -lt 3 ]]; then
    echo "Expected at least three and exactly ${expected} QC inputs, found ${observed}"
    exit 1
fi

"${python_bin}" \
    "${pipeline_root}/proteinqc-seq-composition/bin/seq_composition.py" \
    "${input_dir}" "${report_dir}/seq_composition.csv" --chains A

cp -f \
    "${pipeline_root}/backbone-metrics/bin/pydssp_numpy.py" \
    "${phase_dir}/vendor/pydssp_numpy.py"
docker run --rm \
    -v "${campaign_dir}:${campaign_dir}" \
    -w "${campaign_dir}" \
    ovo-python-structure \
    python3 "${campaign_dir}/scripts/run_r10_pydssp.py" \
    --input-dir "${input_dir}" \
    --pydssp-module "${phase_dir}/vendor/pydssp_numpy.py" \
    --output "${report_dir}/pydssp.csv"

"${python_bin}" "${campaign_dir}/scripts/download_proteinsol.py" \
    --output "${proteinsol_zip}" \
    --report "${report_dir}/proteinsol_download.json"
mkdir -p "${phase_dir}/lib"
unzip -q -o "${proteinsol_zip}" -d "${phase_dir}/lib"
(
    cd "${phase_dir}"
    PYTHONPATH="${pipeline_root}/proteinqc-proteinsol/bin" \
    "${python_bin}" \
    "${pipeline_root}/proteinqc-proteinsol/bin/proteinsol.py" \
    "${input_dir}" "${report_dir}/proteinsol.csv" --chains A
)

cp -f \
    "${pipeline_root}/proteinqc-esmif/bin/esmif.py" \
    "${phase_dir}/vendor/esmif.py"
docker run --rm --gpus all \
    -v "${campaign_dir}:${campaign_dir}" \
    -v "${ovo_home_dir}:${ovo_home_dir}" \
    -w "${campaign_dir}" \
    ovo-esm \
    python3 "${phase_dir}/vendor/esmif.py" \
    "${input_dir}" "${report_dir}/esmif.csv" \
    --chains A \
    --esm_models_dir "${esm_models_dir}" \
    > "${phase_dir}/esmif.log" 2>&1

"${python_bin}" "${campaign_dir}/scripts/summarize_r10_qc.py" \
    --manifest "${manifest}" \
    --seq-composition "${report_dir}/seq_composition.csv" \
    --proteinsol "${report_dir}/proteinsol.csv" \
    --esmif "${report_dir}/esmif.csv" \
    --pydssp "${report_dir}/pydssp.csv" \
    --json "${report_dir}/summary.json" \
    --csv "${report_dir}/metrics.csv"

touch "${phase_dir}/r10_protein_qc.completed"
echo "R10 ProteinQC completed"
