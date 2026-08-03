#!/usr/bin/env bash
set -euo pipefail

ovo_home_dir="${OVO_HOME_DIR:?Set OVO_HOME_DIR to the initialized OVO home directory}"
model_dir="${ovo_home_dir}/reference_files/rfdiffusion_models"
target="${model_dir}/Complex_Fold_base_ckpt.pt"
temporary="${target}.download"
url="https://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt"
expected_size="483626923"
expected_sha256="0ac3b4024aea811078cec41482528291d6d7d7084bf8190ec118f54642fb81a1"

mkdir -p "${model_dir}"
if [[ ! -s "${target}" ]]; then
    wget \
        --https-only \
        --continue \
        --output-document="${temporary}" \
        "${url}"
    observed_size="$(stat --printf='%s' "${temporary}")"
    if [[ "${observed_size}" != "${expected_size}" ]]; then
        echo "Unexpected checkpoint size: ${observed_size}"
        exit 1
    fi
    observed_sha256="$(sha256sum "${temporary}" | awk '{print $1}')"
    if [[ "${observed_sha256}" != "${expected_sha256}" ]]; then
        echo "Unexpected checkpoint SHA-256: ${observed_sha256}"
        exit 1
    fi
    mv "${temporary}" "${target}"
fi

observed_size="$(stat --printf='%s' "${target}")"
observed_sha256="$(sha256sum "${target}" | awk '{print $1}')"
if [[ "${observed_size}" != "${expected_size}" ]]; then
    echo "Installed checkpoint has unexpected size: ${observed_size}"
    exit 1
fi
if [[ "${observed_sha256}" != "${expected_sha256}" ]]; then
    echo "Installed checkpoint has unexpected SHA-256: ${observed_sha256}"
    exit 1
fi

echo "Verified Complex_Fold_base_ckpt.pt: ${observed_sha256}"
