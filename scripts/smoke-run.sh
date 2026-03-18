#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/docker-dev-common.sh"

GTN_SMOKE_WORK_DIR="${GTN_SMOKE_WORK_DIR:-${GTN_REPO_ROOT}/test/data/ddstg_smoke}"
GTN_SMOKE_PAR="${GTN_SMOKE_PAR:-j1141-6545_ddstg_gr_pdfb1_white_noise.par}"
GTN_SMOKE_TIM="${GTN_SMOKE_TIM:-j1141-6545_ddstg_gr_pdfb1_white_noise.tim}"
GTN_SMOKE_OUT="${GTN_SMOKE_OUT:-j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par}"
GTN_SMOKE_JOB_NAME="${GTN_SMOKE_JOB_NAME:-DDSTG_SMOKE}"

if [[ ! -d "${GTN_SMOKE_WORK_DIR}" ]]; then
    echo "Smoke work directory not found: ${GTN_SMOKE_WORK_DIR}" >&2
    exit 1
fi

if [[ ! -f "${GTN_SMOKE_WORK_DIR}/${GTN_SMOKE_PAR}" ]]; then
    echo "Smoke par file not found: ${GTN_SMOKE_WORK_DIR}/${GTN_SMOKE_PAR}" >&2
    exit 1
fi

if [[ ! -f "${GTN_SMOKE_WORK_DIR}/${GTN_SMOKE_TIM}" ]]; then
    echo "Smoke tim file not found: ${GTN_SMOKE_WORK_DIR}/${GTN_SMOKE_TIM}" >&2
    exit 1
fi

echo "Running DDSTG smoke validation in ${GTN_SMOKE_WORK_DIR}" >&2

gtn_docker_exec 0 \
    env \
    GTN_SMOKE_WORK_DIR="${GTN_SMOKE_WORK_DIR}" \
    GTN_SMOKE_PAR="${GTN_SMOKE_PAR}" \
    GTN_SMOKE_TIM="${GTN_SMOKE_TIM}" \
    GTN_SMOKE_OUT="${GTN_SMOKE_OUT}" \
    GTN_SMOKE_JOB_NAME="${GTN_SMOKE_JOB_NAME}" \
    julia --project=. -e 'using Pkg; Pkg.instantiate(); include("docker/dev/smoke_validate.jl")'
