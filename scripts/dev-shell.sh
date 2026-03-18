#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/docker-dev-common.sh"

if [[ "$#" -eq 0 ]]; then
    gtn_docker_exec 1 bash
else
    gtn_docker_exec 1 "$@"
fi
