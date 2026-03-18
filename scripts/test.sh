#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/docker-dev-common.sh"

gtn_docker_exec 0 julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
