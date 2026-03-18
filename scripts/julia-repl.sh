#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/docker-dev-common.sh"

gtn_docker_exec 1 julia --project=. "$@"
