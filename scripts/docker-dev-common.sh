#!/usr/bin/env bash
set -euo pipefail

GTN_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GTN_DEV_IMAGE="${GTN_DEV_IMAGE:-gravitytoolsnext/dev:baseline}"
GTN_DEV_CONTAINER="${GTN_DEV_CONTAINER:-gravitytoolsnext-dev}"
GTN_DEV_PLATFORM="${GTN_DEV_PLATFORM:-linux/arm64}"
GTN_DEV_HOST_MOUNT="${GTN_DEV_HOST_MOUNT:-/Users:/Users}"
GTN_DEV_HOME="${GTN_DEV_HOME:-/tmp/gravitytoolsnext-home}"
GTN_JULIA_DEPOT_PATH="${GTN_JULIA_DEPOT_PATH:-/tmp/gravitytoolsnext-julia-depot}"
GTN_DEV_WORKDIR="${GTN_DEV_WORKDIR:-$GTN_REPO_ROOT}"

gtn_require_docker() {
    command -v docker >/dev/null 2>&1 || {
        echo "docker is required but was not found on PATH" >&2
        exit 1
    }

    docker info >/dev/null 2>&1 || {
        cat >&2 <<EOF
docker is installed but the daemon is not reachable.
Start Docker Desktop or restore access to the Docker socket, then try again.
EOF
        exit 1
    }
}

gtn_require_image() {
    if ! docker image inspect "${GTN_DEV_IMAGE}" >/dev/null 2>&1; then
        cat >&2 <<EOF
Baseline image ${GTN_DEV_IMAGE} was not found.
Build it first from ${GTN_REPO_ROOT}:

  docker build --platform ${GTN_DEV_PLATFORM} -t ${GTN_DEV_IMAGE} -f docker/dev/Dockerfile ${GTN_REPO_ROOT}
EOF
        exit 1
    fi
}

gtn_container_exists() {
    docker container inspect "${GTN_DEV_CONTAINER}" >/dev/null 2>&1
}

gtn_container_running() {
    [[ "$(docker container inspect -f '{{.State.Running}}' "${GTN_DEV_CONTAINER}" 2>/dev/null || true)" == "true" ]]
}

gtn_warn_if_stale_container() {
    local image_id container_image_id

    image_id="$(docker image inspect -f '{{.Id}}' "${GTN_DEV_IMAGE}")"
    container_image_id="$(docker container inspect -f '{{.Image}}' "${GTN_DEV_CONTAINER}")"

    if [[ "${image_id}" != "${container_image_id}" ]]; then
        cat >&2 <<EOF
Warning: container ${GTN_DEV_CONTAINER} was created from an older image.
Recreate it to pick up the rebuilt baseline:

  docker rm -f ${GTN_DEV_CONTAINER}
EOF
    fi
}

gtn_prepare_runtime() {
    docker exec \
        -e "HOME=${GTN_DEV_HOME}" \
        -e "JULIA_DEPOT_PATH=${GTN_JULIA_DEPOT_PATH}" \
        "${GTN_DEV_CONTAINER}" \
        bash -lc 'mkdir -p "$HOME/.julia/logs" "$JULIA_DEPOT_PATH"' >/dev/null
}

gtn_ensure_container() {
    gtn_require_docker
    gtn_require_image

    if ! gtn_container_exists; then
        echo "Creating dev container ${GTN_DEV_CONTAINER}" >&2
        docker run -d \
            --platform "${GTN_DEV_PLATFORM}" \
            --name "${GTN_DEV_CONTAINER}" \
            --hostname "${GTN_DEV_CONTAINER}" \
            --init \
            -u "$(id -u):$(id -g)" \
            -e "HOME=${GTN_DEV_HOME}" \
            -e "JULIA_DEPOT_PATH=${GTN_JULIA_DEPOT_PATH}" \
            -e "JULIA_HISTORY=${GTN_DEV_HOME}/.julia/logs/repl_history.jl" \
            -v "${GTN_DEV_HOST_MOUNT}" \
            -w "${GTN_DEV_WORKDIR}" \
            "${GTN_DEV_IMAGE}" \
            sleep infinity >/dev/null
    elif ! gtn_container_running; then
        echo "Starting dev container ${GTN_DEV_CONTAINER}" >&2
        docker start "${GTN_DEV_CONTAINER}" >/dev/null
    fi

    gtn_warn_if_stale_container
    gtn_prepare_runtime
}

gtn_docker_exec() {
    local interactive="$1"
    shift

    gtn_ensure_container

    local -a docker_args
    docker_args=(exec)

    if [[ "${interactive}" == "1" ]]; then
        if [[ -t 0 && -t 1 ]]; then
            docker_args+=(-it)
        else
            docker_args+=(-i)
        fi
    fi

    docker_args+=(
        -e "HOME=${GTN_DEV_HOME}"
        -e "JULIA_DEPOT_PATH=${GTN_JULIA_DEPOT_PATH}"
        -e "JULIA_HISTORY=${GTN_DEV_HOME}/.julia/logs/repl_history.jl"
        -e "TERM=${TERM:-xterm-256color}"
        -w "${GTN_DEV_WORKDIR}"
        "${GTN_DEV_CONTAINER}"
    )

    docker "${docker_args[@]}" "$@"
}
