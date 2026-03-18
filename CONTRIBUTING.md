# Contributing

Thanks for considering a contribution!

## Quick checklist
- Keep changes focused and well-scoped.
- Add or update tests when the change is testable without external TEMPO/TEMPO2.
- Keep public APIs documented.

## Supported Development Workflow

The supported REPL-oriented baseline is the repository-local Docker workflow in
[docker/dev/README.md](docker/dev/README.md).

First run on a fresh image or fresh Julia depot:

```sh
docker build --platform linux/arm64 -t gravitytoolsnext/dev:baseline -f docker/dev/Dockerfile .
scripts/test.sh
scripts/julia-repl.sh -e 'using GravityToolsNext; println(VERSION); println(Base.active_project())'
scripts/smoke-run.sh
```

Day-to-day commands after that:

- `scripts/dev-shell.sh` opens the long-lived dev container shell
- `scripts/julia-repl.sh` starts `julia --project=.`
- `scripts/test.sh` runs `Pkg.instantiate()` and `Pkg.test()`
- `scripts/smoke-run.sh` runs `Pkg.instantiate()` and the repo-local DDSTG smoke check

`Manifest.toml` is intentionally local-only and is not part of the supported
baseline. See [docker/dev/README.md](docker/dev/README.md) for the rationale,
override knobs, and expected non-bugs such as the reused container and host
UID/GID shell prompt.

## Direct Local Julia Setup

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Notes about TEMPO/TEMPO2

This package can orchestrate external TEMPO/TEMPO2 runs. CI intentionally runs only unit tests that do not require external binaries.
If you add integration tests that require TEMPO/TEMPO2, gate them behind environment checks (so they are skipped on CI).
