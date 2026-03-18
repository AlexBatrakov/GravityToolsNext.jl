# Development Container Baseline

This directory defines the first repository-local baseline container for `GravityToolsNext.jl`.

Scope of the baseline image:

- `linux/arm64` / `aarch64` development container
- pinned `Julia 1.12.4`
- pinned DDSTG-enabled `tempo` fork
- pinned DDSTG-enabled `tempo2` fork
- mandatory DDSTG runtime wiring:
  - `tempo/arm_fix/*` copied into `tempo/src/`
  - archived `data_ddstg` extracted inside the `tempo` tree
  - `${TEMPO2}/data_ddstg -> ${TEMPO}/data_ddstg` symlink recreated
- low-level `pgplot` / `x11` / scientific system libraries retained conservatively until proven unnecessary

Out of scope for this baseline:

- Python
- notebook tooling
- user-facing plotting workflow
- `psrcat`, `psrchive`, `psrxml`, or other optional pulsar tooling layers

## Pinned sources

- `tempo`
  - repo: `https://github.com/AlexBatrakov/tempo-13.103_ddstg`
  - commit: `755a4502825fd49822fccba9bfded556a321b3b0`
- `tempo2`
  - repo: `https://github.com/AlexBatrakov/tempo2`
  - commit: `0dfe8e244949c1e8c6182cdbd2b6abac72d76012`

## Build

From the repository root:

```bash
docker build \
  --platform linux/arm64 \
  -t gravitytoolsnext/dev:baseline \
  -f docker/dev/Dockerfile \
  .
```

The Dockerfile performs two build-time sanity checks:

- `long double` must stay at 16 bytes with 113-bit mantissa
- `tempo2` must see the linked `data_ddstg` tree from `tempo`

## Validation status

This baseline has been validated against `Julia 1.12.4` with:

- a clean `Pkg.instantiate()` from `Project.toml` in a fresh container/depot
- `Pkg.test()` passing (`22 / 22`)
- direct `julia --project=.` package loading
- DDSTG `Tempo()` / `Tempo2()` validation plus one real `Tempo2` smoke run

No `Project.toml` or package-source changes were required for `Julia 1.12.4`.
The clean resolve selected newer dependency versions within the existing compat
ranges.

## First-run bootstrap

For a fresh image or a fresh Julia depot, use this order from the repository
root:

```bash
scripts/test.sh
scripts/julia-repl.sh -e 'using GravityToolsNext; println(VERSION); println(Base.active_project())'
scripts/smoke-run.sh
```

Optional after that:

```bash
scripts/dev-shell.sh
```

Why this order matters:

- `scripts/test.sh` is the first helper that runs `Pkg.instantiate()`
- `scripts/julia-repl.sh` intentionally starts a plain `julia --project=.`
  session and does not modify the environment for you
- `scripts/smoke-run.sh` also runs `Pkg.instantiate()` before the DDSTG smoke
  check

For a clean package-resolution check, point `GTN_DEV_WORKDIR` at a copy of the
repository without `Manifest.toml`.

## Mounted-host workflow defaults

The helper scripts assume the conservative Phase 2 mount strategy:

- mount `/Users:/Users`
- keep the repository on the host
- keep external research data on the host
- preserve existing absolute paths inside the container

Default script behavior:

- image tag: `gravitytoolsnext/dev:baseline`
- container name: `gravitytoolsnext-dev`
- platform: `linux/arm64`
- host mount: `/Users:/Users`
- repo workdir inside the container: the same absolute host path as the repository

The scripts create a long-lived container on first use and then reuse it.

Override knobs:

- `GTN_DEV_IMAGE`
- `GTN_DEV_CONTAINER`
- `GTN_DEV_PLATFORM`
- `GTN_DEV_HOST_MOUNT`
- `GTN_DEV_HOME`
- `GTN_JULIA_DEPOT_PATH`
- `GTN_DEV_WORKDIR`

If you rebuild the image and want the container to pick up the new image, recreate it:

```bash
docker rm -f gravitytoolsnext-dev
```

## Daily commands

- `scripts/dev-shell.sh`
  - open an interactive shell in the long-lived dev container
  - if arguments are passed, run that command inside the same container instead
    of opening `bash`
- `scripts/julia-repl.sh`
  - start `julia --project=.` in the mounted repository
  - assumes the project has already been instantiated in the current depot
- `scripts/test.sh`
  - run `Pkg.instantiate()` and `Pkg.test()` in the container
- `scripts/smoke-run.sh`
  - run `Pkg.instantiate()` and then a small DDSTG-aware `Tempo2` smoke
    validation against the repo-local synthetic fixture

After the first bootstrap, the normal day-to-day loop can be any mix of:

```bash
scripts/dev-shell.sh
scripts/julia-repl.sh
scripts/test.sh
scripts/smoke-run.sh
```

The container stays alive between these commands on purpose so REPL-based work
does not have to recreate the environment each time.

## Manifest.toml policy

`Manifest.toml` is intentionally local-only and is not part of the supported
repository baseline.

The supported baseline is defined by:

- `Project.toml`
- `docker/dev/Dockerfile`
- the pinned `Julia 1.12.4` image build
- the helper scripts in `scripts/`

Rationale:

- `GravityToolsNext.jl` is a package, not an application with a fixed deployed
  environment
- the validated Phase 2 / Phase 3 workflow already reproduces the supported
  baseline from `Project.toml`
- keeping `Manifest.toml` uncommitted avoids resolver churn while the package
  and daily workflow are still being stabilized

Practical consequence:

- a local ignored `Manifest.toml` resolved under another Julia version may emit
  warnings in the container
- `scripts/test.sh` and `scripts/smoke-run.sh` re-instantiate the project in
  the current depot before executing
- `scripts/julia-repl.sh` does not do that bootstrap for you
- if you want to verify a clean resolve, use a repo copy without
  `Manifest.toml`

## Smoke validation defaults

`scripts/smoke-run.sh` defaults to a repo-local synthetic white-noise fixture:

- work dir:
  `/Users/abatrakov/Documents/Work/PhD/GravityToolsNext.jl/test/data/ddstg_smoke`
- `.par`:
  `j1141-6545_ddstg_gr_pdfb1_white_noise.par`
- `.tim`:
  `j1141-6545_ddstg_gr_pdfb1_white_noise.tim`
- output `.par`:
  `j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par`
- job dir:
  `DDSTG_SMOKE`

The default smoke fixture is the shorter single-backend `PDFB1` case. A larger
multi-backend synthetic white-noise fixture is also retained in the same
directory for broader future testing.

These fixture files are synthetic white-noise data derived from a real system,
not raw observational data. The larger copied multi-backend `.tim` fixture also
has its first-column filesystem paths sanitized to basenames so the repository
does not retain private machine paths.

Override with:

- `GTN_SMOKE_WORK_DIR`
- `GTN_SMOKE_PAR`
- `GTN_SMOKE_TIM`
- `GTN_SMOKE_OUT`
- `GTN_SMOKE_JOB_NAME`

## Known non-bugs

The following behavior is expected for the current supported baseline:

- the dev container is long-lived and is reused across shell / REPL / test /
  smoke commands
- rebuilding the image does not replace the existing container automatically;
  recreate it with `docker rm -f gravitytoolsnext-dev`
- the shell prompt may show `I have no name!` because the container runs with
  the host UID/GID to preserve ownership in the mounted repository
- the Docker build currently emits two warnings about undefined
  `LD_LIBRARY_PATH` / `C_INCLUDE_PATH` expansion; these warnings were observed
  during manual validation and are benign for the validated baseline
