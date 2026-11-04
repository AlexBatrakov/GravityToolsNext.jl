# Runtime and experiment environments

The supported Julia floor is **1.12**, within Julia 1.x. The validated Linux
ARM64 development container pins **1.12.4**. CI tests the 1.12 release line and
rolling Julia 1.x; documentation CI uses 1.12. These are support/test targets,
not claims that every patch or a remote CI job has already been verified locally.

`julia = "1.12"` is a compatibility range `[1.12.0, 2.0.0)`, not an exact
runtime selection. Host Julia and its depot are separate from the container.
Use `scripts/julia-repl.sh` for the container REPL and `scripts/test.sh` for
normal tests. Both activate this package environment. The working directory
may instead be an external experiment directory with an explicitly activated
project. After changing environments, restart old REPL sessions to load the
intended package versions and source code.

Direct external dependencies have compat bounds based on reviewed installed
versions and upstream requirements. Bounds describe compatible version families;
they are not exact experiment pins. Significant package-major upgrades require
their own review and validation. Standard libraries follow Julia/Pkg versioning:
their versions need not match the runtime's minor number. For example, the
validated Julia 1.12.4 environment uses Statistics 1.11.1, Dates/Random 1.11.0
and SHA 0.7.0. Do not rewrite those versions to 1.12 by hand.

## Local package lockfile

The package checkout's `Manifest.toml` remains ignored in public commits. A
coherent local manifest is nevertheless used for day-to-day load/test/smoke
checks from the actual checkout. Fresh checkouts resolve from `Project.toml`.
`Pkg.instantiate()` uses an existing manifest; it does not generally repair
stdlib/dependency locks from another Julia environment.

If a local manifest is incompatible, first preserve its exact bytes in a dated
private backup and record its SHA256. Then use the selected container Julia
and Pkg to regenerate the development lockfile from Project.toml. Compare
resolved versions and run package load, the full default suite and the prepared
engine smoke. Do not use a broad Pkg.update, edit stdlib versions manually or
discard the old lockfile. A disposable clean checkout is still useful to test
resolution, but is not a workaround required for ordinary development.

## Long-lived scientific experiments

Give each long-lived experiment its own environment rather than activating only
the package checkout's development project. Inside the prepared container:

```julia
using Pkg
Pkg.activate("/path/to/experiment/env")
Pkg.develop(path="/path/to/GravityToolsNext.jl")
Pkg.instantiate()
using GravityToolsNext
```

Add packages that the experiment directly imports to that experiment Project.
Archive its Project.toml and Manifest.toml with its results (or in the experiment
repository), separately from the public package repository and private data.
Do not assume a path-developed package is frozen: its source can change without
the manifest changing.

Record Julia/platform versions, depot and active project, the exact manifest
hash, package source commit and dirty-state diff/untracked source hashes,
engine/library/resource identities, input/selection hashes, precision/unit and
phase conventions, RNG implementations/states and actual saved realizations.
The environment alone does not reproduce uncommitted source, private timing
resources or random truth. Inspect `Base.active_project()`, `DEPOT_PATH`,
`pathof(GravityToolsNext)` and `Pkg.status()` before recording a run. Keep host
and container environments distinct and use the same mounted source path.

These policies follow the [Pkg compatibility rules](https://pkgdocs.julialang.org/v1/compatibility/)
and [Pkg environment documentation](https://pkgdocs.julialang.org/v1/environments/).
