# Contributing

Thanks for considering a contribution!

## Quick checklist
- Keep changes focused and well-scoped.
- Add or update tests when the change is testable without external TEMPO/TEMPO2.
- Keep public APIs documented.

## Local setup

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
