# Roadmap

This file tracks the project-facing plan for `GravityToolsNext.jl`.

## Phase 0: Preserve The Working Baseline

- [x] Keep the current legacy container available as a known-good fallback.
- [x] Record the essential facts about the legacy environment and why it works.
- [x] Avoid changing the legacy setup before a replacement environment is validated.

## Phase 1: Design A Cleaner Development Environment

- [x] Extract the truly required build/runtime pieces from the legacy Dockerfile.
- [x] Decide which `tempo` and `tempo2` sources should be used in the new setup.
- [x] Decide which parts of the legacy Python and plotting stack are still needed.
- [x] Define the target container workflow for shell, Julia REPL, tests, and example runs.

## Phase 2: Rebuild The Development Container

- [x] Build a cleaner container with `Julia`, `tempo`, and `tempo2`.
- [x] Preserve the working Linux floating-point behavior required by the external tools.
- [x] Verify environment variables, data directories, and executable paths.
- [x] Verify the mounted-host workflow for the repository and external data directories.

## Phase 3: Establish The Daily Workflow

- [x] Document the standard commands for entering the container shell.
- [x] Document the standard commands for starting Julia in the project environment.
- [x] Document the standard commands for running tests and example scripts.
- [x] Keep the workflow convenient for interactive REPL-based development.

## Phase 4: Audit The Package

- [ ] Inventory unfinished or partially implemented features.
- [ ] Review API consistency across `TempoFramework` and `AdaptiveGridFramework`.
- [ ] Review task abstractions, settings, and execution flow for cleanup opportunities.
- [ ] Identify outdated examples and documentation gaps.

## Phase 5: Expand Test Coverage

- [ ] Expand pure-Julia unit test coverage.
- [ ] Add parser tests for `tempo` / `tempo2` output handling.
- [ ] Add tests for settings, materialization, and workspace layout behavior.
- [ ] Add gated integration tests for environments where `tempo` / `tempo2` are available.

## Phase 6: Improve Reliability And Maintainability

- [ ] Improve diagnostics for failed external runs.
- [ ] Reduce hidden assumptions in environment-dependent code.
- [ ] Clarify supported workflows and setup requirements in the documentation.
- [ ] Preserve reproducible job layouts and artifact handling while simplifying the codebase.

## Open Questions

- [ ] Which parts of the legacy pulsar software stack are still required for the current workflow?
- [ ] Which custom forks remain necessary, and which can be replaced with upstream sources?
- [ ] Should Python remain part of the default container image, or become optional?
- [ ] Which end-to-end workflows should become the first supported integration checks?
