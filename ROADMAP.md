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

- [x] Inventory unfinished or partially implemented features.
- [x] Review API consistency across `TempoFramework` and `AdaptiveGridFramework`.
- [x] Review task abstractions, settings, and execution flow for cleanup opportunities.
- [x] Identify outdated examples and documentation gaps.

## Phase 5: Expand Test Coverage

- [x] Expand pure-Julia unit test coverage.
- [x] Add parser tests for `tempo` / `tempo2` output handling.
- [x] Add tests for settings, materialization, and workspace layout behavior.
- [x] Add gated integration tests for environments where `tempo` / `tempo2` are available.

## Phase 6: Improve Reliability And Maintainability

- [ ] Improve diagnostics for failed external runs.
- [ ] Reduce hidden assumptions in environment-dependent code.
- [ ] Clarify supported workflows and setup requirements in the documentation.
- [ ] Preserve reproducible job layouts and artifact handling while simplifying the codebase.

## Phase 7: Release And Distribution Readiness

- [ ] Review versioning and compatibility policy for Julia and direct dependencies.
- [ ] Clarify what belongs to the supported public package surface versus internal-only implementation detail.
- [ ] Improve contributor and onboarding documentation for a clean first setup.
- [ ] Prepare the repository for cleaner tagged releases and changelog discipline.

## Phase 8: Performance And Scale

- [ ] Identify performance bottlenecks in parsing, result assembly, and task execution layers.
- [ ] Review memory and filesystem overhead in large job and grid workflows.
- [ ] Improve the cost of repeated workflow execution where caching or reuse is safe.
- [ ] Preserve correctness and reproducibility while improving throughput.

## Phase 9: API Stabilization And Scientific Expansion

- [ ] Stabilize the most valuable public APIs and extension points.
- [ ] Revisit dormant or deferred task ideas only after the active surface is mature.
- [ ] Expand supported scientific workflows on top of the now-tested and documented core.
- [ ] Reassess architectural cleanup opportunities such as stronger module boundaries if they still provide clear value.

## Scientific Workstream: Synthetic Timing Data

The current research priority is reproducible synthetic timing data for validating
noise-parameter estimation. This workstream can advance alongside the remaining
maintenance phases above; it does not imply their completion.

- [x] S0: Verify fixed-model residual evaluation and relevant TEMPO2 noise conventions.
- [x] S1: Implement precision-preserving TIM reading, editing, and writing.
- [x] S2: Implement reproducible TOA selection and configurable observation schedules.
- [ ] S3: Produce reusable noiseless templates with explicit convergence checks.
- [ ] S4: Generate white-noise datasets with known EFAC/EQUAD and saved truth.
- [ ] S5: Add Fourier red-noise simulation with verified normalization and saved coefficients.
- [ ] S6: Establish ensemble recovery experiments, starting with white noise and then temporal red-noise prediction.

The first end-to-end milestone is a white-noise dataset derived from an existing
TIM and a fixed reference PAR. Generated observation schedules are also in scope;
if they require additional investigation, complete the template-based milestone
first and track schedule completion explicitly.

Generation remains independent of estimation. Preserve original observations,
nominal uncertainties, TOA metadata, and the distinction between generation truth
and analysis settings. Use small reproducible examples before larger research
campaigns.

## Open Questions

- [ ] Which parts of the legacy pulsar software stack are still required for the current workflow?
- [ ] Which custom forks remain necessary, and which can be replaced with upstream sources?
- [ ] Should Python remain part of the default container image, or become optional?
- [ ] Which end-to-end workflows should become the first supported integration checks?

## Possible Further Directions

- [ ] Add broader CI layering for fast unit tests, fixture-driven tests, and gated integration checks.
- [ ] Introduce lightweight benchmark or profiling workflows for hot execution paths.
- [ ] Expand white-noise and result-analysis validation once the reliability layer is mature.
- [ ] Explore whether some conceptual package layers should eventually become stronger module boundaries.
- [ ] Consider release engineering improvements such as more formal changelog, release notes, or packaging strategy if the project becomes more broadly reused.
