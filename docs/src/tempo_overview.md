# Tempo Framework — Overview

The Tempo Framework orchestrates TEMPO/TEMPO2 runs in a predictable, reproducible way.
It is built from small, typed settings that describe paths, engine options, input modifications,
what to capture/keep, runtime workspace layout, and logging.

Core ideas:
- You configure a run using `TempoRunSettings` (composed of `RunPaths`, `EngineOptions`,
  `InputModifiers`, `CaptureOptions`, `RetentionOptions`, `WhiteNoiseOptions`,
  `WorkspaceOptions`, `LoggingOptions`).
- Files are materialized into a job workspace before execution, according to `WorkspaceOptions`.
- The engine runs with consistent flags and working directory selection.
- Outputs are parsed into structured Julia results.

## Data flow (typical)
1. Build `TempoRunSettings` (or copy/modify with `copy_with`).
2. Optionally `validate(settings)` to check inputs and `par_output` filename.
3. Clean old artifacts in the job area (optional, `cleanup_before_run=true`).
4. Materialize a job workspace (flat or split layout).
5. Execute TEMPO/TEMPO2 with derived flags (capture residuals if requested).
6. Parse engine outputs into result types.
7. Optionally save artifacts/manifest and clean temporary files per policy.

## Paths and layout
- `par_input` and `tim_input` are file names relative to `work_dir`.
- `par_output` is a filename-only token; the runner writes it under:
  - `job_root/<par_output>` for `layout=:flat`
  - `job_root/output/<par_output>` for `layout=:split`

## WorkspaceOptions (highlights)
- `work_mode = :inplace | :jobdir` and `job_name` control where the job root lives.
- `layout = :flat | :split` and optional `temp_dir` control on-disk structure and execution cwd.
- `link_tim`, `snapshot_par` control how inputs are staged into the job root.
- `cleanup_before_run`, `keep_tmp_on_success`, `keep_tmp_on_error` control cleanup behavior.
- `write_manifest` adds a small record of what was executed.
- Optional I/O mirroring controls (e.g., `io_mirror = :none | :full | (:depth_minus, N)` if enabled in your build) can mirror job directories or outputs up the directory tree.

## Results (high level)
Parsing builds structured results that can include:
- Basic fit/engine outputs and internal iterations
- Residual statistics
- Optional white noise diagnostics (if enabled)

See the dedicated Results page for details.
