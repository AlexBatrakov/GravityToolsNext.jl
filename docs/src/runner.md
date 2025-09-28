# Runner

This page explains how a run is materialized and executed.

## Materialization
- A job root is chosen based on `WorkspaceOptions`:
  - `:inplace`: job root = `work_dir`
  - `:jobdir`: job root = `work_dir/<job_name>` (auto-generated if missing)
- For `layout = :split`, subfolders `input/`, `output/`, and `tmp/` are created.
- `par_input` and `tim_input` are resolved relative to `work_dir`.
- `par_output` remains a filename-only token; the actual write path is:
  - `job_root/output/<par_output>` for `layout=:split`
  - `job_root/<par_output>` for `:flat`

## Execution
- The engine is invoked in `temp_dir` if provided; otherwise in job root (or `tmp/` for split).
- Capture options add flags (e.g., `-write_residuals`).

## Cleanup and retention
- Old artifacts can be cleared before run.
- Temporary directories may be kept on success or error based on options.

## Manifests
- A minimal manifest can be written with absolute paths of staged inputs and the resolved output path.
