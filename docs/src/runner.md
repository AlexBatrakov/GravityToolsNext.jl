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

## Preparation and TIM identity

`MaterializedJob` owns workspace layout and staged paths. Internal preparation
then writes the transformed `_init.par` and freezes the command and input hashes;
these are separate responsibilities, with no new task or public execution API.

The prepared TIM path is the absolute **logical pathname** of the command's TIM
argument relative to its CWD. The runner hashes the file reached through that
pathname, rather than an earlier resolved symlink target. Before launch it checks
that content again, alongside the prepared PAR hash. The staged TIM archive must
have the same expected bytes, because parsed and Basic consumers read that archive.
An equal-byte regular copy is valid; shared inodes or symlinks are not required.

Missing/unreadable or changed TIM inputs, including replaced/retargeted CWD links
and archive/CWD divergence, raise a validation error before launching the engine.
The runner also rechecks both TIM paths before returning captured results, so an
engine-side input edit cannot become successful analysis against a different
archive. It never silently resynchronizes or overwrites either input. Existing
workspace reuse and cleanup policies still apply during materialization.
This is content validation at preparation/execution boundaries, not filesystem
locking or protection against every concurrent edit after a check.

## Execution
- The engine is invoked in `temp_dir` if provided; otherwise in job root (or `tmp/` for split).
- Capture options add flags (e.g., `-write_residuals`).
- A positive `timeout_s` enforces a POSIX engine deadline, captures partial output,
  and terminates the owned process group with a one-second kill grace. The raw
  return remains `(stdout, stderr, coarse_exit_code)`; timeout is a nonzero result.
  `nothing` keeps the existing no-deadline behavior. Finite deadlines are rejected
  on platforms without POSIX process-group support.

## Cleanup and retention
- Old artifacts can be cleared before run.
- Temporary directories may be kept on success or error based on options.

## Manifests
- `write_manifest` and `manifest_style` remain reserved; the runner does not emit a manifest.
