# Settings

This page documents the core configuration types used to run TEMPO/TEMPO2.

## RunPaths

Holds paths for a TEMPO run.
- `work_dir::String` — absolute working directory
- `par_input::String` — input .par file name (relative to work_dir)
- `tim_input::String` — input .tim file name (relative to work_dir)
- `par_output::String` — output .par filename (no directories)

Use `default_par_output(par_input)` to derive a default out name.

## EngineOptions

Low-level engine configuration:
- `tempo_version::AbstractTempoVersion` (Tempo() or Tempo2())
- `flags::String` (additional CLI flags)
- `nits::Int` (internal iterations ≥ 1)
- `gain::Float64` (> 0)

## InputModifiers
- `override_params::Vector{TempoParameter}`
- `time_start::Union{Nothing,Float64}`
- `time_finish::Union{Nothing,Float64}`

## CaptureOptions
- `write_output::Bool`
- `write_residuals::Bool`

## RetentionOptions
- `save_internal_iterations::Bool`
- `save_residuals::Bool`

## WhiteNoiseOptions
- `enabled::Bool`
- `scope::Symbol` (`:final` | `:all`)

## WorkspaceOptions
Runtime/materialization settings:
- `work_mode` (`:inplace` | `:jobdir`)
- `job_name::Union{Nothing,String}`
- `overwrite` (`:error` | `:overwrite` | `:unique`)
- `layout` (`:flat` | `:split`)
- `temp_dir::Union{Nothing,String}`
- `link_tim::Bool`, `snapshot_par::Bool`
- Cleanup: `cleanup_before_run`, `keep_tmp_on_success`, `keep_tmp_on_error`
- Manifest: `timeout_s`, `write_manifest`, `manifest_style` (`:json` | `:toml`)

## LoggingOptions
- `verbosity` (`:silent | :warn | :info | :debug` or `0..3`)
- `with_timestamps::Bool`

## TempoRunSettings
Top-level settings composed from the types above with a keyword constructor.

```julia
s = TempoRunSettings(
    work_dir   = "/abs/workdir",
    par_input  = "a.par",
    tim_input  = "a.tim",
    par_output = default_par_output("a.par"),
    tempo_version = Tempo2(), nits=2, gain=1.0,
    white_noise_enabled = false,
    work_mode = :jobdir, layout = :split,
)
```
