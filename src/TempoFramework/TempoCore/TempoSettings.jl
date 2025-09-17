# TempoSettings.jl
# Settings and options for running TEMPO/TEMPO2.

# Assumes AbstractTempoVersion, Tempo, Tempo2 are defined in AbstractTempo.jl

# --------------------------------------------------------------------------------------------------------------
# Base type
# --------------------------------------------------------------------------------------------------------------

"""
    AbstractTempoSettings

Abstract base type for all TEMPO settings structures.
Used for dispatching and type hierarchies.
"""
abstract type AbstractTempoSettings end

# --------------------------------------------------------------------------------------------------------------
# Run Paths
# --------------------------------------------------------------------------------------------------------------

# NOTE: Inner constructor normalizes all inputs to `String` so callers can pass
# `AbstractString` (e.g. `SubString`, `OS-specific path` types). Keeping concrete
# `String` fields helps with type stability elsewhere.

"""
    RunPaths

Holds paths for a TEMPO run.

Fields
- `work_dir::String`   : Absolute working directory for the run
- `par_input::String`  : Input `.par` file name or **relative** path (resolved against `work_dir`; absolute paths are rejected)
- `par_output::String` : Output `.par` file **name** (filename only; no directories)
- `tim_input::String`  : `.tim` file name or **relative** path (resolved against `work_dir`; absolute paths are rejected)

Notes
- A keyword constructor `RunPaths(; work_dir, par_input, tim_input, par_output=default_par_output(par_input))` is available.
"""
struct RunPaths
    work_dir::String
    par_input::String
    par_output::String
    tim_input::String

    function RunPaths(work_dir::AbstractString,
                      par_input::AbstractString,
                      par_output::AbstractString,
                      tim_input::AbstractString)
        wd = String(work_dir)
        pi = String(par_input)
        po = String(par_output)
        ti = String(tim_input)
        isabspath(wd)                               || error("RunPaths: work_dir must be an absolute path; got: $wd")
        !occursin(r"\.[Pp][Aa][Rr]$", basename(pi)) && error("RunPaths: par_input must end with .par (case-insensitive); got $pi")
        isabspath(pi)                               && error("RunPaths: par_input must be relative to work_dir (absolute paths are not allowed); got: $pi")
        !occursin(r"\.[Pp][Aa][Rr]$", basename(po)) && error("RunPaths: par_output must end with .par (case-insensitive); got $po")
        basename(po) == po                          || error("RunPaths: par_output must be a file name (no directories); got: $po")
        !occursin(r"\.[Tt][Ii][Mm]$", basename(ti)) && error("RunPaths: tim_input must end with .tim (case-insensitive); got $ti")
        isabspath(ti)                               && error("RunPaths: tim_input must be relative to work_dir (absolute paths are not allowed); got: $ti")
        isempty(pi)                                 && error("RunPaths: par_input must be non-empty")
        isempty(ti)                                 && error("RunPaths: tim_input must be non-empty")
        isempty(po)                                 && error("RunPaths: par_output must be non-empty")
        new(wd, pi, po, ti)
    end
end

# Case-insensitive `.par` handling
_par_stripped(name::AbstractString) = replace(name, r"\.[Pp][Aa][Rr]$" => "")

"""
    default_par_output(par_input::AbstractString) -> String
    
Generate a default output **file name** (no directories) by taking `basename(par_input)`
and replacing its `.par` extension with `*_out.par`.

- Accepts either a bare filename or any path.
- Still enforces that the basename ends with `.par` (case‑insensitive).
"""
function default_par_output(par_input::AbstractString)
    # allow paths: derive name from basename
    name = basename(par_input)
    # enforce .par extension (case-insensitive)
    if !occursin(r"\.[Pp][Aa][Rr]$", name)
        error("RunPaths: par_input must end with .par (case-insensitive); got basename '$name' (from '$par_input')")
    end
    return string(_par_stripped(name), "_out.par")
end

function RunPaths(; work_dir::AbstractString,
                    par_input::AbstractString,
                    tim_input::AbstractString,
                    par_output::AbstractString = default_par_output(par_input))
    return RunPaths(work_dir, par_input, par_output, tim_input)
end

function Base.show(io::IO, ::MIME"text/plain", f::RunPaths)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "RunPaths")
    println(io, spad, "work_dir:        ", f.work_dir)
    println(io, spad, "par_input:       ", f.par_input)
    println(io, spad, "par_output:      ", f.par_output)
    println(io, spad, "tim_input:       ", f.tim_input)
end

# --------------------------------------------------------------------------------------------------------------
# Engine options
# --------------------------------------------------------------------------------------------------------------

# normalize flags into a single string; reject newlines for safety
@inline _normalize_flags(flags::AbstractString) = strip(String(flags))

"""
    EngineOptions

Defines low-level options for a TEMPO run.

Fields
- `tempo_version::AbstractTempoVersion` : Which TEMPO flavor to use
- `flags::String`                       : Additional command-line flags (normalized string)
- `nits::Int`                           : Number of internal iterations (must be ≥ 1)
- `gain::Float64`                       : GAIN parameter controlling convergence damping (must be > 0)

Notes
- `flags` should be a single string; it is normalized (trimmed) and must not contain newlines.
- `gain` is validated to be strictly positive.
"""
struct EngineOptions
    tempo_version::AbstractTempoVersion
    flags::String
    nits::Int
    gain::Float64
    function EngineOptions(tempo_version::AbstractTempoVersion,
                           flags::AbstractString,
                           nits::Integer,
                           gain::Real)
        nits < 1          && error("EngineOptions: nits must be ≥ 1, got $nits")
        gain <= 0         && error("EngineOptions: gain must be > 0, got $(gain)")
        f = _normalize_flags(flags)
        occursin('\n', f) && error("EngineOptions: flags must not contain newlines")
        return new(tempo_version, f, Int(nits), Float64(gain))
    end
end

function EngineOptions(; tempo_version::AbstractTempoVersion,
                         flags::Union{AbstractString,AbstractVector{<:AbstractString}}="",
                         nits::Integer=1, gain::Real=1.0)
    fl = flags isa AbstractVector ? join(flags, " ") : flags
    return EngineOptions(tempo_version, fl, nits, gain)
end

function Base.show(io::IO, ::MIME"text/plain", opts::EngineOptions)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    println(io, pad, "EngineOptions")

    # inline tempo_version (one line, no parentheses)
    print(io, spad, "tempo_version: ")
    show(IOContext(io, :inline => true), MIME"text/plain"(), opts.tempo_version)
    println(io)

    if isempty(opts.flags)
        println(io, spad, "flags:         (none)")
    else
        println(io, spad, "flags:         ", opts.flags)
    end
    println(io, spad, "nits:          ", opts.nits)
    println(io, spad, "gain:          ", opts.gain)
end

# --------------------------------------------------------------------------------------------------------------
# Modifiers
# --------------------------------------------------------------------------------------------------------------

"""
    InputModifiers

Describes transformations applied to the *inputs* of a TEMPO run.

Fields
- `override_params::Vector{TempoParameter}` : parameters to upsert into the `.par`
- `time_start::Union{Nothing,Float64}`      : optional lower MJD bound (inclusive)
- `time_finish::Union{Nothing,Float64}`     : optional upper MJD bound (inclusive)
- `couple_f1_to_ddot::Bool`                 : if true, auto-adjust `F1` when `DDOT` is overridden (applied during materialization)

Notes
- If both `time_start` and `time_finish` are provided, `time_start ≤ time_finish` is enforced.
"""
struct InputModifiers
    override_params::Vector{TempoParameter}
    time_start::Union{Nothing, Float64}
    time_finish::Union{Nothing, Float64}
    couple_f1_to_ddot::Bool
    function InputModifiers(override_params::Vector{TempoParameter},
                            time_start::Union{Nothing, Real},
                            time_finish::Union{Nothing, Real},
                            couple_f1_to_ddot::Bool)
        ts = time_start  === nothing ? nothing : Float64(time_start)
        tf = time_finish === nothing ? nothing : Float64(time_finish)
        if ts !== nothing && tf !== nothing && ts > tf
            error("InputModifiers: time_start ($ts) must be ≤ time_finish ($tf)")
        end
        return new(override_params, ts, tf, couple_f1_to_ddot)
    end
end

function InputModifiers(; override_params::Vector{TempoParameter}=TempoParameter[],
                          time_start::Union{Nothing, Real}=nothing,
                          time_finish::Union{Nothing, Real}=nothing,
                          couple_f1_to_ddot::Bool=false)
    InputModifiers(override_params, time_start, time_finish, couple_f1_to_ddot)
end

function Base.show(io::IO, ::MIME"text/plain", mods::InputModifiers)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    tpad   = repeat(" ", indent + 4)

    println(io, pad,  "InputModifiers")
    println(io, spad, "override_params: ", length(mods.override_params), " parameter(s)")
    for p in mods.override_params
        println(io, tpad, p)
    end
    println(io, spad, "time_start:      ", isnothing(mods.time_start)  ? "nothing" : mods.time_start)
    println(io, spad, "time_finish:     ", isnothing(mods.time_finish) ? "nothing" : mods.time_finish)
    println(io, spad, "couple_f1_to_ddot: ", mods.couple_f1_to_ddot)
end

# --------------------------------------------------------------------------------------------------------------
# Capture options (what we ask TEMPO/TEMPO2 to emit to files/stdout)
# --------------------------------------------------------------------------------------------------------------

"""
    CaptureOptions

Controls what is captured/emitted by the TEMPO engine itself (flags / process I/O).

Fields
- `write_output::Bool`    : capture stdout/stderr from the engine
- `write_residuals::Bool` : for Tempo2 uses `-write_residuals` (per internal iteration);
                            for Tempo uses `-residuals` (final only)
"""
struct CaptureOptions
    write_output::Bool
    write_residuals::Bool
end

CaptureOptions(; write_output::Bool=true,
                 write_residuals::Bool=true) =
    CaptureOptions(write_output, write_residuals)

function Base.show(io::IO, ::MIME"text/plain", c::CaptureOptions)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "CaptureOptions")
    println(io, spad, "write_output:    ", c.write_output)
    println(io, spad, "write_residuals: ", c.write_residuals)
end

# --------------------------------------------------------------------------------------------------------------
# Retention options (what we keep inside Julia results)
# --------------------------------------------------------------------------------------------------------------

"""
    RetentionOptions

Controls what data we keep in the in-memory/result structures after a run.

Fields
- `save_internal_iterations::Bool` : retain all intermediate iteration results
- `save_residuals::Bool`           : retain residuals arrays inside results (can be large)
"""
struct RetentionOptions
    save_internal_iterations::Bool
    save_residuals::Bool
end

RetentionOptions(; save_internal_iterations::Bool=false,
                   save_residuals::Bool=false) =
    RetentionOptions(save_internal_iterations, save_residuals)

function Base.show(io::IO, ::MIME"text/plain", r::RetentionOptions)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    println(io, pad, "RetentionOptions")
    println(io, spad, "save_internal_iterations: ", r.save_internal_iterations)
    println(io, spad, "save_residuals:           ", r.save_residuals)
end

# --------------------------------------------------------------------------------------------------------------
# White-noise analysis options
# --------------------------------------------------------------------------------------------------------------

const _WN_ALLOWED_SCOPES = (:final, :all)

"""
    WhiteNoiseOptions(; enabled=false, scope=:final)

Options for white-noise analysis:
- `enabled` — toggle analysis;
- `scope`   — `:final` (last internal iteration) or `:all` (every internal iteration).
"""
struct WhiteNoiseOptions
    enabled::Bool
    scope::Symbol   # :final | :all
    function WhiteNoiseOptions(enabled::Bool, scope::Symbol)
        scope in _WN_ALLOWED_SCOPES || error("white_noise_scope must be one of $(_WN_ALLOWED_SCOPES); got :$scope")
        return new(enabled, scope)
    end
end

WhiteNoiseOptions(; enabled::Bool=false, scope::Symbol=:final) =
    WhiteNoiseOptions(enabled, scope)

function Base.show(io::IO, ::MIME"text/plain", wn::WhiteNoiseOptions)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    println(io, pad, "WhiteNoiseOptions")
    println(io, spad, "enabled: ", wn.enabled)
    println(io, spad, "scope:   ", wn.scope)
end

# --------------------------------------------------------------------------------------------------------------
# Workspace Options
# --------------------------------------------------------------------------------------------------------------

"""
    WorkspaceOptions

Workspace/runtime controls for a single TEMPO run. The options are grouped by intent:

Lifecycle & isolation
- `work_mode::Symbol`   : `:inplace` | `:jobdir`
    * `:inplace` — run directly inside `work_dir` (job root = `work_dir`).
    * `:jobdir`  — create/use a subdirectory of `work_dir` as an isolated **job root**.
- `job_name::Union{Nothing,String}` : subdirectory name when `work_mode=:jobdir`.
    * If `nothing`, an auto name like `job-YYYYmmdd-HHMMSS-<rand>` is generated.
- `overwrite::Symbol`   : `:error` | `:reuse` | `:unique` | `:clean`
    * Behavior when the target job directory already exists (`:unique` appends `-001`, `-002`, ...).
    * `:error`  — fail if exists; `:reuse` — use as-is; `:clean` — purge contents then use.
    * `:overwrite` is accepted for backward compatibility and treated as `:clean` (a warning is emitted).

Layout of files (inside the job root)
- `layout::Symbol`      : `:flat` | `:split`
    * `:flat`  — files live directly in the job root.
    * `:split` — create `input/`, `output/`, `tmp/` under the job root.
- `temp_dir::Union{Nothing,String}`
    * If `nothing`, the execution `cwd` is the job root for `:flat`, or `tmp/` for `:split`.
    * If non-`nothing`, use this subdirectory as the execution `cwd` (created if missing).

**I/O mirroring**
- `io_mirror::Union{Symbol,Int,Tuple{Symbol,Int}} = :none`
    * Controls how the `temp_dir` (or job root) directory structure is mirrored for I/O.
    * `:none` (flat), `:full` (mirror full `temp_dir`), `Int` (first N segments), `(:depth_minus, k)` (first depth-k segments).
    * Ignored when `layout = :flat`.

Inputs capture into job root
- `link_tim::Bool`     : symlink the `.tim` file instead of copying (falls back to copy if symlink fails).
- `snapshot_par::Bool` : copy the input `.par` for reproducibility.

Cleanup & retention
- `cleanup_before_run::Bool`   : remove engine artifacts in the run directory before launch.
- `keep_tmp_on_success::Bool`  : keep `tmp/` after a successful run (useful for debugging).
- `keep_tmp_on_error::Bool`    : keep `tmp/` if the run failed.

Runtime & manifest
- `timeout_s::Union{Nothing,Int}` : kill the run after this many seconds (`nothing` = unlimited).
- `write_manifest::Bool`          : write a small manifest file in the job root.
- `manifest_style::Symbol`        : `:json` | `:toml`.
"""
struct WorkspaceOptions
    # Lifecycle & isolation
    work_mode::Symbol
    job_name::Union{Nothing,String}
    overwrite::Symbol

    # Layout of files
    layout::Symbol
    temp_dir::Union{Nothing,String}
    io_mirror::Union{Symbol,Int,Tuple{Symbol,Int}}

    # Inputs capture into job root
    link_tim::Bool
    snapshot_par::Bool

    # Cleanup & retention
    cleanup_before_run::Bool
    keep_tmp_on_success::Bool
    keep_tmp_on_error::Bool

    # Runtime & manifest
    timeout_s::Union{Nothing,Int}
    write_manifest::Bool
    manifest_style::Symbol
end

const _WORK_MODES     = (:inplace, :jobdir)
const _LAYOUTS        = (:flat, :split)
const _OVERWRITE_POL  = (:error, :reuse, :unique, :clean, :overwrite)  # :overwrite is a legacy alias for :clean
const _MANIFEST_STYLE = (:json, :toml)



# Validator for io_mirror
function _validate_io_mirror(x)
    if x isa Symbol
        x in (:none, :full) || error("WorkspaceOptions: io_mirror Symbol must be :none or :full, got :$x")
    elseif x isa Int
        x >= 0 || error("WorkspaceOptions: io_mirror Int must be >= 0, got $x")
    elseif x isa Tuple{Symbol,Int}
        x[1] == :depth_minus || error("WorkspaceOptions: io_mirror tuple must have first element :depth_minus, got $(x[1])")
        x[2] >= 0 || error("WorkspaceOptions: io_mirror tuple second element must be >= 0, got $(x[2])")
    else
        error("WorkspaceOptions: io_mirror must be Symbol (:none/:full), Int (>=0), or (:depth_minus, Int>=0); got $(x)")
    end
    return x
end

WorkspaceOptions(; 
    # Lifecycle & isolation
    work_mode::Symbol = :inplace,
    job_name::Union{Nothing,AbstractString} = nothing,
    overwrite::Symbol = :error,

    # Layout of files
    layout::Symbol = :flat,
    temp_dir::Union{Nothing,AbstractString} = nothing,
    io_mirror::Union{Symbol,Int,Tuple{Symbol,Int}} = :none,

    # Inputs capture into job root
    link_tim::Bool = false,
    snapshot_par::Bool = true,

    # Cleanup & retention
    cleanup_before_run::Bool = true,
    keep_tmp_on_success::Bool = false,
    keep_tmp_on_error::Bool = true,

    # Runtime & manifest
    timeout_s::Union{Nothing,Integer} = nothing,
    write_manifest::Bool = false,
    manifest_style::Symbol = :json,
) = begin
    work_mode in _WORK_MODES           || error("WorkspaceOptions: work_mode must be one of $(_WORK_MODES)")
    layout    in _LAYOUTS              || error("WorkspaceOptions: layout must be one of $(_LAYOUTS)")
    manifest_style in _MANIFEST_STYLE  || error("WorkspaceOptions: manifest_style must be one of $(_MANIFEST_STYLE)")
    overwrite in _OVERWRITE_POL        || error("WorkspaceOptions: overwrite must be one of $(_OVERWRITE_POL)")
    _validate_io_mirror(io_mirror)

    WorkspaceOptions(
        # Lifecycle & isolation
        work_mode,
        job_name === nothing ? nothing : String(job_name),
        overwrite,
        # Layout of files
        layout,
        temp_dir === nothing ? nothing : String(temp_dir),
        io_mirror,
        # Inputs capture into job root
        link_tim,
        snapshot_par,
        # Cleanup & retention
        cleanup_before_run,
        keep_tmp_on_success,
        keep_tmp_on_error,
        # Runtime & manifest
        timeout_s === nothing ? nothing : Int(timeout_s),
        write_manifest,
        manifest_style,
    )
end

function Base.show(io::IO, ::MIME"text/plain", p::WorkspaceOptions)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "WorkspaceOptions")

    # Lifecycle & isolation
    println(io, spad, "work_mode:            ", p.work_mode)
    println(io, spad, "job_name:             ", p.job_name)
    println(io, spad, "overwrite:            ", p.overwrite)

    # Layout & temp
    println(io, spad, "layout:               ", p.layout)
    println(io, spad, "temp_dir:             ", p.temp_dir)
    println(io, spad, "io_mirror:            ", p.io_mirror)

    # Inputs capture
    println(io, spad, "link_tim:             ", p.link_tim)
    println(io, spad, "snapshot_par:         ", p.snapshot_par)

    # Cleanup & retention
    println(io, spad, "cleanup_before_run:   ", p.cleanup_before_run)
    println(io, spad, "keep_tmp_on_success:  ", p.keep_tmp_on_success)
    println(io, spad, "keep_tmp_on_error:    ", p.keep_tmp_on_error)

    # Runtime & manifest
    println(io, spad, "timeout_s:            ", p.timeout_s)
    println(io, spad, "write_manifest:       ", p.write_manifest)
    println(io, spad, "manifest_style:       ", p.manifest_style)
end

# --------------------------------------------------------------------------------------------------------------
# Logging Options
# --------------------------------------------------------------------------------------------------------------

# Allowed verbosity levels
const _VERBOSITY_LEVELS = (:silent, :warn, :info, :debug)

# Normalize verbosity from Int or Symbol into a Symbol from _VERBOSITY_LEVELS
# Mapping (for convenience when users pass integers):
#   0=>:silent, 1=>:warn, 2=>:info, 3=>:debug (values <0 -> :silent, >3 -> :debug)
@inline function _normalize_verbosity(v::Integer)
    v <= 0 && return :silent
    v == 1 && return :warn
    v == 2 && return :info
    return :debug
end
@inline function _normalize_verbosity(v::Symbol)
    v in _VERBOSITY_LEVELS || error("verbosity must be one of $(_VERBOSITY_LEVELS); got :$v")
    return v
end

"""
    LoggingOptions

Lightweight logging verbosity controls for the runner.
- `verbosity` — one of `:silent | :warn | :info | :debug` (also accepts `0..3` in the constructor)
- `with_timestamps` — include timestamps in log lines
"""
struct LoggingOptions
    verbosity::Symbol        # :silent | :warn | :info | :debug
    with_timestamps::Bool
    function LoggingOptions(verbosity::Union{Integer,Symbol}, with_timestamps::Bool)
        verbosity = _normalize_verbosity(verbosity)
        return new(verbosity, with_timestamps)
    end
end

LoggingOptions(; verbosity::Union{Integer,Symbol}=:info, with_timestamps::Bool=true) =
    LoggingOptions(verbosity, with_timestamps)

function Base.show(io::IO, ::MIME"text/plain", l::LoggingOptions)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    println(io, pad, "LoggingOptions")
    println(io, spad, "verbosity:       ", l.verbosity)
    println(io, spad, "with_timestamps: ", l.with_timestamps)
end

# Lightweight predicates for convenience across the runner/diagnostics
"""
    is_silent(l::LoggingOptions) -> Bool
"""
is_silent(l::LoggingOptions) = l.verbosity === :silent
"""
    is_warn(l::LoggingOptions) -> Bool
"""
is_warn(l::LoggingOptions)   = l.verbosity === :warn
"""
    is_info(l::LoggingOptions) -> Bool

Treats `:debug` as informational as well (superset).
"""
is_info(l::LoggingOptions)   = (l.verbosity === :info) || (l.verbosity === :debug)
"""
    is_debug(l::LoggingOptions) -> Bool
"""
is_debug(l::LoggingOptions)  = l.verbosity === :debug

# --------------------------------------------------------------------------------------------------------------
# Basic Tempo settings
# --------------------------------------------------------------------------------------------------------------

"""
    TempoRunSettings

Top-level settings for a single TEMPO run. Composed from smaller structures:

- `paths::RunPaths`
- `engine::EngineOptions`
- `modifiers::InputModifiers`
- `capture::CaptureOptions`
- `retention::RetentionOptions`
- `analysis::WhiteNoiseOptions`
- `workspace::WorkspaceOptions`
- `logging::LoggingOptions`
"""
struct TempoRunSettings <: AbstractTempoSettings
    paths::RunPaths
    engine::EngineOptions
    modifiers::InputModifiers
    capture::CaptureOptions
    retention::RetentionOptions
    analysis::WhiteNoiseOptions
    workspace::WorkspaceOptions
    logging::LoggingOptions
end

function Base.show(io::IO, ::MIME"text/plain", s::TempoRunSettings)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    iop    = IOContext(io, :indent => indent + 4)

    println(io, pad, "TempoRunSettings")
    println(io, spad, "paths:")
    show(iop, MIME"text/plain"(), s.paths)
    println(io, spad, "engine:")
    show(iop, MIME"text/plain"(), s.engine)
    println(io, spad, "modifiers:")
    show(iop, MIME"text/plain"(), s.modifiers)
    println(io, spad, "capture:")
    show(iop, MIME"text/plain"(), s.capture)
    println(io, spad, "retention:")
    show(iop, MIME"text/plain"(), s.retention)
    println(io, spad, "analysis:")
    show(iop, MIME"text/plain"(), s.analysis)
    println(io, spad, "workspace:")
    show(iop, MIME"text/plain"(), s.workspace)
    println(io, spad, "logging:")
    show(iop, MIME"text/plain"(), s.logging)
end

# --------------------------------------------------------------------------------------------------------------
# Convenience constructor
# --------------------------------------------------------------------------------------------------------------

"""
    TempoRunSettings(; kwargs...) -> TempoRunSettings

Keyword-style construction with sensible defaults.
Raises an error if `nits < 1`.

Keyword arguments
- Paths: `work_dir`, `par_input`, `par_output`, `tim_input`
- Execution: `tempo_version`, `flags`, `nits`, `gain`
- Modifiers: `override_params`, `time_start`, `time_finish`, `couple_f1_to_ddot`
- Capture: `write_output`, `write_residuals`
- Retention: `save_internal_iterations`, `save_residuals`
- White-noise analysis (convenience): `white_noise_enabled`, `white_noise_scope` (`:final` | `:all`)
  or pass `analysis::WhiteNoiseOptions`
- Workspace (convenience): `work_mode`, `job_name`, `overwrite` (`:error|:reuse|:unique|:clean`), `layout`, `temp_dir`,
  `link_tim`, `snapshot_par`, `cleanup_before_run`, `keep_tmp_on_success`, `keep_tmp_on_error`,
  `timeout_s`, `write_manifest`, `manifest_style`
  or pass `workspace::WorkspaceOptions`
- Logging (convenience): `verbosity`, `with_timestamps`
  or pass `logging::LoggingOptions`
"""
function TempoRunSettings(;
    # Paths
    work_dir::AbstractString,
    par_input::AbstractString,
    par_output::AbstractString = default_par_output(par_input),
    tim_input::AbstractString,
    paths::Union{Nothing,RunPaths} = nothing,

    # Engine
    tempo_version::AbstractTempoVersion,
    flags::AbstractString = "",
    nits::Integer = 1,
    gain::Real = 1.0,
    engine::Union{Nothing,EngineOptions} = nothing,

    # Modifiers
    override_params::Vector{TempoParameter} = TempoParameter[],
    time_start::Union{Nothing, Real} = nothing,
    time_finish::Union{Nothing, Real} = nothing,
    couple_f1_to_ddot::Bool = false,
    modifiers::Union{Nothing,InputModifiers} = nothing,

    # Capture
    write_output::Bool = true,
    write_residuals::Bool = true,
    capture::Union{Nothing,CaptureOptions} = nothing,

    # Retention
    save_internal_iterations::Bool = false,
    save_residuals::Bool = false,
    retention::Union{Nothing,RetentionOptions} = nothing,

    # White-noise analysis (either convenience keys or full struct)
    white_noise_enabled::Bool = false,
    white_noise_scope::Symbol = :final,
    analysis::Union{Nothing,WhiteNoiseOptions} = nothing,

    # Workspace (either convenience keys or full struct)
    work_mode::Symbol = :inplace,
    job_name::Union{Nothing,AbstractString} = nothing,
    overwrite::Symbol = :error,
    layout::Symbol = :flat,
    temp_dir::Union{Nothing,AbstractString} = nothing,
    io_mirror::Union{Symbol,Int,Tuple{Symbol,Int}} = :none,
    link_tim::Bool = false,
    snapshot_par::Bool = true,
    cleanup_before_run::Bool = true,
    keep_tmp_on_success::Bool = false,
    keep_tmp_on_error::Bool = true,
    timeout_s::Union{Nothing,Integer} = nothing,
    write_manifest::Bool = false,
    manifest_style::Symbol = :json,
    workspace::Union{Nothing,WorkspaceOptions} = nothing,

    # Logging (either convenience keys or full struct)
    verbosity::Union{Integer,Symbol} = :info,
    with_timestamps::Bool = true,
    logging::Union{Nothing,LoggingOptions} = nothing,
)
    paths_new = paths !== nothing ? paths : RunPaths(
        work_dir, 
        par_input, 
        par_output, 
        tim_input
    )

    engine_new = engine !== nothing ? engine : EngineOptions(
        tempo_version,
        flags,
        nits,
        gain
    )

    modifiers_new = modifiers !== nothing ? modifiers : InputModifiers(
        override_params,
        time_start,
        time_finish,
        couple_f1_to_ddot,
    )

    capture_new = capture !== nothing ? capture : CaptureOptions(
        write_output,
        write_residuals
    )

    retention_new = retention !== nothing ? retention : RetentionOptions(
        save_internal_iterations,
        save_residuals
    )

    analysis_new = analysis !== nothing ? analysis : WhiteNoiseOptions(
        white_noise_enabled,
        white_noise_scope
    )

    workspace_new = workspace !== nothing ? workspace : WorkspaceOptions(
        work_mode           = work_mode,
        job_name            = job_name,
        overwrite           = overwrite,
        layout              = layout,
        temp_dir            = temp_dir,
        io_mirror           = io_mirror,
        link_tim            = link_tim,
        snapshot_par        = snapshot_par,
        cleanup_before_run  = cleanup_before_run,
        keep_tmp_on_success = keep_tmp_on_success,
        keep_tmp_on_error   = keep_tmp_on_error,
        timeout_s           = timeout_s,
        write_manifest      = write_manifest,
        manifest_style      = manifest_style,
    )

    logging_new = logging !== nothing ? logging : LoggingOptions(
        verbosity, 
        with_timestamps
    )

    return TempoRunSettings(paths_new, engine_new, modifiers_new, capture_new, retention_new, analysis_new, workspace_new, logging_new)
end

# --------------------------------------------------------------------------------------------------------------
# Copy with overrides (extended)
# --------------------------------------------------------------------------------------------------------------

"""
    copy_with(s::TempoRunSettings; kwargs...) -> TempoRunSettings

Create a modified copy of `s` with keyword overrides.
You can override either whole sub-structs (`analysis`, `workspace`, `logging`)
or individual convenience keys (e.g. `white_noise_enabled`, `timeout_s`, `verbosity`, ...).

Additional override-params controls:
- `override_params_clear::Bool=false`: start from an empty override list.
- `override_params_delete`: names (Symbol/String or a collection of them) to remove from overrides.
- `override_params_upsert::Vector{TempoParameter}=TP[]`: add/replace overrides by name.
If `override_params` is provided, delete/upsert are applied on top of it.
"""
function copy_with(s::TempoRunSettings; kwargs...)
    # Paths
    paths_kw = get(kwargs, :paths, nothing)
    paths = paths_kw !== nothing ? paths_kw : RunPaths(
        get(kwargs, :work_dir, s.paths.work_dir), 
        get(kwargs, :par_input, s.paths.par_input), 
        get(kwargs, :par_output, s.paths.par_output), 
        get(kwargs, :tim_input, s.paths.tim_input)  
    )

    # Engine
    engine_kw = get(kwargs, :engine, nothing)
    engine = engine_kw !== nothing ? engine_kw : EngineOptions(
        get(kwargs, :tempo_version, s.engine.tempo_version),
        get(kwargs, :flags, s.engine.flags),
        get(kwargs, :nits, s.engine.nits),
        get(kwargs, :gain, s.engine.gain)
    )

    # Modifiers — base vector and extended controls
    base_overrides  = get(kwargs, :override_params, s.modifiers.override_params)
    clear_overrides = get(kwargs, :override_params_clear, false)
    del_names_any   = get(kwargs, :override_params_delete, nothing)
    upsert_vec      = get(kwargs, :override_params_upsert, TempoParameter[])::Vector{TempoParameter}

    overrides = clear_overrides ? TempoParameter[] : copy(base_overrides)

    if del_names_any !== nothing
        del_list = del_names_any isa AbstractVector ? del_names_any : (del_names_any,)
        del_syms = Symbol.(String.(del_list))
        overrides = without_params(overrides, del_syms)
    end

    if !isempty(upsert_vec)
        overrides = with_upserted_params(overrides, upsert_vec)
    end

    time_start      = get(kwargs, :time_start, s.modifiers.time_start)
    time_finish     = get(kwargs, :time_finish, s.modifiers.time_finish)
    couple_f1_to_ddot = get(kwargs, :couple_f1_to_ddot, s.modifiers.couple_f1_to_ddot)

    # Rebuild modifiers struct
    modifiers_kw = get(kwargs, :modifiers, nothing)
    modifiers = modifiers_kw !== nothing ? modifiers_kw : InputModifiers(overrides, time_start, time_finish, couple_f1_to_ddot)

    # Capture
    capture_kw = get(kwargs, :capture, nothing)
    capture = capture_kw !== nothing ? capture_kw : CaptureOptions(
        get(kwargs, :write_output, s.capture.write_output),
        get(kwargs, :write_residuals, s.capture.write_residuals)
    )

    # Retention
    retention_kw = get(kwargs, :retention, nothing)
    retention = retention_kw !== nothing ? retention_kw : RetentionOptions(
        get(kwargs, :save_internal_iterations, s.retention.save_internal_iterations),
        get(kwargs, :save_residuals, s.retention.save_residuals)
    )

    # Analysis (prefer full struct if provided)
    analysis_kw = get(kwargs, :analysis, nothing)
    analysis    = analysis_kw !== nothing ? analysis_kw : WhiteNoiseOptions(
        get(kwargs, :white_noise_enabled, s.analysis.enabled),
        get(kwargs, :white_noise_scope, s.analysis.scope)
    )

    # Workspace (prefer full struct if provided)
    workspace_kw = get(kwargs, :workspace, nothing)
    workspace    = workspace_kw !== nothing ? workspace_kw : WorkspaceOptions(
        timeout_s           = get(kwargs, :timeout_s, s.workspace.timeout_s),
        cleanup_before_run  = get(kwargs, :cleanup_before_run, s.workspace.cleanup_before_run),
        temp_dir            = get(kwargs, :temp_dir, s.workspace.temp_dir),
        work_mode           = get(kwargs, :work_mode, s.workspace.work_mode),
        job_name            = get(kwargs, :job_name, s.workspace.job_name),
        layout              = get(kwargs, :layout, s.workspace.layout),
        overwrite           = get(kwargs, :overwrite, s.workspace.overwrite),
        io_mirror           = get(kwargs, :io_mirror, s.workspace.io_mirror),
        link_tim            = get(kwargs, :link_tim, s.workspace.link_tim),
        snapshot_par        = get(kwargs, :snapshot_par, s.workspace.snapshot_par),
        keep_tmp_on_success = get(kwargs, :keep_tmp_on_success, s.workspace.keep_tmp_on_success),
        keep_tmp_on_error   = get(kwargs, :keep_tmp_on_error, s.workspace.keep_tmp_on_error),
        write_manifest      = get(kwargs, :write_manifest, s.workspace.write_manifest),
        manifest_style      = get(kwargs, :manifest_style, s.workspace.manifest_style),
    )

    # Logging (prefer full struct if provided)
    logging_kw = get(kwargs, :logging, nothing)
    logging    = logging_kw !== nothing ? logging_kw : LoggingOptions(
        get(kwargs, :verbosity, s.logging.verbosity),
        get(kwargs, :with_timestamps, s.logging.with_timestamps),
    )

    return TempoRunSettings(paths, engine, modifiers, capture, retention, analysis, workspace, logging)
end


# --------------------------------------------------------------------------------------------------------------
# Optional helpers
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# Path helpers (absolute input paths resolved against work_dir)
# --------------------------------------------------------------------------------------------------------------
"""
    par_input_path(s::TempoRunSettings) -> String

Absolute path to the input `.par` file, resolving a relative name against `s.paths.work_dir`.
"""
par_input_path(s::TempoRunSettings) = joinpath(s.paths.work_dir, s.paths.par_input)

"""
    tim_input_path(s::TempoRunSettings) -> String

Absolute path to the input `.tim` file, resolving a relative name against `s.paths.work_dir`.
"""
tim_input_path(s::TempoRunSettings) = joinpath(s.paths.work_dir, s.paths.tim_input)

"""
    validate(s::TempoRunSettings) -> Bool

Alias for [`validate_inputs_exist`] — quick pre-materialization check that input `.par`/`.tim` exist
relative to `work_dir` and that `par_output` is a bare file name. Does not create directories.
"""
validate(s::TempoRunSettings) = validate_inputs_exist(s)

"""
    validate_inputs_exist(s::TempoRunSettings) -> Bool

Validate inputs before materialization:
- checks that `par_input` and `tim_input` exist relative to `work_dir`;
- ensures `par_output` is a file name (no directories).
Does not create directories. Returns true if valid, otherwise throws.
"""
function validate_inputs_exist(s::TempoRunSettings)
    par_in = joinpath(s.paths.work_dir, s.paths.par_input)
    tim_in = joinpath(s.paths.work_dir, s.paths.tim_input)
    basename(s.paths.par_output) == s.paths.par_output ||
        error("par_output must be a file name (no directories): $(s.paths.par_output)")
    isfile(par_in) || error("Input par file not found: $par_in")
    isfile(tim_in) || error("TIM file not found: $tim_in")
    return true
end
