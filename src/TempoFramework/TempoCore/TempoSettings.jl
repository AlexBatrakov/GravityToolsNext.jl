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
# Files
# --------------------------------------------------------------------------------------------------------------

"""
    TempoRunFiles

Holds paths for a TEMPO run.

Fields
- `work_dir::String`        : Absolute working directory for the run
- `par_file_input::String`  : Input `.par` **file name or relative path** (resolved against `work_dir`)
- `par_file_output::String` : Output `.par` **file name** to be written in `work_dir`
- `tim_file::String`        : `.tim` **file name or relative path** (resolved against `work_dir`)
"""
struct TempoRunFiles
    work_dir::String
    par_file_input::String
    par_file_output::String
    tim_file::String
end

function TempoRunFiles(work_dir::AbstractString,
                       par_file_input::AbstractString,
                       par_file_output::AbstractString,
                       tim_file::AbstractString)
    return TempoRunFiles(String(work_dir),
                         String(par_file_input),
                         String(par_file_output),
                         String(tim_file))
end

function Base.show(io::IO, ::MIME"text/plain", f::TempoRunFiles)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "TempoRunFiles")
    println(io, spad, "work_dir:        ", f.work_dir)
    println(io, spad, "par_file_input:  ", f.par_file_input)
    println(io, spad, "par_file_output: ", f.par_file_output)
    println(io, spad, "tim_file:        ", f.tim_file)
end
# Case-insensitive `.par` handling
_par_stripped(name::AbstractString) = replace(name, r"\.[Pp][Aa][Rr]$" => "")

"""
    default_par_file_output(par_file_input::AbstractString) -> String

Generate a default output filename by replacing the extension with `_out.par`.
Throws an error if the input does not look like a `.par` file.
"""
function default_par_file_output(par_file_input::AbstractString)
    if !occursin(r"\.[Pp][Aa][Rr]$", par_file_input)
        error("Expected input par file to end with .par (case-insensitive), got: $par_file_input")
    end
    return string(_par_stripped(par_file_input), "_out.par")
end

# --------------------------------------------------------------------------------------------------------------
# Execution options
# --------------------------------------------------------------------------------------------------------------

"""
    TempoExecutionOptions

Defines low-level options for a TEMPO run.

Fields
- `tempo_version::AbstractTempoVersion` : Which TEMPO flavor to use
- `flags::String`                       : Additional command-line flags for TEMPO
- `nits::Int`                           : Number of internal iterations (must be ≥ 1)
- `gain::Float64`                       : GAIN parameter controlling convergence damping
"""
struct TempoExecutionOptions
    tempo_version::AbstractTempoVersion
    flags::String
    nits::Int
    gain::Float64
end

function TempoExecutionOptions(tempo_version::AbstractTempoVersion,
                               flags::AbstractString,
                               nits::Integer,
                               gain::Real)
    nits < 1 && error("TempoExecutionOptions: nits must be ≥ 1, got $nits")
    return TempoExecutionOptions(tempo_version, String(flags), Int(nits), Float64(gain))
end

function Base.show(io::IO, ::MIME"text/plain", opts::TempoExecutionOptions)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    println(io, pad, "TempoExecutionOptions")

    # inline tempo_version (one line, no parentheses)
    print(io, pad, "  tempo_version: ")
    show(IOContext(io, :inline => true), MIME"text/plain"(), opts.tempo_version)
    println(io)

    println(io, pad, "  flags:         ", opts.flags)
    println(io, pad, "  nits:          ", opts.nits)
    println(io, pad, "  gain:          ", opts.gain)
end

# --------------------------------------------------------------------------------------------------------------
# Modifiers
# --------------------------------------------------------------------------------------------------------------

"""
    TempoRunModifiers

Describes modifications applied to the run.

Fields
- `override_params::Vector{TempoParameter}` : Parameters to override/inject into the `.par`
- `time_start::Union{Nothing,Float64}`             : Optional lower bound on TOAs (MJD)
- `time_finish::Union{Nothing,Float64}`            : Optional upper bound on TOAs (MJD)
"""
struct TempoRunModifiers
    override_params::Vector{TempoParameter}
    time_start::Union{Nothing, Float64}
    time_finish::Union{Nothing, Float64}
end

function TempoRunModifiers(override_params::Vector{TempoParameter}=TempoParameter[];
                           time_start::Union{Nothing, Real}=nothing,
                           time_finish::Union{Nothing, Real}=nothing)
    ts = time_start === nothing ? nothing : Float64(time_start)
    tf = time_finish === nothing ? nothing : Float64(time_finish)
    if ts !== nothing && tf !== nothing && ts > tf
        error("TempoRunModifiers: time_start ($ts) must be ≤ time_finish ($tf)")
    end
    return TempoRunModifiers(override_params, ts, tf)
end

function Base.show(io::IO, ::MIME"text/plain", mods::TempoRunModifiers)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    println(io, pad, "TempoRunModifiers")
    println(io, pad, "  override_params: ", length(mods.override_params), " parameter(s)")
    for p in mods.override_params
        println(io, pad, "    ", p)
    end
    println(io, pad, "  time_start:      ", isnothing(mods.time_start) ? "nothing" : mods.time_start)
    println(io, pad, "  time_finish:     ", isnothing(mods.time_finish) ? "nothing" : mods.time_finish)
end

# --------------------------------------------------------------------------------------------------------------
# Behavior
# --------------------------------------------------------------------------------------------------------------

"""
    TempoRunBehavior

Controls what is captured/saved during the run.

Fields
- `write_output::Bool`            : Capture stdout/stderr from TEMPO
- `write_residuals::Bool`         : For Tempo2 uses `-write_residuals` (per internal iteration);
                                    for Tempo uses `-residuals` (final only)
- `save_internal_iterations::Bool`: Retain all intermediate iteration results
- `save_residuals::Bool`          : Retain residuals arrays inside results (can be large)
"""
struct TempoRunBehavior
    write_output::Bool
    write_residuals::Bool
    save_internal_iterations::Bool
    save_residuals::Bool
end

TempoRunBehavior(; write_output::Bool=true,
                   write_residuals::Bool=true,
                   save_internal_iterations::Bool=false,
                   save_residuals::Bool=false) =
    TempoRunBehavior(write_output, write_residuals, save_internal_iterations, save_residuals)

function Base.show(io::IO, ::MIME"text/plain", b::TempoRunBehavior)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    println(io, pad, "TempoRunBehavior")
    println(io, pad, "  write_output:             ", b.write_output)
    println(io, pad, "  write_residuals:          ", b.write_residuals)
    println(io, pad, "  save_internal_iterations: ", b.save_internal_iterations)
    println(io, pad, "  save_residuals:           ", b.save_residuals)
end

# --------------------------------------------------------------------------------------------------------------
# White-noise analysis options
# --------------------------------------------------------------------------------------------------------------

"""
    WhiteNoiseAnalysisOptions(; enabled=false, scope=:final)

Options for white-noise analysis:
- `enabled` — toggle analysis;
- `scope`   — `:final` (last internal iteration) or `:all` (every internal iteration).
"""
struct WhiteNoiseAnalysisOptions
    enabled::Bool
    scope::Symbol   # :final | :all
end

const _WN_ALLOWED_SCOPES = (:final, :all)

@inline function _validate_wn_scope(scope::Symbol)
    scope in _WN_ALLOWED_SCOPES || error("white_noise_scope must be one of $(_WN_ALLOWED_SCOPES); got :$scope")
    return scope
end

WhiteNoiseAnalysisOptions(; enabled::Bool=false, scope::Symbol=:final) =
    WhiteNoiseAnalysisOptions(enabled, _validate_wn_scope(scope))

function Base.show(io::IO, ::MIME"text/plain", wn::WhiteNoiseAnalysisOptions)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    println(io, pad, "WhiteNoiseAnalysisOptions")
    println(io, pad, "  enabled: ", wn.enabled)
    println(io, pad, "  scope:   ", wn.scope)
end

# --------------------------------------------------------------------------------------------------------------
# Process / Logging
# --------------------------------------------------------------------------------------------------------------

"""
    ProcessOptions

Minimal process/runtime controls for a single run.

Fields
- `timeout_s`          : kill the run after this many seconds (or `nothing`)
- `cleanup_before_run` : delete/clean run directory artifacts before launch
- `temp_dir`           : if non-`nothing`, run TEMPO in this temp directory
"""
struct ProcessOptions
    timeout_s::Union{Nothing,Int}
    cleanup_before_run::Bool
    temp_dir::Union{Nothing,String}
end

ProcessOptions(; timeout_s=nothing, cleanup_before_run::Bool=true, temp_dir::Union{Nothing,AbstractString}=nothing) =
    ProcessOptions(timeout_s === nothing ? nothing : Int(timeout_s),
                   cleanup_before_run,
                   temp_dir === nothing ? nothing : String(temp_dir))

function Base.show(io::IO, ::MIME"text/plain", p::ProcessOptions)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    println(io, pad, "ProcessOptions")
    println(io, pad, "  timeout_s:          ", p.timeout_s)
    println(io, pad, "  cleanup_before_run: ", p.cleanup_before_run)
    println(io, pad, "  temp_dir:           ", p.temp_dir)
end

"""
    LoggingOptions

Lightweight logging verbosity controls for the runner.
"""
struct LoggingOptions
    verbosity::Int           # 0=silent, 1=info, 2=debug, 3=trace
    with_timestamps::Bool
end

LoggingOptions(; verbosity::Int=1, with_timestamps::Bool=true) =
    LoggingOptions(verbosity, with_timestamps)

function Base.show(io::IO, ::MIME"text/plain", l::LoggingOptions)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    println(io, pad, "LoggingOptions")
    println(io, pad, "  verbosity:       ", l.verbosity)
    println(io, pad, "  with_timestamps: ", l.with_timestamps)
end

# --------------------------------------------------------------------------------------------------------------
# Basic Tempo settings
# --------------------------------------------------------------------------------------------------------------

"""
    BasicTempoSettings

Top-level settings for a single TEMPO run. Composed from smaller structures:

- `files::TempoRunFiles`
- `options::TempoExecutionOptions`
- `modifiers::TempoRunModifiers`
- `behavior::TempoRunBehavior`
- `analysis::WhiteNoiseAnalysisOptions`
- `process::ProcessOptions`
- `logging::LoggingOptions`
"""
struct BasicTempoSettings
    files::TempoRunFiles
    options::TempoExecutionOptions
    modifiers::TempoRunModifiers
    behavior::TempoRunBehavior
    analysis::WhiteNoiseAnalysisOptions
    process::ProcessOptions
    logging::LoggingOptions
end

function Base.show(io::IO, ::MIME"text/plain", s::BasicTempoSettings)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    iop    = IOContext(io, :indent => indent + 4)

    println(io, pad, "BasicTempoSettings")
    println(io, spad, "files:")
    show(iop, MIME"text/plain"(), s.files)
    println(io, spad, "options:")
    show(iop, MIME"text/plain"(), s.options)
    println(io, spad, "modifiers:")
    show(iop, MIME"text/plain"(), s.modifiers)
    println(io, spad, "behavior:")
    show(iop, MIME"text/plain"(), s.behavior)
    println(io, spad, "analysis:")
    show(iop, MIME"text/plain"(), s.analysis)
    println(io, spad, "process:")
    show(iop, MIME"text/plain"(), s.process)
    println(io, spad, "logging:")
    show(iop, MIME"text/plain"(), s.logging)
end

# --------------------------------------------------------------------------------------------------------------
# Convenience constructor
# --------------------------------------------------------------------------------------------------------------

"""
    BasicTempoSettings(; kwargs...) -> BasicTempoSettings

Keyword-style construction with sensible defaults.
Raises an error if `nits < 1`.

Keyword arguments
- Files: `work_dir`, `par_file_input`, `par_file_output`, `tim_file`
- Execution: `tempo_version`, `flags`, `nits`, `gain`
- Modifiers: `override_params`, `time_start`, `time_finish`
- Behavior: `write_output`, `write_residuals`, `save_internal_iterations`, `save_residuals`
- White-noise analysis (convenience): `white_noise_enabled`, `white_noise_scope` (`:final` | `:all`)
  or pass `analysis::WhiteNoiseAnalysisOptions`
- Process (convenience): `timeout_s`, `cleanup_before_run`, `temp_dir`
  or pass `process::ProcessOptions`
- Logging (convenience): `verbosity`, `with_timestamps`
  or pass `logging::LoggingOptions`
"""
function BasicTempoSettings(;
    # Files
    work_dir::AbstractString,
    par_file_input::AbstractString,
    par_file_output::AbstractString = default_par_file_output(par_file_input),
    tim_file::AbstractString,

    # Execution
    tempo_version::AbstractTempoVersion,
    flags::AbstractString = "",
    nits::Integer = 1,
    gain::Real = 1.0,

    # Modifiers
    override_params::Vector{TempoParameter} = TempoParameter[],
    time_start::Union{Nothing, Real} = nothing,
    time_finish::Union{Nothing, Real} = nothing,

    # Behavior
    write_output::Bool = true,
    write_residuals::Bool = true,
    save_internal_iterations::Bool = false,
    save_residuals::Bool = false,

    # White-noise analysis (either convenience keys or full struct)
    white_noise_enabled::Bool = false,
    white_noise_scope::Symbol = :final,
    analysis::Union{Nothing,WhiteNoiseAnalysisOptions} = nothing,

    # Process (either convenience keys or full struct)
    timeout_s::Union{Nothing,Integer} = nothing,
    cleanup_before_run::Bool = true,
    temp_dir::Union{Nothing,AbstractString} = nothing,
    process::Union{Nothing,ProcessOptions} = nothing,

    # Logging (either convenience keys or full struct)
    verbosity::Integer = 1,
    with_timestamps::Bool = true,
    logging::Union{Nothing,LoggingOptions} = nothing,
)
    nits < 1 && error("Number of iterations (nits) must be ≥ 1, got $nits")

    files     = TempoRunFiles(work_dir, par_file_input, par_file_output, tim_file)
    options   = TempoExecutionOptions(tempo_version, flags, nits, Float64(gain))
    modifiers = TempoRunModifiers(override_params, time_start, time_finish)
    behavior  = TempoRunBehavior(write_output, write_residuals, save_internal_iterations, save_residuals)

    analysis_new = analysis !== nothing ?
        analysis :
        WhiteNoiseAnalysisOptions(enabled = white_noise_enabled,
                                  scope   = _validate_wn_scope(white_noise_scope))

    process_new = process !== nothing ?
        process :
        ProcessOptions(timeout_s = timeout_s,
                       cleanup_before_run = cleanup_before_run,
                       temp_dir = temp_dir)

    logging_new = logging !== nothing ?
        logging :
        LoggingOptions(verbosity = Int(verbosity), with_timestamps = with_timestamps)

    return BasicTempoSettings(files, options, modifiers, behavior, analysis_new, process_new, logging_new)
end

# --------------------------------------------------------------------------------------------------------------
# Copy with overrides (extended)
# --------------------------------------------------------------------------------------------------------------

"""
    copy_with(s::BasicTempoSettings; kwargs...) -> BasicTempoSettings

Create a modified copy of `s` with keyword overrides.
You can override either whole sub-structs (`analysis`, `process`, `logging`)
or individual convenience keys (e.g. `white_noise_enabled`, `timeout_s`, `verbosity`, ...).

Additional override-params controls:
- `override_params_clear::Bool=false`: start from an empty override list.
- `override_params_delete`: names (Symbol/String or a collection of them) to remove from overrides.
- `override_params_upsert::Vector{TempoParameter}=TP[]`: add/replace overrides by name.
If `override_params` is provided, delete/upsert are applied on top of it.
"""
function copy_with(s::BasicTempoSettings; kwargs...)
    # Files
    work_dir        = get(kwargs, :work_dir, s.files.work_dir)
    par_file_input  = get(kwargs, :par_file_input, s.files.par_file_input)
    par_file_output = get(kwargs, :par_file_output, s.files.par_file_output)
    tim_file        = get(kwargs, :tim_file, s.files.tim_file)

    # Execution
    tempo_version   = get(kwargs, :tempo_version, s.options.tempo_version)
    flags           = get(kwargs, :flags, s.options.flags)
    nits            = get(kwargs, :nits, s.options.nits)
    gain            = Float64(get(kwargs, :gain, s.options.gain))

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

    # Behavior
    write_output             = get(kwargs, :write_output, s.behavior.write_output)
    write_residuals          = get(kwargs, :write_residuals, s.behavior.write_residuals)
    save_internal_iterations = get(kwargs, :save_internal_iterations, s.behavior.save_internal_iterations)
    save_residuals           = get(kwargs, :save_residuals, s.behavior.save_residuals)

    # Analysis (prefer full struct if provided)
    analysis_kw = get(kwargs, :analysis, nothing)
    analysis    = analysis_kw === nothing ? WhiteNoiseAnalysisOptions(
        enabled = get(kwargs, :white_noise_enabled, s.analysis.enabled),
        scope   = _validate_wn_scope(get(kwargs, :white_noise_scope, s.analysis.scope)),
    ) : analysis_kw

    # Process (prefer full struct if provided)
    process_kw = get(kwargs, :process, nothing)
    process    = process_kw === nothing ? ProcessOptions(
        timeout_s           = get(kwargs, :timeout_s, s.process.timeout_s),
        cleanup_before_run  = get(kwargs, :cleanup_before_run, s.process.cleanup_before_run),
        temp_dir            = get(kwargs, :temp_dir, s.process.temp_dir),
    ) : process_kw

    # Logging (prefer full struct if provided)
    logging_kw = get(kwargs, :logging, nothing)
    logging    = logging_kw === nothing ? LoggingOptions(
        verbosity       = Int(get(kwargs, :verbosity, s.logging.verbosity)),
        with_timestamps = get(kwargs, :with_timestamps, s.logging.with_timestamps),
    ) : logging_kw

    files     = TempoRunFiles(work_dir, par_file_input, par_file_output, tim_file)
    options   = TempoExecutionOptions(tempo_version, flags, nits, gain)
    modifiers = TempoRunModifiers(overrides, time_start, time_finish)
    behavior  = TempoRunBehavior(write_output, write_residuals, save_internal_iterations, save_residuals)

    return BasicTempoSettings(files, options, modifiers, behavior, analysis, process, logging)
end

# --------------------------------------------------------------------------------------------------------------
# Optional helpers
# --------------------------------------------------------------------------------------------------------------

"""
    validate(s::BasicTempoSettings) -> Bool

Validate settings: nits ≥ 1; input `.par` and `.tim` exist.
Resolves relative paths against `work_dir`.
Does not create directories. Returns true if valid, otherwise throws.
"""
function validate(s::BasicTempoSettings)
    s.options.nits >= 1 || error("nits must be ≥ 1 (got $(s.options.nits))")

    par_in = isabspath(s.files.par_file_input) ? s.files.par_file_input : joinpath(s.files.work_dir, s.files.par_file_input)
    tim_in = isabspath(s.files.tim_file)       ? s.files.tim_file       : joinpath(s.files.work_dir, s.files.tim_file)

    isfile(par_in) || error("Input par file not found: $par_in")
    isfile(tim_in) || error("Tim file not found: $tim_in")
    return true
end

"""
    ensure_work_dir!(s::BasicTempoSettings) -> Nothing

Create `s.files.work_dir` if it does not exist.
"""
function ensure_work_dir!(s::BasicTempoSettings)
    isdir(s.files.work_dir) || mkpath(s.files.work_dir)
    nothing
end