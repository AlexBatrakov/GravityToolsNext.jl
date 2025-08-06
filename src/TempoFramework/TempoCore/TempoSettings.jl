#--------------------------------------------------------------------------------------------------------------
"""
Abstract base type for all TEMPO settings structures.
Used for dispatching and type hierarchies.
"""
abstract type AbstractTempoSettings end

#------------------------- FILES -------------------------

"""
Holds paths to files involved in a TEMPO run.

# Fields
- `work_dir::String`: Directory where the run will be executed.
- `par_file_input::String`: Path to the input `.par` file (initial model).
- `par_file_output::String`: Path to the output `.par` file (will be generated).
- `tim_file::String`: Path to the `.tim` file with TOAs.
"""
struct TempoRunFiles
    work_dir::String
    par_file_input::String      # начальный шаблон
    par_file_output::String     # будет создан, возможно, в work_dir
    tim_file::String
end

function Base.show(io::IO, ::MIME"text/plain", f::TempoRunFiles)
    indent = get(io, :indent, 0)
    pad = " "^indent
    println(io, pad, "TempoRunFiles")
    println(io, pad, "  work_dir:         ", f.work_dir)
    println(io, pad, "  par_file_input:   ", f.par_file_input)
    println(io, pad, "  par_file_output:  ", f.par_file_output)
    println(io, pad, "  tim_file:         ", f.tim_file)
end

"""
Generates a default output filename for a `.par` file by appending `_out.par`.

Throws an error if the input file does not end with `.par`.

# Arguments
- `par_file_input::String`: Path to the original `.par` file.

# Returns
- `String`: Output file name with `_out.par` suffix.
"""
function default_par_file_output(par_file_input::String)
    if !endswith(par_file_input, ".par")
        error("Expected input par file to end with .par, got: $par_file_input")
    end
    return par_file_input[1:end-4] * "_out.par"
end

#------------------------- OPTIONS -------------------------

"""
Defines low-level options for a TEMPO run.

# Fields
- `tempo_version::AbstractTempoVersion`: Which TEMPO version to use.
- `flags::String`: Additional command-line flags to pass to TEMPO.
- `nits::Int`: Number of internal iterations (must be ≥ 1).
- `gain::Float64`: GAIN parameter controlling convergence damping.
"""
struct TempoExecutionOptions
    tempo_version::AbstractTempoVersion
    flags::String
    nits::Int
    gain::Float64
end

function Base.show(io::IO, ::MIME"text/plain", opts::TempoExecutionOptions)
    indent = get(io, :indent, 0)
    pad = " "^indent
    println(io, pad, "TempoExecutionOptions")
    println(io, pad, "  tempo_version: ", opts.tempo_version)
    println(io, pad, "  flags:         ", opts.flags)
    println(io, pad, "  nits:          ", opts.nits)
    println(io, pad, "  gain:          ", opts.gain)
end

#------------------------- MODIFIERS -------------------------

"""
Describes modifications applied to the run.

# Fields
- `override_params::Vector{GeneralTempoParameter}`: Parameters to override or inject into the `.par` file.
- `time_start::Union{Nothing, Float64}`: Optional lower bound on TOAs (MJD).
- `time_finish::Union{Nothing, Float64}`: Optional upper bound on TOAs (MJD).
"""
struct TempoRunModifiers
    override_params::Vector{GeneralTempoParameter}
    time_start::Union{Nothing, Float64}
    time_finish::Union{Nothing, Float64}
end

function Base.show(io::IO, ::MIME"text/plain", mods::TempoRunModifiers)
    indent = get(io, :indent, 0)
    pad = " "^indent
    println(io, pad, "TempoRunModifiers")

    println(io, pad, "  override_params: ", length(mods.override_params), " parameter(s)")
    for param in mods.override_params
        println(io, pad, param)
    end

    println(io, pad, "  time_start:      ", isnothing(mods.time_start) ? "nothing" : mods.time_start)
    println(io, pad, "  time_finish:     ", isnothing(mods.time_finish) ? "nothing" : mods.time_finish)
end

#------------------------- BEHAVIOR -------------------------

"""
Controls what is captured or saved during the run.

# Fields
- `write_output::Bool`: Whether to capture stdout/stderr from TEMPO.
- `write_residuals::Bool`: Whether to pass `-residuals` flag to TEMPO to save residuals at each iteration.
- `save_internal_iterations::Bool`: Whether to retain all intermediate iteration results.
- `save_residuals::Bool`: Whether to retain residual files for each iteration.
"""
struct TempoRunBehavior
    write_output::Bool
    write_residuals::Bool
    save_internal_iterations::Bool
    save_residuals::Bool
end

function Base.show(io::IO, ::MIME"text/plain", b::TempoRunBehavior)
    indent = get(io, :indent, 0)
    pad = " "^indent
    println(io, pad, "TempoRunBehavior")
    println(io, pad, "  write_output:             ", b.write_output)
    println(io, pad, "  write_residuals:          ", b.write_residuals)
    println(io, pad, "  save_internal_iterations: ", b.save_internal_iterations)
    println(io, pad, "  save_residuals:           ", b.save_residuals)
end

#------------------------- BASIC SETTINGS -------------------------

"""
Encapsulates all settings needed for a basic TEMPO run.

This structure combines file paths, execution options, modifications,
and behavior control flags.

# Fields
- `files::TempoRunFiles`: Input/output file locations.
- `options::TempoExecutionOptions`: Core run-time options.
- `modifiers::TempoRunModifiers`: Pre-execution modifications.
- `behavior::TempoRunBehavior`: Output control settings.
"""
struct BasicTempoSettings <: AbstractTempoSettings
    files::TempoRunFiles
    options::TempoExecutionOptions
    modifiers::TempoRunModifiers
    behavior::TempoRunBehavior
end

function Base.show(io::IO, ::MIME"text/plain", s::BasicTempoSettings)
    indent = get(io, :indent, 0)
    pad = " " ^ indent
    next_indent = indent + 4

    println(io, pad, "BasicTempoSettings")

    println(io, pad, "  files:")
    show(IOContext(io, :indent => next_indent), MIME"text/plain"(), s.files)

    println(io, pad, "  options:")
    show(IOContext(io, :indent => next_indent), MIME"text/plain"(), s.options)

    println(io, pad, "  modifiers:")
    show(IOContext(io, :indent => next_indent), MIME"text/plain"(), s.modifiers)

    println(io, pad, "  behavior:")
    show(IOContext(io, :indent => next_indent), MIME"text/plain"(), s.behavior)
end

#------------------------- CONVENIENT CONSTRUCTOR -------------------------

"""
Convenience constructor for `BasicTempoSettings`.

Allows keyword-style construction with defaults for optional fields.
Raises an error if `nits < 1`.

# Keyword Arguments
- `work_dir`, `par_file_input`, `par_file_output`, `tim_file` — see `TempoRunFiles`
- `tempo_version`, `flags`, `nits`, `gain` — see `TempoExecutionOptions`
- `override_params`, `time_start`, `time_finish` — see `TempoRunModifiers`
- `write_output`, `write_residuals`, `save_internal_iterations`, `save_residuals` — see `TempoRunBehavior`
"""
function BasicTempoSettings(;
    work_dir::String,
    par_file_input::String,
    par_file_output::String = default_par_file_output(par_file_input),
    tim_file::String,
    tempo_version::AbstractTempoVersion,
    flags::String = "",
    nits::Int = 1,
    gain::Real = 1.0,
    override_params::Vector{GeneralTempoParameter} = GeneralTempoParameter[],
    time_start::Union{Nothing, Float64} = nothing,
    time_finish::Union{Nothing, Float64} = nothing,
    write_output::Bool = true,
    write_residuals::Bool = true,
    save_internal_iterations::Bool = false,
    save_residuals::Bool = false
)
    if nits < 1
        error("Number of iterations (nits) must be ≥ 1, got $nits")
    end

    files = TempoRunFiles(work_dir, par_file_input, par_file_output, tim_file)
    options = TempoExecutionOptions(tempo_version, flags, nits, gain)
    modifiers = TempoRunModifiers(override_params, time_start, time_finish)
    behavior = TempoRunBehavior(write_output, write_residuals, save_internal_iterations, save_residuals)
    return BasicTempoSettings(files, options, modifiers, behavior)
end

#------------------------- SHOW -------------------------

# """
# Pretty-print the settings in a human-readable format.
# """
# function Base.show(io::IO, settings::BasicTempoSettings)
#     println(io, "BasicTempoSettings:")
#     println(io, "  Work dir: ", settings.files.work_dir)
#     println(io, "  Input par file: ", settings.files.par_file_input)
#     println(io, "  Output par file: ", settings.files.par_file_output)
#     println(io, "  Tim file: ", settings.files.tim_file)
#     println(io, "  TEMPO version: ", settings.options.tempo_version)
#     println(io, "  Flags: ", settings.options.flags)
#     println(io, "  Number of iterations: ", settings.options.nits)
#     println(io, "  Gain parameter: ", settings.options.gain)
#     println(io, "  Override tempo parameters: ", length(settings.modifiers.override_params))
#     println(io, "  Time start: ", settings.modifiers.time_start)
#     println(io, "  Time finish: ", settings.modifiers.time_finish)
#     println(io, "  Write output: ", settings.behavior.write_output)
#     println(io, "  Write residuals: ", settings.behavior.write_residuals)
#     println(io, "  Save iterations: ", settings.behavior.save_internal_iterations)
#     println(io, "  Save residuals: ", settings.behavior.save_residuals)
# end

#------------------------- COPY -------------------------

"""
Create a modified copy of `BasicTempoSettings` with optional keyword overrides.

Any subset of fields may be changed via keyword arguments.

# Example
```julia
new_settings = copy_settings(old_settings; nits=10, gain=0.5)
```
"""
function copy_settings(s::BasicTempoSettings; kwargs...)
    work_dir        = get(kwargs, :work_dir, s.files.work_dir)
    par_file_input  = get(kwargs, :par_file_input, s.files.par_file_input)
    par_file_output = get(kwargs, :par_file_output, s.files.par_file_output)
    tim_file        = get(kwargs, :tim_file, s.files.tim_file)

    tempo_version   = get(kwargs, :tempo_version, s.options.tempo_version)
    flags           = get(kwargs, :flags, s.options.flags)
    nits            = get(kwargs, :nits, s.options.nits)
    gain            = get(kwargs, :gain, s.options.gain)

    override_params = get(kwargs, :override_params, s.modifiers.override_params)
    time_start      = get(kwargs, :time_start, s.modifiers.time_start)
    time_finish     = get(kwargs, :time_finish, s.modifiers.time_finish)

    write_output             = get(kwargs, :write_output, s.behavior.write_output)
    write_residuals          = get(kwargs, :write_residuals, s.behavior.write_residuals)
    save_internal_iterations = get(kwargs, :save_internal_iterations, s.behavior.save_internal_iterations)
    save_residuals           = get(kwargs, :save_residuals, s.behavior.save_residuals)

    files     = TempoRunFiles(work_dir, par_file_input, par_file_output, tim_file)
    options   = TempoExecutionOptions(tempo_version, flags, nits, gain)
    modifiers = TempoRunModifiers(override_params, time_start, time_finish)
    behavior  = TempoRunBehavior(write_output, write_residuals, save_internal_iterations, save_residuals)

    return BasicTempoSettings(files, options, modifiers, behavior)
end