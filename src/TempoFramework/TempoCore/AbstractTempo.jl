# AbstractTempo.jl
# Minimal interface and concrete types for Tempo/Tempo2 configuration.

# --------------------------------------------------------------------------------------------------------------
# Abstract type and minimal interface
# --------------------------------------------------------------------------------------------------------------

"""
    AbstractTempoVersion

Abstract base for different TEMPO flavors (Tempo / Tempo2).
Implementations must provide:

- `tempo_data_dir(v) :: String`    — data directory path (TEMPO* env target)
- `tempo_cmd_path(v) :: String`    — resolved executable path (or bare name)
- `tempo_cmd(v)      :: Cmd`       — ready-to-run command
- `tempo_env(v)      :: Dict{String,String}` — env vars to pass to the process
"""
abstract type AbstractTempoVersion end

# Default empty env; concrete versions override it
tempo_env(::AbstractTempoVersion) = Dict{String,String}()

# --------------------------------------------------------------------------------------------------------------
# Tempo (original)
# --------------------------------------------------------------------------------------------------------------

struct Tempo <: AbstractTempoVersion
    custom_data_directory::String
    Tempo() = new(default_tempo_dir())
    Tempo(custom_dir::String) = new(custom_dir)
end

function default_tempo_dir()
    dd = get(ENV, "TEMPO", nothing)
    dd === nothing && throw(ArgumentError("TEMPO environment variable is not set"))
    return dd
end

# Executable: TEMPO_CMD override → Sys.which("tempo") → "tempo"
function tempo_cmd_path(::Tempo)
    if (p = get(ENV, "TEMPO_CMD", nothing)) !== nothing
        return p
    end
    w = Sys.which("tempo")
    return w === nothing ? "tempo" : w
end

tempo_data_dir(t::Tempo) = t.custom_data_directory
tempo_env(t::Tempo) = Dict("TEMPO" => tempo_data_dir(t))
tempo_cmd(t::Tempo) = Cmd(`$(tempo_cmd_path(t))`; ignorestatus=false)

# Back-compat string API (instance-based)
get_tempo_command(v::Tempo) = tempo_cmd_path(v)

# --------------------------------------------------------------------------------------------------------------
# Tempo2
# --------------------------------------------------------------------------------------------------------------

struct Tempo2 <: AbstractTempoVersion
    custom_data_directory::String
    Tempo2() = new(default_tempo2_dir())
    Tempo2(custom_dir::String) = new(custom_dir)
end

function default_tempo2_dir()
    dd = get(ENV, "TEMPO2", nothing)
    dd === nothing && throw(ArgumentError("TEMPO2 environment variable is not set"))
    return dd
end

# Executable: TEMPO2_CMD override → Sys.which("tempo2") → "tempo2"
function tempo_cmd_path(::Tempo2)
    if (p = get(ENV, "TEMPO2_CMD", nothing)) !== nothing
        return p
    end
    w = Sys.which("tempo2")
    return w === nothing ? "tempo2" : w
end

tempo_data_dir(t::Tempo2) = t.custom_data_directory
tempo_env(t::Tempo2) = Dict("TEMPO2" => tempo_data_dir(t))
tempo_cmd(t::Tempo2) = Cmd(`$(tempo_cmd_path(t))`; ignorestatus=false)

# Back-compat string API (instance-based)
get_tempo_command(v::Tempo2) = tempo_cmd_path(v)

# --------------------------------------------------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------------------------------------------------

"""
    validate(v::AbstractTempoVersion) -> Bool

Check that data directory exists and the executable is available (or at least named).
Returns `true` if OK, otherwise throws.
"""
function validate(v::AbstractTempoVersion)
    dd = tempo_data_dir(v)
    isdir(dd) || throw(ArgumentError("Tempo data directory does not exist: $dd"))

    cp = tempo_cmd_path(v)

    # If it's an explicit path, warn if file missing; if it's a bare name, warn if not on PATH.
    if occursin('/' , cp) || occursin('\\', cp)
        isfile(cp) || @warn "Tempo executable path does not exist" path=cp
    else
        Sys.which(cp) === nothing && @warn "Tempo executable not found on PATH" exe=cp
    end
    return true
end

# --------------------------------------------------------------------------------------------------------------
# Pretty printing
# --------------------------------------------------------------------------------------------------------------

_tempo_version_name(v) = v isa Tempo ? "Tempo" : v isa Tempo2 ? "Tempo2" : "Tempo<?>"

"""
Pretty-print `Tempo`/`Tempo2`. Supports inline one-liner via IOContext key `:inline => true`.
"""
function Base.show(io::IO, ::MIME"text/plain", v::AbstractTempoVersion)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    inline = get(io, :inline, false)

    name = _tempo_version_name(v)
    data_dir = try
        normpath(tempo_data_dir(v))
    catch
        getfield(v, :custom_data_directory)
    end

    exe = try
        normpath(tempo_cmd_path(v))
    catch
        get_tempo_command(v)
    end

    if inline
        # single line; no parentheses
        print(io, name, "  data_dir=", data_dir, "  exe=", exe)
    else
        # multi-line; no parentheses
        println(io, pad, name)
        println(io, pad, "  data_dir: ", data_dir)
        print(  io, pad, "  exe     : ", exe)
    end
end

# Fallback show delegates to text/plain
function Base.show(io::IO, v::AbstractTempoVersion)
    show(IOContext(io, :indent => get(io, :indent, 0)), MIME"text/plain"(), v)
end