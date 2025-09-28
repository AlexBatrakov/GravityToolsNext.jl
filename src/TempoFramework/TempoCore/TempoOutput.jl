# TempoOutput.jl
# Parsing of TEMPO2 console output (per-iteration basic stats, fit table, and errors).

# --------------------------------------------------------------------------------------------------------------
# Basic structures
# --------------------------------------------------------------------------------------------------------------

"""
    BasicTempoOutput

Basic statistical metrics from a single TEMPO2 fit iteration.

Fields:
- `rms_pre_fit_residual_us::Float64`     : RMS of pre-fit residuals (µs)
- `rms_post_fit_residual_us::Float64`    : RMS of post-fit residuals (µs)
- `rms_tn_post_fit_residual_us::Float64` : RMS of post-fit residuals after TN plugin (µs), if present (NaN otherwise)
- `chisqr::Float64`                      : Total chi-square
- `nfree::Int`                           : Degrees of freedom
- `chisqr_red::Float64`                  : Reduced chi-square (`chisqr / nfree`)
- `pre_post::Float64`                    : `pre-fit RMS / post-fit RMS`
- `number_of_fit_parameters::Int`        : Number of fitted parameters
- `number_of_points_in_fit::Int`         : Number of TOAs in fit
- `offset_value::Float64`                : Best-fit offset value
- `offset_error::Float64`                : Uncertainty of the offset
- `offset_e_sqrt_n::Float64`             : `offset_error * sqrt(n)`
"""
struct BasicTempoOutput
    rms_pre_fit_residual_us::Float64
    rms_post_fit_residual_us::Float64
    rms_tn_post_fit_residual_us::Float64
    chisqr::Float64
    nfree::Int
    chisqr_red::Float64
    pre_post::Float64
    number_of_fit_parameters::Int
    number_of_points_in_fit::Int
    offset_value::Float64
    offset_error::Float64
    offset_e_sqrt_n::Float64
end

function Base.show(io::IO, basic::BasicTempoOutput)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    
    println(io, pad,  "Basic Tempo Output:")
    println(io, spad, "RMS Pre-fit Residual (µs): ", basic.rms_pre_fit_residual_us)
    println(io, spad, "RMS Post-fit Residual (µs): ", basic.rms_post_fit_residual_us)
    println(io, spad, "RMS TN Post-fit Residual (µs): ", basic.rms_tn_post_fit_residual_us)
    println(io, spad, "Chisq: ", basic.chisqr)
    println(io, spad, "Nfree: ", basic.nfree)
    println(io, spad, "Chisq/Nfree: ", basic.chisqr_red)
    println(io, spad, "Pre/Post: ", basic.pre_post)
    println(io, spad, "Number of Fit Parameters: ", basic.number_of_fit_parameters)
    println(io, spad, "Number of Points in Fit: ", basic.number_of_points_in_fit)
    println(io, spad, "Offset value: ", basic.offset_value)
    println(io, spad, "Offset error: ", basic.offset_error)
    print(io,   spad, "Offset E * sqrt(n): ", basic.offset_e_sqrt_n)
end

# --------------------------------------------------------------------------------------------------------------

"""
    FitParameter

Represents a single model parameter reported in the TEMPO2 fit table.

Fields:
- `name::String`         : parameter name (e.g., "F0", "PB", ...)
- `name_symbol::Symbol`  : cached symbol for fast lookup
- `pre_fit::Float64`     : pre-fit value
- `post_fit::Float64`    : post-fit value
- `uncertainty::Float64` : 1σ uncertainty
- `difference::Float64`  : `post_fit - pre_fit`
- `fit_flag::Bool`       : whether the parameter was fitted (`true`) or held fixed (`false`)
"""
struct FitParameter
    name::String
    name_symbol::Symbol
    pre_fit::Float64
    post_fit::Float64
    uncertainty::Float64
    difference::Float64
    fit_flag::Bool

    function FitParameter(name::String, pre_fit::Float64, post_fit::Float64,
                          uncertainty::Float64, difference::Float64, fit_flag::Bool)
        new(name, Symbol(name), pre_fit, post_fit, uncertainty, difference, fit_flag)
    end
end

function Base.show(io::IO, p::FitParameter)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)

    println(io, pad, p.name,
        "  Pre-fit: ", p.pre_fit,
        "  Post-fit: ", p.post_fit,
        "  Unc: ", p.uncertainty,
        "  Diff: ", p.difference,
        "  fit: ", p.fit_flag)
end

# --------------------------------------------------------------------------------------------------------------

"""
    TempoOutputError

Describes an error encountered while parsing one iteration section.

Fields:
- `error_type::Symbol` : `:none` if no error; otherwise a specific error tag (e.g., `:missing_rms`, `:crash`, ...)
- `message::String`    : diagnostic message
"""
struct TempoOutputError
    error_type::Symbol
    message::String
end

TempoOutputError() = TempoOutputError(:none, "")
iserror(e::TempoOutputError) = e.error_type != :none

Base.:(==)(a::TempoOutputError, b::TempoOutputError) =
    a.error_type == b.error_type && a.message == b.message

function Base.show(io::IO, e::TempoOutputError)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Tempo Output Error:")
    println(io, spad, "Error type: ", e.error_type)
    println(io, spad, "Error message: ", e.message)
end

# --------------------------------------------------------------------------------------------------------------

"""
    InternalIterationOutput

Container for one internal fit iteration (or pre-iteration preamble if the first iteration
did not actually run but produced logs/errors).

Fields:
- `basic::Union{BasicTempoOutput,Nothing}`          : basic fit metrics, or `nothing` if unavailable
- `fit_parameters::Union{Vector{FitParameter},Nothing}` : parameter table, or `nothing`
- `error::TempoOutputError`                         : error status
"""
struct InternalIterationOutput
    basic::Union{BasicTempoOutput, Nothing}
    fit_parameters::Union{Vector{FitParameter}, Nothing}
    error::TempoOutputError
end

InternalIterationOutput() = InternalIterationOutput(nothing, nothing, TempoOutputError())

function Base.show(io::IO, it::InternalIterationOutput)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Internal iteration:")

    if it.basic !== nothing
        show(IOContext(io, :indent => indent+4), it.basic)
    else
        println(io, ' '^(indent+4), "No basic output.")
    end

    if it.fit_parameters !== nothing
        println(io, ' '^(indent+4), "Fit Parameters:")
        for p in it.fit_parameters
            show(IOContext(io, :indent => indent+8), p)
        end
    else
        println(io, ' '^(indent+4), "No fit parameters.")
    end

    if iserror(it.error)
        show(IOContext(io, :indent => indent+4), it.error)
    end
end

"""
    getindex(output::InternalIterationOutput, key::Union{String,Symbol})

Convenience lookup. If `key` matches a field in `BasicTempoOutput`, returns it.
Otherwise, if it matches a parameter name, returns the corresponding `FitParameter`.
Throws `KeyError` if nothing matches.
"""
function Base.getindex(output::InternalIterationOutput, key::Union{String, Symbol})
    key_sym = key isa Symbol ? key : Symbol(key)
    key_str = key isa String ? key : String(key)

    if output.basic !== nothing && hasfield(BasicTempoOutput, key_sym)
        return getfield(output.basic, key_sym)
    end

    if output.fit_parameters !== nothing
        for p in output.fit_parameters
            if p.name_symbol == key_sym || p.name == key_str
                return p
            end
        end
    end

    throw(KeyError("No such key '$key' in BasicTempoOutput fields or fit parameters"))
end

# --------------------------------------------------------------------------------------------------------------
# Parsing utilities and regexes
# --------------------------------------------------------------------------------------------------------------

# Precompiled regex patterns (do not allocate on each call)
const _NAN_RE               = r"(?i)^\s*nan\s*$"
const _ITER_SPLIT_ANCHOR_RE = r"(?m)^(?:Complete fit|RMS pre-fit residual|ss/fs)"

const _RMS_RE        = r"RMS pre-fit residual = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\), RMS post-fit residual = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\)"
const _RMS_TN_RE     = r"RMS post-fit residual TN = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\)"
const _CHISQ_RE      = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])\s+Chisqr/nfree = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])/(0|\d+)\s*=\s*(\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])\s+pre/post = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])"
const _NFIT_RE       = r"Number of fit parameters:\s+(\d+)"
const _NPTS_RE       = r"Number of points in fit\s*=\s*(\d+)"
const _OFFSET_RE     = r"Offset:\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+offset_e\*sqrt\(n\)\s*=\s*([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+n\s*=\s*(\d+)"

# 1st group: name (kept verbatim, minus unit annotations)
# Columns: pre, post, unc, diff, fitted(Y/N)
const _PARAM_ROW_RE  = r"^(.*?)\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+(Y|N)$"

# Safe numeric parsing: accept NaN (case-insensitive) as NaN
safe_parse(x::AbstractString) = occursin(_NAN_RE, x) ? NaN : parse(Float64, x)
safe_parse(x) = x === nothing ? NaN : parse(Float64, x)

# --------------------------------------------------------------------------------------------------------------
# Iteration splitting and top-level parse
# --------------------------------------------------------------------------------------------------------------

"""
    _split_iterations(output::AbstractString) -> Vector{SubString{String}}

Split the full output into iteration sections by anchor lines, BUT keep any preamble
BEFORE the first anchor as part of the FIRST section (so early errors are captured).
Anchors supported: lines starting with `Complete fit`, `RMS pre-fit residual`, or `ss/fs`.
"""
function _split_iterations(output::AbstractString)::Vector{SubString{String}}
    anchors = collect(eachmatch(_ITER_SPLIT_ANCHOR_RE, output))

    # No anchors at all -> single section = entire output
    isempty(anchors) && return [SubString(String(output))]

    sections = SubString{String}[]
    for i in 1:length(anchors)
        start_i = (i == 1) ? firstindex(output) : anchors[i].offset
        end_i   = (i < length(anchors)) ? (anchors[i+1].offset - 1) : lastindex(output)
        push!(sections, SubString(String(output), start_i, end_i))
    end
    return sections
end

"""
    parse_tempo_output(output::String, ::Type{Tempo2}) -> Vector{InternalIterationOutput}

Split full TEMPO2 output into iteration sections, keeping any preamble BEFORE the first
anchor as part of the FIRST section, and parse each section.
"""
function parse_tempo_output(output::String, ::Type{Tempo2})::Vector{InternalIterationOutput}
    sections = _split_iterations(output)
    results = InternalIterationOutput[]
    for sec in sections
        push!(results, parse_internal_iteration_tempo_output(String(sec), Tempo2))
    end
    return results
end

# --------------------------------------------------------------------------------------------------------------
# Section parsing
# --------------------------------------------------------------------------------------------------------------

"""
    parse_internal_iteration_tempo_output(section::String, ::Type{Tempo2}) -> InternalIterationOutput

Parse a single iteration section:
1. detect known fatal errors;
2. parse basic fit statistics;
3. parse fit-parameter table.

If any step fails, returns an output with `error` set, others `nothing`.
"""
function parse_internal_iteration_tempo_output(section::String, ::Type{Tempo2})::InternalIterationOutput
    err_err = parse_tempo_output_error(section, Tempo2)
    # if iserror(err)
    #     return InternalIterationOutput(nothing, nothing, err)
    # end

    basic, err_basic = parse_basic_tempo_output(section, Tempo2)
    # if iserror(err)
    #     return InternalIterationOutput(nothing, nothing, err)
    # end

    params, err_params = parse_fit_parameters(section, Tempo2)
    # if iserror(err)
    #     return InternalIterationOutput(nothing, nothing, err)
    # end

    err_total = iserror(err_err) ? err_err : (iserror(err_basic) ? err_basic : (iserror(err_params) ? err_params : TempoOutputError()))

    return InternalIterationOutput(basic, params, err_total)
end

"""
    parse_basic_tempo_output(section::String, ::Type{Tempo2})
        -> (BasicTempoOutput | nothing, TempoOutputError)

Parse RMS, chi-square, number of points/params, and offset block.
"""
function parse_basic_tempo_output(section::String, ::Type{Tempo2})::Tuple{Union{BasicTempoOutput, Nothing}, TempoOutputError}
    rms_match    = match(_RMS_RE, section)
    rms_tn_match = match(_RMS_TN_RE, section)
    chisq_match  = match(_CHISQ_RE, section)
    nfit_match   = match(_NFIT_RE, section)
    npts_match   = match(_NPTS_RE, section)
    offset_match = match(_OFFSET_RE, section)

    if rms_match === nothing
        return nothing, TempoOutputError(:missing_rms, "RMS pre/post block not found")
    end
    if chisq_match === nothing
        return nothing, TempoOutputError(:missing_chisq, "Chisq/nfree block not found")
    end
    if nfit_match === nothing
        return nothing, TempoOutputError(:missing_fit_params, "Number of fit parameters not found")
    end
    if npts_match === nothing
        return nothing, TempoOutputError(:missing_points, "Number of points in fit not found")
    end
    if offset_match === nothing
        return nothing, TempoOutputError(:missing_offset, "Offset block not found")
    end

    try
        return BasicTempoOutput(
            safe_parse(rms_match[1]),
            safe_parse(rms_match[2]),
            rms_tn_match !== nothing ? safe_parse(rms_tn_match[1]) : NaN,
            safe_parse(chisq_match[1]),
            parse(Int, chisq_match[3]),
            safe_parse(chisq_match[4]),
            safe_parse(chisq_match[5]),
            parse(Int, nfit_match[1]),
            parse(Int, npts_match[1]),
            safe_parse(offset_match[1]),
            safe_parse(offset_match[2]),
            safe_parse(offset_match[3])
        ), TempoOutputError()
    catch err
        return nothing, TempoOutputError(:parse_failure, "Failed to parse numeric fields: $err")
    end
end

"""
    parse_fit_parameters(section::String, ::Type{Tempo2})
        -> (Vector{FitParameter} | nothing, TempoOutputError)

Extract the parameter table between dashed lines. If dashed separators are missing,
attempt a fallback search using a header-like line.
"""
function parse_fit_parameters(section::String, ::Type{Tempo2})::Tuple{Union{Vector{FitParameter}, Nothing}, TempoOutputError}
    lines = split(section, '\n')

    # find dashed separator lines (>= 8 dashes)
    stripdash = l -> strip(l)
    isdash    = l -> !isempty(l) && all(==('-'), l) && length(l) >= 8

    dashed_idx = findall(i -> isdash(stripdash(lines[i])), eachindex(lines))

    param_lines = nothing
    if length(dashed_idx) >= 2
        param_lines = lines[(dashed_idx[1]+1):(dashed_idx[2]-1)]
    else
        # Fallback: try to locate a header line and the following dashed separators
        hdr_idx = findfirst(l -> occursin(r"(?i)parameter\s+pre[-\s]?fit", l), lines)
        if hdr_idx === nothing
            return nothing, TempoOutputError(:missing_fit_table, "Parameter table not found")
        end
        dash_after = findfirst(i -> isdash(stripdash(lines[i])), (hdr_idx+1):length(lines))
        dash_after === nothing && return nothing, TempoOutputError(:missing_fit_table, "Parameter table separator not found")
        dash_after2 = findfirst(i -> isdash(stripdash(lines[i])), (dash_after+1):length(lines))
        dash_after2 === nothing && return nothing, TempoOutputError(:missing_fit_table, "Parameter table end separator not found")
        param_lines = lines[(dash_after+1):(dash_after2-1)]
    end

    params = FitParameter[]
    for line in param_lines
        m = match(_PARAM_ROW_RE, line)
        m === nothing && continue

        name_raw = strip(m[1])
        name = replace(name_raw, r"\s*\(.*\)" => "")  # drop unit annotations if present

        pre_fit     = safe_parse(m[2])
        post_fit    = safe_parse(m[3])
        uncertainty = safe_parse(m[4])
        difference  = safe_parse(m[5])
        fit_flag    = m[6] == "Y"

        if any(isnan, (pre_fit, post_fit, uncertainty, difference))
            return nothing, TempoOutputError(:nan_in_fit_param, "NaN encountered in parameter row: '$line'")
        end

        push!(params, FitParameter(name, pre_fit, post_fit, uncertainty, difference, fit_flag))
    end

    return params, TempoOutputError()
end

"""
    parse_tempo_output_error(section::String, ::Type{Tempo2}) -> TempoOutputError

Scan a section for known error patterns (lost connection, segfaults, panic, etc.).
"""
function parse_tempo_output_error(section::String, ::Type{Tempo2})::TempoOutputError
    for raw in split(section, '\n')
        s = strip(raw)

        if occursin("Error: lost connection", s)
            return TempoOutputError(:lost_connection, s)
        end
        if occursin("Segmentation fault", s) || occursin("Aborted", s)
            return TempoOutputError(:crash, s)
        end
        if occursin(r"(?i)(failed|cannot|unable to|panic)", s)
            return TempoOutputError(:runtime_error, s)
        end
        if occursin(r"(?i)^error", s)
            return TempoOutputError(:explicit_error, s)
        end
    end
    return TempoOutputError()
end

# --------------------------------------------------------------------------------------------------------------