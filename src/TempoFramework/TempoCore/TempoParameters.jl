# TempoParameters.jl
# Utilities for representing and parsing TEMPO/TEMPO2 .par parameters.

using Printf
using Statistics

# Keep higher precision for BigFloat arithmetic/printing (consistent with your previous code)
setprecision(BigFloat, 80)

# ------------------------------------------------------------------------------
# Type
# ------------------------------------------------------------------------------

"""
    GeneralTempoParameter

Parameter from a TEMPO/TEMPO2 `.par` file.

Fields
- `name::String`        — canonical parameter name (may include spaces, e.g. "JUMP -be SCAMP")
- `name_symbol::Symbol` — cached symbol version of `name`
- `value::Union{Int64, BigFloat, String, Nothing}`
- `flag::Union{Int64, Nothing}`           — fit flag (typically -1, 0, 1), optional
- `uncertainty::Union{BigFloat, Nothing}` — 1σ uncertainty, optional
"""
mutable struct GeneralTempoParameter
    name::String
    name_symbol::Symbol
    value::Union{Int64, BigFloat, String, Nothing}
    flag::Union{Int64, Nothing}
    uncertainty::Union{BigFloat, Nothing}
end

const TP = GeneralTempoParameter

"""
Convenient constructor. Accepts `AbstractString` for name (e.g. `SubString`), stores as `String`.
Promotes any `Real` value/uncertainty to `BigFloat`.
"""
function GeneralTempoParameter(
    name::AbstractString,
    value::Union{Int64, Real, BigFloat, String, Nothing}=nothing;
    flag::Union{Int64, Nothing}=nothing,
    uncertainty::Union{Real, BigFloat, Nothing}=nothing
)
    promoted = value isa Integer      ? Int(value) :
               value isa AbstractFloat ? BigFloat(value) :
               value
    return GeneralTempoParameter(
        String(name),
        Symbol(String(name)),
        promoted,
        flag,
        uncertainty === nothing ? nothing : BigFloat(uncertainty)
    )
end

# Keep `name`/`name_symbol` in sync; promote numerics for `value`/`uncertainty`.
function Base.setproperty!(p::GeneralTempoParameter, f::Symbol, v)
    if f === :name
        setfield!(p, :name, String(v))
        setfield!(p, :name_symbol, Symbol(p.name))
        return v
    elseif f === :name_symbol
        setfield!(p, :name_symbol, Symbol(v))
        setfield!(p, :name, String(p.name_symbol))
        return v
    elseif f === :value
        v2 = v isa Integer ? Int(v) :
             v isa AbstractFloat ? BigFloat(v) : v
        setfield!(p, :value, v2)
        return v
    elseif f === :uncertainty
        setfield!(p, :uncertainty, v === nothing ? nothing : BigFloat(v))
        return v
    else
        return setfield!(p, f, v)
    end
end

function Base.show(io::IO, param::GeneralTempoParameter)
    indent = get(io, :indent, 0)
    print(io, " "^indent, param.name)

    val = param.value
    if val isa BigFloat
        # compact but without losing precision to Float64
        if abs(val) < big"1e-3" || abs(val) > big"1e6"
            print(io, " ", @sprintf("%.10e", val))
        else
            print(io, " ", @sprintf("%.10g", val))
        end
    else
        print(io, " ", val)
    end

    if param.flag !== nothing
        print(io, " ", param.flag)
    end

    if param.uncertainty !== nothing
        u = param.uncertainty
        if u isa BigFloat
            if abs(u) < big"1e-3" || abs(u) > big"1e6"
                print(io, " ±", @sprintf("%.10e", u))
            else
                print(io, " ±", @sprintf("%.10f", u))
            end
        else
            print(io, " ±", u)
        end
    end
    return nothing
end

# ------------------------------------------------------------------------------
# Parsing
# ------------------------------------------------------------------------------

"""
    parse_tempo_parameter_field(str) -> Int64 | BigFloat | String

Try `Int64`, then `BigFloat`, otherwise return `String`.
"""
function parse_tempo_parameter_field(value_str::AbstractString)
    (v = tryparse(Int64, value_str))    !== nothing && return v
    (v = tryparse(BigFloat, value_str)) !== nothing && return v
    return String(value_str)
end

"""
    extract_tempo_parameter_from_line(line) -> GeneralTempoParameter

Parse one `.par` line using a simple grammar:

- `name` is either 1 token, or **3 tokens** if the first three tokens are non-numeric strings
  (e.g. "JUMP -be SCAMP", "TNEF -be SCAMP").
- Then go: `value` (String/Int64/BigFloat), optional `flag` (Int64), optional `uncertainty` (BigFloat).
- If there is one trailing token after value:
    * Int64 → it's a flag
    * BigFloat → it's an uncertainty
- If there are two trailing tokens after value and they look like `Int64` then `BigFloat`,
  treat them as `flag` and `uncertainty` respectively.
"""
function extract_tempo_parameter_from_line(line::AbstractString)
    toks = split(strip(line))
    n = length(toks)
    n == 0 && error("Empty .par line")

    parsed = parse_tempo_parameter_field.(toks)
    types  = typeof.(parsed)

    # name: 3-string tokens → 3-word name, else 1-word name
    if n >= 3 && types[1] === String && types[2] === String && types[3] === String
        n_name = 3
        name = String(join(toks[1:3], " "))
    else
        n_name = 1
        name = String(toks[1])
    end

    # value
    value = n > n_name ? parsed[n_name + 1] : nothing

    flag = nothing
    unc  = nothing
    rest = n - n_name - 1  # tokens after value

    if rest == 1
        t2 = parsed[n_name + 2]
        if t2 isa BigFloat
            unc = t2
        elseif t2 isa Int64
            flag = t2
        end
    elseif rest >= 2
        t2 = parsed[n_name + 2]
        t3 = parsed[n_name + 3]
        if (t2 isa Int64) && (t3 isa BigFloat)
            flag = t2
            unc  = t3
        elseif t2 isa BigFloat
            # value + uncertainty (no flag)
            unc = t2
        end
        # extra tokens (if any) are ignored deliberately
    end

    return GeneralTempoParameter(name, value; flag=flag, uncertainty=unc)
end

# ------------------------------------------------------------------------------
# Formatting to .par
# ------------------------------------------------------------------------------

const PAR_NAME_W   = 23
const PAR_VALUE_W  = 32
const PAR_FLAG_W   = 6
const PAR_UNCERT_W = 32

_pad(s::AbstractString, w::Int) = (length(s) >= w ? String(s) : String(s) * " "^(w - length(s)))

# Back-compat alias for any old calls
align_str(s::String, w::Int) = _pad(s, w)

"""
    get_par_file_representation(param) -> String

Return a fixed-width `.par` line: name, value, optional flag, optional uncertainty.
No precision is lost: `BigFloat` values are formatted directly as BigFloat.
"""
function get_par_file_representation(param::GeneralTempoParameter)
    # Negative numbers take one extra column; mirror your old alignment trick
    n_name  = (param.value isa BigFloat && param.value < 0) ? 22 : PAR_NAME_W
    n_value = (param.value isa BigFloat && param.value < 0) ? 33 : PAR_VALUE_W
    n_flag  = PAR_FLAG_W
    n_unc   = PAR_UNCERT_W

    line = align_str(param.name, n_name)

    # value (keep BigFloat, no downcast)
    if param.value isa BigFloat
        if abs(param.value) < big"1e-3" || abs(param.value) > big"1e6"
            line *= align_str(@sprintf("%.21e", param.value), n_value)
        else
            line *= align_str(@sprintf("%.21g", param.value), n_value)
        end
    else
        line *= align_str(string(param.value), n_value)
    end

    # flag
    if param.flag !== nothing
        line *= align_str(string(param.flag), n_flag)
    else
        line *= " "^n_flag
    end

    # uncertainty (keep BigFloat, no downcast; decide e/f by magnitude of uncertainty)
    if param.uncertainty !== nothing
        u = param.uncertainty
        if u isa BigFloat && (abs(u) < big"1e-3" || abs(u) > big"1e6")
            line *= align_str(@sprintf("%.21e", u), n_unc)
        else
            line *= align_str(@sprintf("%.21f", u), n_unc)
        end
    end

    return line
end

# ------------------------------------------------------------------------------
# Updates
# ------------------------------------------------------------------------------

"""
    update_or_add_tempo_parameter!(params, param) -> Vector{GeneralTempoParameter}

Update by name or append if not present. Keeps order. Modifies `params` in place.
"""
function update_or_add_tempo_parameter!(params::Vector{GeneralTempoParameter}, param::GeneralTempoParameter)
    @inbounds for i in eachindex(params)
        if params[i].name == param.name
            params[i] = param
            return params
        end
    end
    push!(params, param)
    return params
end

"""
    update_many_tempo_parameters!(params, new_params) -> Vector{GeneralTempoParameter}

Apply `update_or_add_tempo_parameter!` for each element of `new_params`.
"""
function update_many_tempo_parameters!(params::Vector{GeneralTempoParameter}, new_params::Vector{GeneralTempoParameter})
    @inbounds for p in new_params
        update_or_add_tempo_parameter!(params, p)
    end
    return params
end

"""
    update_params_by_dict!(params, dict::Dict{Symbol,GeneralTempoParameter}) -> Vector{GeneralTempoParameter}

Batch update by `name_symbol`. Preserves original order; missing keys are appended.
"""
function update_params_by_dict!(params::Vector{GeneralTempoParameter}, dict::Dict{Symbol,GeneralTempoParameter})
    idx = Dict{Symbol,Int}()
    @inbounds for (i, p) in pairs(params)
        idx[p.name_symbol] = i
    end
    for (k, p) in dict
        if haskey(idx, k)
            params[idx[k]] = p
        else
            push!(params, p)
            idx[k] = length(params)
        end
    end
    return params
end