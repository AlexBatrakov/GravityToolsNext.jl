# TempoParameters.jl
# Utilities for representing and parsing TEMPO/TEMPO2 .par parameters.

# Keep higher precision for BigFloat arithmetic/printing (consistent with your previous code)
setprecision(BigFloat, 80)

# ------------------------------------------------------------------------------
# Type
# ------------------------------------------------------------------------------

"""
    TempoParameter

Parameter from a TEMPO/TEMPO2 `.par` file.

Fields
- `name::String`        — canonical parameter name (may include spaces, e.g. "JUMP -be SCAMP")
- `name_symbol::Symbol` — cached symbol version of `name`
- `value::Union{Int64, BigFloat, String, Nothing}`
- `flag::Union{Int64, Nothing}`           — fit flag (typically -1, 0, 1), optional
- `uncertainty::Union{BigFloat, Nothing}` — 1σ uncertainty, optional
"""
mutable struct TempoParameter
    name::String
    name_symbol::Symbol
    value::Union{Int64, BigFloat, String, Nothing}
    flag::Union{Int64, Nothing}
    uncertainty::Union{BigFloat, Nothing}
end

const TP = TempoParameter

const _ALLOWED_FLAGS = (-1, 0, 1)

validate_flag(flag::Union{Nothing,Int}) =
    (flag === nothing || flag in _ALLOWED_FLAGS) ||
    throw(ArgumentError("flag must be -1, 0, or 1; got $flag"))

"""
Convenient constructor. Accepts `AbstractString` for name (e.g. `SubString`), stores as `String`.
Promotes any `Real` value/uncertainty to `BigFloat`.
"""
function TempoParameter(
    name::AbstractString,
    value::Union{Real,String,Nothing}=nothing;
    flag::Union{Int64, Nothing}=nothing,
    uncertainty::Union{Real,Nothing}=nothing,
)
    s = String(name)
    v2 = value === nothing ? nothing :
         value isa Integer ? Int(value) :
         value isa String  ? value :
         BigFloat(value)
    validate_flag(flag)
    u2 = uncertainty === nothing ? nothing : BigFloat(uncertainty)
    return TempoParameter(s, Symbol(s), v2, flag, u2)
end

# Keep `name`/`name_symbol` in sync; promote numerics for `value`/`uncertainty`.
function Base.setproperty!(p::TempoParameter, f::Symbol, v)
    if f === :name
        s = String(v)
        setfield!(p, :name, s)
        setfield!(p, :name_symbol, Symbol(s))
        return v
    elseif f === :name_symbol
        sym = Symbol(v)
        setfield!(p, :name_symbol, sym)
        setfield!(p, :name, String(sym))
        return v
    elseif f === :value
        v2 = v === nothing ? nothing :
             v isa Integer ? Int(v) :
             v isa String  ? v :
             BigFloat(v)
        setfield!(p, :value, v2)
        return v
    elseif f === :flag
        validate_flag(v)
        setfield!(p, :flag, v === nothing ? nothing : Int(v))
        return v
    elseif f === :uncertainty
        setfield!(p, :uncertainty, v === nothing ? nothing : BigFloat(v))
        return v
    else
        return setfield!(p, f, v)
    end
end

Base.summary(io, param::TempoParameter) = print(io, param.name)

function Base.show(io::IO, param::TempoParameter)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    print(io, pad, param.name)

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
# Helpers
# ------------------------------------------------------------------------------

value_as_float(p::TempoParameter)::Union{Nothing,Float64} =
    p.value === nothing      ? nothing :
    p.value isa Int64        ? Float64(p.value) :
    p.value isa BigFloat     ? Float64(p.value) :
    nothing

uncertainty_as_float(p::TempoParameter)::Union{Nothing,Float64} =
    p.uncertainty === nothing ? nothing : Float64(p.uncertainty)

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
    extract_tempo_parameter_from_line(line) -> TempoParameter

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
        if t2 isa Int64
            if t2 in _ALLOWED_FLAGS
                flag = t2
            else
                unc = t2
            end
        elseif t2 isa BigFloat
            unc = t2
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

    return TempoParameter(name, value; flag=flag, uncertainty=unc)
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
Format numerics with up to 21 significant digits; increase if needed.
"""
function get_par_file_representation(param::TempoParameter)
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

# ---------------------------
# Utils for Vector{TempoParameter}
# ---------------------------

# Accept either Symbol or String for lookups
_to_sym(name::Union{Symbol,AbstractString}) = name isa Symbol ? name : Symbol(String(name))

"""
    param_index(params, name) -> Union{Int,Nothing}

Find index of parameter by name (compared with `p.name_symbol`). Returns `nothing` if missing.
"""
function param_index(params::AbstractVector{TempoParameter}, name::Union{Symbol,AbstractString})
    s = _to_sym(name)
    @inbounds for i in eachindex(params)
        if params[i].name_symbol === s
            return i
        end
    end
    return nothing
end

"""
    has_param(params, name) -> Bool
"""
has_param(params::AbstractVector{TempoParameter}, name::Union{Symbol,AbstractString}) =
    param_index(params, name) !== nothing

"""
    get_param(params, name; default=nothing) -> Union{TempoParameter,Nothing}

Return a parameter by name or `default` if missing.
"""
function get_param(params::AbstractVector{TempoParameter}, name::Union{Symbol,AbstractString}; default=nothing)
    i = param_index(params, name)
    return i === nothing ? default : params[i]
end

# ---------------------------
# Upsert (mutating / non-mutating)
# ---------------------------

"""
    upsert_param!(params, p) -> Vector{TempoParameter}

Update by name or append if not present. Keeps order. Mutates `params`.

Rules:
- If parameter with this name is new → push `p`.
- If exists → replace value/flag intelligently:
    * If `p.value === nothing`, keep the old value.
    * If `p.flag === nothing`, keep the old flag.
    * Otherwise take the provided values.
"""
function upsert_param!(params::Vector{T}, p::T)::Vector{T} where {T<:TempoParameter}
    i = param_index(params, p.name_symbol)
    if i === nothing
        push!(params, p)
    else
        old = params[i]
        params[i] = TempoParameter(
            p.name,
            p.name_symbol,
            p.value === nothing ? old.value : p.value,
            p.flag  === nothing ? old.flag  : p.flag,
            p.uncertainty === nothing ? old.uncertainty : p.uncertainty
        )
    end
    return params
end

"""
    upsert_params!(params, new_params) -> Vector{TempoParameter}

Apply `upsert_param!` for each element of `new_params`. Order of new elements is preserved.
"""
function upsert_params!(params::Vector{T}, new_params::AbstractVector{T})::Vector{T} where {T<:TempoParameter}
    @inbounds for p in new_params
        upsert_param!(params, p)
    end
    return params
end

"""
    with_upserted_params(params, new_params) -> Vector{TempoParameter}

Non-mutating version of `upsert_params!`.
"""
function with_upserted_params(params::AbstractVector{T}, new_params::AbstractVector{T})::Vector{T} where {T<:TempoParameter}
    out = collect(params)
    upsert_params!(out, new_params)
    return out
end

upsert!(params::Vector{TempoParameter}, name, value; flag=nothing, uncertainty=nothing) =
    upsert_param!(params, TP(String(name), value; flag=flag, uncertainty=uncertainty))

# ---------------------------
# Delete by name (mutating / non-mutating)
# ---------------------------

"""
    delete_param!(params, name) -> Bool

Delete parameter by name if present. Returns `true` if deleted.
"""
function delete_param!(params::Vector{TempoParameter}, name::Union{Symbol,AbstractString})::Bool
    i = param_index(params, name)
    i === nothing && return false
    deleteat!(params, i)
    return true
end


"""
    delete_params!(params, names) -> Int

Delete all parameters whose names are in `names`. Returns count deleted.
Preserves relative order of remaining items.
"""
function delete_params!(params::Vector{TempoParameter}, names)::Int
    want = Set(_to_sym.(collect(names)))
    before = length(params)
    filter!(p -> !(p.name_symbol in want), params)
    return before - length(params)
end

"""
    without_params(params, names) -> Vector{TempoParameter}

Non-mutating bulk deletion.
"""
function without_params(params::AbstractVector{T}, names)::Vector{T} where {T<:TempoParameter}
    want = Set(_to_sym.(collect(names)))
    out = Vector{T}()
    sizehint!(out, length(params))
    @inbounds for p in params
        p.name_symbol in want || push!(out, p)
    end
    return out
end
