# --------------------------------------------------------------------------------------------------------------
# GeneralTempoParameter: Representation of parameters in TEMPO/TEMPO2 .par files
#
# This module defines the main parameter type used to represent individual timing model parameters in
# TEMPO/TEMPO2 .par files. Each parameter can carry a value, flag, and uncertainty.
#
# Provides:
# - A structured, optionally mutable type to store parameter values
# - Parsing logic to read from .par files and convert values correctly
# - Pretty-printing and serialization to .par file format
# - Utility functions for parameter update and merge

# --------------------------------------------------------------------------------------------------------------

# Set precision for BigFloat to 80 digits
setprecision(BigFloat, 80)

"""
    GeneralTempoParameter

Represents a parameter from a TEMPO `.par` file. Each parameter consists of:
- `name`: full name of the parameter as a string,
- `name_symbol`: the corresponding symbol version of the name,
- `value`: parameter value (can be `Int`, `BigFloat`, `String`, or `nothing`),
- `flag`: optional integer flag used by TEMPO to fix/free the parameter,
- `uncertainty`: optional measurement uncertainty (stored as `BigFloat`).
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
    GeneralTempoParameter(name::String, value::Union{...}; flag=..., uncertainty=...)

Convenient constructor for `GeneralTempoParameter`, converting `Float64` values to `BigFloat`
and assigning optional `flag` and `uncertainty` fields.
"""
function GeneralTempoParameter(name::String, value::Union{Int64, Float64, BigFloat, String, Nothing}=nothing; flag::Union{Int64, Nothing}=nothing, uncertainty::Union{Float64, BigFloat, Nothing}=nothing)
    name_symbol = Symbol(name)  # Convert name to symbol
    big_value = isa(value, Float64) ? BigFloat(value) : value
    big_uncertainty = isnothing(uncertainty) ? nothing : BigFloat(uncertainty)
    return GeneralTempoParameter(name, name_symbol, big_value, flag, big_uncertainty)
end

"""
    Base.setproperty!(x::GeneralTempoParameter, f::Symbol, v::Float64)

Automatically converts assigned `Float64` values to `BigFloat` when modifying `value` or `uncertainty`.
"""
function Base.setproperty!(x::GeneralTempoParameter, f::Symbol, v::Float64)
    if f === :value || f === :uncertainty
        v = BigFloat(v)
    end
    setfield!(x, f, v)
end


# Constructor to convert from ValueVariable type
# GeneralTempoParameter(var::ValueVariable) = GeneralTempoParameter(var.name, var.value)

# Display function for GeneralTempoParameter
function Base.show(io::IO, param::GeneralTempoParameter)
    indent = get(io, :indent, 0)
    print(io, " "^indent, param.name)

    # Formatted output for BigFloat values
    if param.value isa BigFloat
        print(io, " ", @sprintf("%.10g", param.value))
    else
        print(io, " ", param.value)
    end

    if param.flag !== nothing
        print(io, " ", param.flag)
    end

    # Formatted output for BigFloat uncertainties
    if param.uncertainty !== nothing
        if param.uncertainty isa BigFloat
            print(io, " ±", @sprintf("%.10g", param.uncertainty))
        else
            print(io, " ±", param.uncertainty)
        end
    end

    return nothing
end

"""
    parse_tempo_parameter_field(str)

Parses a string value from a `.par` line field and returns an `Int64`, `BigFloat`, or `String`.
Used during line parsing for type inference.
"""
function parse_tempo_parameter_field(value_str)
    value_int64 = tryparse(Int64, value_str)
    if value_int64 !== nothing
        return value_int64
    end

    value_bigfloat = tryparse(BigFloat, value_str)
    if value_bigfloat !== nothing
        return value_bigfloat
    end
    
    return String(value_str)
end

"""
    extract_tempo_parameter_from_line(line::String) -> GeneralTempoParameter

Parses a line from a `.par` file into a `GeneralTempoParameter`. Supports parameters with:
- single-word or three-word names,
- optional flag and uncertainty,
- scientific or fixed-point formats.
"""
function extract_tempo_parameter_from_line(line::String)
    line_split = split(line)
    n = length(line_split)
    line_parsed = parse_tempo_parameter_field.(line_split)
    line_parsed_types = typeof.(line_parsed)

    # Название параметра: 1 или 3 слова
    if n >= 3 && line_parsed_types[1:3] == [String, String, String]
        n_name = 3
        name = join(line_split[1:3], " ")
    else
        n_name = 1
        name = String(line_split[1])
    end

    value = n > n_name ? line_parsed[n_name + 1] : nothing

    flag = nothing
    uncertainty = nothing

    # остались ≥ 2 слова? → value + flag
    if n > n_name + 1
        maybe_flag_or_uncertainty = line_parsed[n_name + 2]
        if isa(maybe_flag_or_uncertainty, BigFloat)
            uncertainty = maybe_flag_or_uncertainty
        else
            flag = maybe_flag_or_uncertainty
        end
    end

    # осталась ещё одна штука после флага? → uncertainty
    if n > n_name + 2 && flag !== nothing
        maybe_uncertainty = line_parsed[n_name + 3]
        if isa(maybe_uncertainty, BigFloat)
            uncertainty = maybe_uncertainty
        end
    end

    return GeneralTempoParameter(name, value, flag=flag, uncertainty=uncertainty)
end


"""
    align_str(s::String, width::Int) -> String

Pads a string with spaces to a fixed width. Used for alignment in `.par` file formatting.
"""
function align_str(s::String, n::Int)
    return s * " "^max(n - length(s), 1)
end

"""
    get_par_file_representation(param::GeneralTempoParameter) -> String

Returns a formatted string representation of a parameter suitable for writing to a `.par` file.
Includes name, value, flag (if present), and uncertainty (if present), using scientific notation when needed.
"""
function get_par_file_representation(param::GeneralTempoParameter)
    n_name  = (param.value isa BigFloat && param.value < 0) ? 22 : 23
    n_value = (param.value isa BigFloat && param.value < 0) ? 33 : 32
    n_flag = 6
    n_uncertainty = 32

    line = align_str(param.name, n_name)


    # If the value is a BigFloat, use @sprintf for formatting
    if param.value isa BigFloat
        if abs(param.value) < 1e-3 || abs(param.value) > 1e6
            line *= align_str(@sprintf("%.21e", param.value), n_value)
        else 
            line *= align_str(@sprintf("%.21g", param.value), n_value)
        end
    else
        line *= align_str(string(param.value), n_value)
    end

    if param.flag !== nothing
        line *= string(param.flag) * "  "
    else
        line *= "   "
    end

    # If uncertainty is a BigFloat, use @sprintf for formatting
    if param.uncertainty !== nothing
        if param.uncertainty isa BigFloat && (abs(param.uncertainty) < 1e-3 || abs(param.uncertainty) > 1e6)
            line *= align_str(@sprintf("%.21e", param.uncertainty), n_uncertainty)
        else
            line *= align_str(@sprintf("%.21f", param.uncertainty), n_uncertainty)
        end
    end

    return line
end

"""
    update_or_add_tempo_parameter!(params, param)

Updates an existing parameter in a vector `params` if found by name. Otherwise, appends it.
Returns the modified vector (same reference).
"""
function update_or_add_tempo_parameter!(params::Vector{GeneralTempoParameter}, param::GeneralTempoParameter)
    # Поиск параметра в списке
    found = false
    for (i, existing_param) in enumerate(params)
        if existing_param.name == param.name
            params[i] = param # Обновление значения параметра
            found = true
            break
        end
    end

    # Добавление параметра, если он не был найден
    if !found
        push!(params, param)
    end

    return params
end

"""
    update_many_tempo_parameters!(params, new_params)

Applies `update_or_add_tempo_parameter!` to each element of `new_params`.
Modifies `params` in-place.
"""
function update_many_tempo_parameters!(params::Vector{GeneralTempoParameter}, new_params::Vector{GeneralTempoParameter})
    for p in new_params
        update_or_add_tempo_parameter(params, p)
    end
end
