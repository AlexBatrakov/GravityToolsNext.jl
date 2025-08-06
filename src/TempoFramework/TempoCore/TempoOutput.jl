#--------------------------------------------------------------------------------------------------------------

"""
    struct BasicTempoOutput

Contains basic statistical metrics from a single TEMPO2 fit iteration.

Fields:
- `rms_pre_fit_residual_us::Float64`: RMS of pre-fit residuals in microseconds.
- `rms_post_fit_residual_us::Float64`: RMS of post-fit residuals in microseconds.
- `rms_tn_post_fit_residual_us::Float64`: RMS of post-fit residuals after TN plugin, if present.
- `chisqr::Float64`: Total chi-square of the fit.
- `nfree::Int`: Number of degrees of freedom (data points minus fitted parameters).
- `chisqr_red::Float64`: Reduced chi-square (chisqr / nfree).
- `pre_post::Float64`: Ratio of pre-fit to post-fit RMS.
- `number_of_fit_parameters::Int`: Number of fitted parameters in the model.
- `number_of_points_in_fit::Int`: Number of TOAs included in the fit.
- `offset_value::Float64`: Best-fit offset value applied to TOAs.
- `offset_error::Float64`: Uncertainty of the offset value.
- `offset_e_sqrt_n::Float64`: Offset error scaled by sqrt(N).
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

function Base.show(io::IO, basic_output::BasicTempoOutput)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Basic Tempo Output:")
    println(io, ' '^(indent+4), "RMS Pre-fit Residual (us): ", basic_output.rms_pre_fit_residual_us)
    println(io, ' '^(indent+4), "RMS Post-fit Residual (us): ", basic_output.rms_post_fit_residual_us)
    println(io, ' '^(indent+4), "RMS TN Post-fit Residual (us): ", basic_output.rms_tn_post_fit_residual_us)
    println(io, ' '^(indent+4), "Chisqr: ", basic_output.chisqr)
    println(io, ' '^(indent+4), "Nfree: ", basic_output.nfree)
    println(io, ' '^(indent+4), "Chisqr/Nfree: ", basic_output.chisqr_red)
    println(io, ' '^(indent+4), "Pre/Post: ", basic_output.pre_post)
    println(io, ' '^(indent+4), "Number of Fit Parameters: ", basic_output.number_of_fit_parameters)
    println(io, ' '^(indent+4), "Number of Points in Fit: ", basic_output.number_of_points_in_fit)
    println(io, ' '^(indent+4), "Offset value: ", basic_output.offset_value)
    println(io, ' '^(indent+4), "Offset error: ", basic_output.offset_error)
    println(io, ' '^(indent+4), "Offset E * Sqrt(n): ", basic_output.offset_e_sqrt_n)
end

#--------------------------------------------------------------------------------------------------------------

"""
    struct FitParameter

Represents a single model parameter that was fitted by TEMPO2.

Fields:
- `name::String`: Full parameter name (e.g., "F0", "PB", etc).
- `name_symbol::Symbol`: Same name but as a symbol for fast lookup.
- `pre_fit::Float64`: Value of the parameter before fitting.
- `post_fit::Float64`: Value after fitting.
- `uncertainty::Float64`: Estimated 1σ uncertainty.
- `difference::Float64`: Difference between post-fit and pre-fit values.
- `fit_flag::Bool`: Whether the parameter was fitted (`true`) or held fixed (`false`).
"""
struct FitParameter
    name::String
    name_symbol::Symbol
    pre_fit::Float64
    post_fit::Float64
    uncertainty::Float64
    difference::Float64
    fit_flag::Bool

    function FitParameter(name::String, pre_fit::Float64, post_fit::Float64, uncertainty::Float64, difference::Float64, fit_flag::Bool)
        new(name, Symbol(name), pre_fit, post_fit, uncertainty, difference, fit_flag)
    end
end

function Base.show(io::IO, param::FitParameter)
    indent = get(io, :indent, 0)
    println(io, ' '^(indent), param.name, "  Pre-fit: ", param.pre_fit, "  Post-fit: ", param.post_fit, "  Unc: ", param.uncertainty, "  Diff: ", param.difference, "  fit: ", param.fit_flag)
    # и так далее для остальных полей
end



#--------------------------------------------------------------------------------------------------------------
"""
    struct TempoOutputError

Describes an error encountered while parsing TEMPO2 output.

Fields:
- `error_type::Symbol`: Type of the error (`:none`, `:missing_rms`, `:crash`, etc.).
- `message::String`: Explanation of the error.

Helper:
- `iserror(::TempoOutputError)`: Returns `true` if `error_type ≠ :none`.
"""
struct TempoOutputError
    error_type::Symbol
    message::String
end

TempoOutputError() = TempoOutputError(:none, "")

iserror(error::TempoOutputError) = error.error_type != :none

Base.:(==)(a::TempoOutputError, b::TempoOutputError) = 
    a.error_type == b.error_type && a.message == b.message

function Base.show(io::IO, error::TempoOutputError)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Tempo Output Error:")
    println(io, ' '^(indent+4), "Error type: ", error.error_type)
    println(io, ' '^(indent+4), "Error message: ", error.message)
end

#--------------------------------------------------------------------------------------------------------------

"""
    struct InternalIterationOutput

Stores the result of a single internal fit iteration from TEMPO2.

Fields:
- `basic::Union{BasicTempoOutput, Nothing}`: Overall fit metrics, or `nothing` if parsing failed.
- `fit_parameters::Union{Vector{FitParameter}, Nothing}`: Vector of fit parameters, or `nothing` if parsing failed.
- `error::TempoOutputError`: Error status of the parsing.
"""
struct InternalIterationOutput
    basic::Union{BasicTempoOutput, Nothing}
    fit_parameters::Union{Vector{FitParameter}, Nothing}
    error::TempoOutputError
end

InternalIterationOutput() = InternalIterationOutput(nothing, nothing, TempoOutputError())

function Base.show(io::IO, int_iter::InternalIterationOutput)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Internal iteration:")

    if !isnothing(int_iter.basic)
        show(IOContext(io, :indent => indent+4), int_iter.basic)
    else
        println(io, ' '^(indent+4), "No basic output.")
    end

    if !isnothing(int_iter.fit_parameters)
        println(io, ' '^(indent+4), "Fit Parameters:")
        for param in int_iter.fit_parameters
            show(IOContext(io, :indent => indent+8), param)
        end
    else
        println(io, ' '^(indent+4), "No fit parameters.")
    end

    if iserror(int_iter.error)
        show(IOContext(io, :indent => indent+4), int_iter.error)
    end
end

"""
    getindex(output::InternalIterationOutput, key::Union{String, Symbol})

Allows convenient lookup into `InternalIterationOutput`. 

If the key is a field of `BasicTempoOutput`, it returns the value directly.
If the key matches any parameter name, returns the corresponding `FitParameter`.

Throws a `KeyError` if no match is found.
"""
function Base.getindex(output::InternalIterationOutput, key::Union{String, Symbol})
    # Приведение к строке и символу
    key_sym = key isa Symbol ? key : Symbol(key)
    key_str = key isa String ? key : String(key)

    # Поиск среди полей basic
    if !isnothing(output.basic) && hasfield(BasicTempoOutput, key_sym)
        return getfield(output.basic, key_sym)
    end

    # Поиск среди параметров
    if !isnothing(output.fit_parameters)
        for param in output.fit_parameters
            if param.name_symbol == key_sym || param.name == key_str
                return param
            end
        end
    end

    throw(KeyError("Нет такого ключа '$key' ни в BasicTempoOutput, ни среди параметров"))
end

#--------------------------------------------------------------------------------------------------------------

"""
    parse_tempo_output(output::String, ::Type{Tempo2}) -> Vector{InternalIterationOutput}

Splits the full TEMPO2 output into iteration sections using the "Complete fit" delimiter,
then parses each section with `parse_internal_iteration_tempo_output`.

Returns a vector of `InternalIterationOutput`, each representing one iteration.
"""
function parse_tempo_output(output::String, ::Type{Tempo2})::Vector{InternalIterationOutput}
    # Разделение на блоки итераций
    output_split = split(output, "ss/fs")

    sections = ones(String, max(length(output_split)-1,1))
    sections[1] *= String(output_split[1])
    
    for iter in 2:length(output_split)
        sections[iter - 1] *= "ss/fs" * String(output_split[iter])
    end

    internal_iterations = InternalIterationOutput[]

    for (iter, section) in enumerate(sections)
        section = String(section)
        internal_iteration_output = parse_internal_iteration_tempo_output(section, Tempo2)

        # Проверка на наличие ошибок в секции

        # if internal_iteration_result.error !== TempoOutputError()
        #     break
        # end

        # Если ошибок нет, парсим результаты итерации

        push!(internal_iterations, internal_iteration_output)
    end

    return internal_iterations
end

"""
    parse_internal_iteration_tempo_output(section::String, ::Type{Tempo2}) -> InternalIterationOutput

Parses a single iteration block. Extracts error messages, then attempts to parse:
- basic metrics,
- fit parameters.

Returns `InternalIterationOutput`, containing result or error info.
"""
function parse_internal_iteration_tempo_output(section::String, ::Type{Tempo2})::InternalIterationOutput
    err = parse_tempo_output_error(section, Tempo2)
    if iserror(err)
        return InternalIterationOutput(nothing, nothing, err)
    end

    basic, err = parse_basic_tempo_output(section, Tempo2)
    if iserror(err)
        return InternalIterationOutput(nothing, nothing, err)
    end

    fit_params, err = parse_fit_parameters(section, Tempo2)
    if iserror(err)
        return InternalIterationOutput(nothing, nothing, err)
    end

    return InternalIterationOutput(basic, fit_params, TempoOutputError())  # без ошибки
end

const NAN_REGEX = r"(?i)nan"

"""
    safe_parse(x) -> Float64

Converts a string to Float64, returning `NaN` if the value is a case-insensitive "nan" or `nothing`.
"""
safe_parse(x) = x === nothing || occursin(NAN_REGEX, x) ? NaN : parse(Float64, x)

"""
    parse_basic_tempo_output(section::String, ::Type{Tempo2}) -> (BasicTempoOutput, TempoOutputError)

Parses basic fit statistics (RMS, chi-square, offset, etc.) from the output section.

Returns a tuple: parsed `BasicTempoOutput` or `nothing`, and a `TempoOutputError`.
"""
function parse_basic_tempo_output(section::String, ::Type{Tempo2})::Tuple{Union{BasicTempoOutput, Nothing}, TempoOutputError}
    
    rms_regex        = r"RMS pre-fit residual = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\), RMS post-fit residual = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\)"
    rms_tn_regex     = r"RMS post-fit residual TN = (\d+\.\d+|[-+]?[Nn][Aa][Nn]) \(us\)"
    chisq_regex      = r"Fit Chisq = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])\s+Chisqr/nfree = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])/(0|\d+)\s*=\s*(\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])\s+pre/post = (\d+\.?\d*[eE]?[-+]?\d*|[-+]?[Nn][Aa][Nn])"
    fit_params_regex = r"Number of fit parameters: (\d+)"
    points_regex     = r"Number of points in fit = (\d+)"
    offset_regex     = r"Offset: ([\d\.\-eE]+|[-+]?[Nn][Aa][Nn]) ([\d\.\-eE]+|[-+]?[Nn][Aa][Nn]) offset_e\*sqrt\(n\) = ([\d\.\-eE]+|[-+]?[Nn][Aa][Nn]) n = (\d+)"

    rms_match = match(rms_regex, section)
    rms_tn_match = match(rms_tn_regex, section)
    chisq_match = match(chisq_regex, section)
    fit_params_match = match(fit_params_regex, section)
    points_match = match(points_regex, section)
    offset_match = match(offset_regex, section)

    if rms_match === nothing
        return nothing, TempoOutputError(:missing_rms, "Не найден RMS pre/post блок")
    end
    if chisq_match === nothing
        return nothing, TempoOutputError(:missing_chisq, "Не найден блок chisqr/nfree и pre/post")
    end
    if fit_params_match === nothing
        return nothing, TempoOutputError(:missing_fit_params, "Не найдено число фитируемых параметров")
    end
    if points_match === nothing
        return nothing, TempoOutputError(:missing_points, "Не найдено число точек в fit-е")
    end
    if offset_match === nothing
        return nothing, TempoOutputError(:missing_offset, "Не найден offset блок")
    end

    # Теперь, когда всё найдено — безопасно строим структуру
    try
        return BasicTempoOutput(
            safe_parse(rms_match[1]),
            safe_parse(rms_match[2]),
            rms_tn_match !== nothing ? safe_parse(rms_tn_match[1]) : NaN,
            safe_parse(chisq_match[1]),
            parse(Int, chisq_match[3]),
            safe_parse(chisq_match[4]),
            safe_parse(chisq_match[5]),
            parse(Int, fit_params_match[1]),
            parse(Int, points_match[1]),
            safe_parse(offset_match[1]),
            safe_parse(offset_match[2]),
            safe_parse(offset_match[3])
        ), TempoOutputError()
    catch err
        return nothing, TempoOutputError(:parse_failure, "Не удалось распарсить числовые поля: $err")
    end
end

"""
    parse_fit_parameters(section::String, ::Type{Tempo2}) -> (Vector{FitParameter}, TempoOutputError)

Extracts the table of individual fitted parameters from a TEMPO2 iteration output.

Returns a tuple of parsed `FitParameter` vector and error status.
"""
function parse_fit_parameters(section::String, ::Type{Tempo2})::Tuple{Union{Vector{FitParameter}, Nothing}, TempoOutputError}
    lines = split(section, '\n')

    # Поиск границ таблицы параметров
    param_block_indices = findall(x -> all(c -> c == '-', x) && !isempty(x), lines)
    if length(param_block_indices) < 2
        return nothing, TempoOutputError(:missing_fit_table, "Не найдены границы таблицы параметров")
    end

    param_lines = lines[param_block_indices[1]+1 : param_block_indices[2]-1]

    regex = r"^(.*?)\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+([\d\.\-eE]+|[-+]?[Nn][Aa][Nn])\s+(Y|N)$"

    fit_parameters = FitParameter[]
    for line in param_lines
        m = match(regex, line)
        if m === nothing
            continue  # пропускаем строки, не соответствующие формату
        end

        name_raw = strip(m[1])
        name = replace(name_raw, r"\s*\(.*\)" => "")  # удаляем юниты

        try
            pre_fit     = safe_parse(m[2])
            post_fit    = safe_parse(m[3])
            uncertainty = safe_parse(m[4])
            difference  = safe_parse(m[5])
            fit_flag    = m[6] == "Y"

            if any(isnan, (pre_fit, post_fit, uncertainty, difference))
                return nothing, TempoOutputError(:nan_in_fit_param, "В строке '$line' содержится NaN")
            end

            push!(fit_parameters, FitParameter(name, pre_fit, post_fit, uncertainty, difference, fit_flag))

        catch err
            return nothing, TempoOutputError(:parse_failure, "Не удалось распарсить строку: '$line' ($err)")
        end
    end

    return fit_parameters, TempoOutputError()
end

"""
    parse_tempo_output_error(section::String, ::Type{Tempo2}) -> TempoOutputError

Scans an output section for known error patterns (segfaults, lost connection, panic, etc.).

Returns a `TempoOutputError` describing what was found, or `:none` if no issues.
"""
function parse_tempo_output_error(section::String, ::Type{Tempo2})::TempoOutputError
    lines = split(section, '\n')

    for line in lines
        s = strip(line)

        # Потеря соединения
        if occursin("Error: lost connection", s)
            return TempoOutputError(:lost_connection, s)
        end

        # Сегфолт или падение
        if occursin("Segmentation fault", s) || occursin("Aborted", s)
            return TempoOutputError(:crash, s)
        end

        # Прочие сообщения TEMPO2 об ошибке
        if occursin(r"(?i)(failed|cannot|unable to|panic)", s)
            return TempoOutputError(:runtime_error, s)
        end

        # Явная ошибка
        if occursin(r"(?i)^error", s)
            return TempoOutputError(:explicit_error, s)
        end
    end

    return TempoOutputError()
end

#--------------------------------------------------------------------------------------------------------------
