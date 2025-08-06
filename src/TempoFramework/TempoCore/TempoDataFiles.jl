
# Тип для одной точки TOA
struct TimTOAEntry
    index::Int                # индекс после удаления комментариев
    timfile_line::Int         # исходный номер строки в .tim файле
    toa::Float64              # MJD
    freq::Float64             # Частота
    uncertainty::Float64      # В микросекундах
    backend::String           # Название бекенда
end

"""
    get_tim_line_number(clean_index::Int, discarded_lines::Vector{Int}) -> Int

Возвращает оригинальный номер строки в .tim файле для записи с индексом `clean_index`
(то есть после удаления комментариев и заголовка).
"""
function get_tim_line_number(clean_index::Int, discarded_lines::Vector{Int})
    line_number = clean_index + 2  # компенсируем заголовок из 2 строк
    for skipped in discarded_lines
        if line_number >= skipped
            line_number += 1
        end
    end
    return line_number
end

"""
    read_tim_file(path::String) -> Vector{TimTOAEntry}

Reads a .tim file and extracts clean TOA entries, skipping comment lines and preserving original line numbering.
"""
function read_tim_file(path::String)::Vector{TimTOAEntry}
    tim_file_data = readdlm(path, String)

    # Найдём закомментированные строки
    discarded_lines = filter(i -> tim_file_data[i, 1] == "C", 1:size(tim_file_data, 1))

    # Определим строки, не являющиеся комментариями
    is_not_comment(i) = tim_file_data[i, 1] != "C"
    valid_rows = filter(is_not_comment, 1:size(tim_file_data, 1))

    # Перебираем "чистые" строки данных
    N_TOAs = length(valid_rows) - 2  # первые две строки — заголовок и метаданные

    entries = TimTOAEntry[]
    for i in 1:N_TOAs
        timfile_line = get_tim_line_number(i, discarded_lines)
        row = tim_file_data[valid_rows[i + 2], :]  # сдвиг на 2 строки заголовка

        try
            toa         = parse(Float64, row[3])
            freq        = parse(Float64, row[2])
            uncertainty = parse(Float64, row[4])
            be_index    = findfirst(x -> x == "-be", row)
            backend     = isnothing(be_index) ? "unknown" : row[be_index + 1]

            push!(entries, TimTOAEntry(i, timfile_line, toa, freq, uncertainty, backend))
        catch err
            @warn "Не удалось распарсить строку $timfile_line: $row" exception=err
        end
    end

    return entries
end


# Тип для одной строки из residuals.dat
struct TempoResidualEntry
    time::Float64              
    signal::Float64
    residual::Float64
    residual_tn::Float64
    uncertainty::Float64
end

function read_residual_file(path::String)::Vector{TempoResidualEntry}
    residuals_data = open(path, "r") do io
        readlines(io)
    end

    entries = TempoResidualEntry[]

    for (i, line) in enumerate(residuals_data)
        fields = split(strip(line))
        if length(fields) != 5
            @warn "Пропущена строка $i: ожидалось 5 элементов, получено $(length(fields))"
            continue
        end

        try
            time         = parse(Float64, fields[1])
            signal       = parse(Float64, fields[2]) * 1e6
            residual     = parse(Float64, fields[3]) * 1e6
            residual_tn  = parse(Float64, fields[4]) * 1e6
            uncertainty  = parse(Float64, fields[5]) * 1e6

            push!(entries, TempoResidualEntry(time, signal, residual, residual_tn, uncertainty))
        catch err
            @warn "Ошибка при разборе строки $i: $line" exception=err
        end
    end

    return entries
end

struct CombinedTOAEntry
    index::Int          
    timfile_line::Int
    in_fit::Bool
    toa::Float64     
    freq::Float64        
    signal::Float64      
    red_noise::Float64
    residual::Float64     
    residual_tn::Float64    
    uncertainty_orig::Float64  
    uncertainty::Float64     
    backend::String        
    weight::Float64
end

function combine_tim_and_residuals(
    tim_entries::Vector{TimTOAEntry},
    residuals::Vector{TempoResidualEntry};
    time_start::Union{Nothing, Float64} = nothing,
    time_finish::Union{Nothing, Float64} = nothing
    )::Vector{CombinedTOAEntry}

    N = length(tim_entries)
    if length(residuals) != N
        error("Размеры TOA и резидуалов не совпадают: $(length(tim_entries)) ≠ $(length(residuals))")
    end

    entries = Vector{CombinedTOAEntry}(undef, N)

    for i in 1:N
        toa  = tim_entries[i]
        res  = residuals[i]

        in_fit = true
        if !isnothing(time_start) && toa.toa < time_start
            in_fit = false
        end
        if !isnothing(time_finish) && toa.toa > time_finish
            in_fit = false
        end

        weight = 1.0 / res.uncertainty^2

        entries[i] = CombinedTOAEntry(
            toa.index,
            toa.timfile_line,
            in_fit,
            toa.toa,
            toa.freq,
            res.signal,
            res.residual - res.residual_tn,
            res.residual,
            res.residual_tn,
            toa.uncertainty,
            res.uncertainty,
            toa.backend,
            weight
        )
    end

    return entries
end

