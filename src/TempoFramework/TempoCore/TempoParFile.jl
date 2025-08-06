#--------------------------------------------------------------------------------------------------------------
# tempo par files

"""
    TempoParFile

A representation of a TEMPO `.par` file, containing:
- `name`: file name (e.g., `"B1937.par"`),
- `dir`: directory path,
- `path`: full path (`joinpath(dir, name)`),
- `params`: dictionary of parameters (`Symbol => GeneralTempoParameter`),
- `order`: ordered list of parameter symbols to preserve file layout.
"""
struct TempoParFile
    name::String
    dir::String
    path::String
    params::Dict{Symbol,GeneralTempoParameter}
    order::Vector{Symbol}

    function TempoParFile(name::String, dir::String, params::Dict{Symbol,GeneralTempoParameter} = Dict{Symbol,GeneralTempoParameter}(), order::Vector{Symbol} = Symbol[])
        path = joinpath(dir, name)
        return new(name, dir, path, params, order)
    end
end


"""
    Base.show(io::IO, par_file::TempoParFile)

Pretty-prints the `.par` file name, directory, and all parameters in order, using formatted output.
"""
function Base.show(io::IO, par_file::TempoParFile)
    println(io, "Tempo par file: ")
    println(io, "   name: ", par_file.name)
    println(io, "   dir: ", par_file.dir)
    for (i, name_symbol) in enumerate(par_file.order)
        print("    ", get_par_file_representation(par_file.params[name_symbol]))
        if i < length(par_file.order)
            print("\n")
        end
    end
	return nothing
end

"""
    read_par_file!(par_file::TempoParFile) -> TempoParFile

Reads and parses the `.par` file on disk.
Ignores comments and blank lines.
Clears and fills the `params` dictionary and `order` vector in-place.
"""
function read_par_file!(par_file::TempoParFile)
    empty!(par_file.params)
    empty!(par_file.order)

    open(par_file.path, "r") do file_in
        for line in eachline(file_in)
            if startswith(line, "C ") || startswith(line, "c ") || startswith(line, " ")
                continue
            end
            param = extract_tempo_parameter_from_line(line)
            par_file.params[param.name_symbol] = param
            push!(par_file.order, param.name_symbol)
        end
    end
    return par_file
end

"""
    write_par_file!(par_file::TempoParFile) -> TempoParFile

Writes the parameter list in `par_file` to its associated file on disk,
preserving the order defined in `order`.
"""
function write_par_file!(par_file::TempoParFile)
    open(par_file.path, "w") do file_out
        for name_symbol in par_file.order
            param = par_file.params[name_symbol]
            println(file_out, get_par_file_representation(param))
        end
    end
    return par_file
end

"""
    generate_par_file_name(base_name::String, suffix::String) -> String

Generates a new `.par` filename by removing the `.par` extension from `base_name`,
appending the `suffix`, and restoring the `.par` extension.
"""
function generate_par_file_name(base_name::String, suffix::String)
    base_name_wihout_extension = base_name[1:end-4]
    return base_name_wihout_extension * "_" * suffix * ".par"
end

"""
    copy_par_file(par_file::TempoParFile; new_name=nothing, suffix=nothing, new_dir=...) -> TempoParFile

Returns a copy of the `TempoParFile`, optionally changing:
- the filename (`new_name` or `suffix`),
- the directory (`new_dir`).

The parameter values and order are reused as-is (not deep-copied).
"""
function copy_par_file(par_file::TempoParFile; new_name::Union{Nothing, String}=nothing, suffix::Union{Nothing, String}=nothing, new_dir::String = par_file.dir)
    if isnothing(new_name) && !isnothing(suffix)
        new_name = generate_par_file_name(par_file.name, suffix)
    end
    return TempoParFile(new_name, new_dir, par_file.params, par_file.order)
end

"""
    update_par_file!(par_file::TempoParFile;
                     override_params=Vector{GeneralTempoParameter},
                     nits=nothing, gain=nothing,
                     time_start=nothing, time_finish=nothing) -> TempoParFile

Applies a list of parameter overrides and optionally sets:
- number of iterations (`NITS`),
- fitting gain factor (`GAIN`),
- start time (`START`) and finish time (`FINISH`).

All changes are applied in-place to `par_file.params`.
"""
function update_par_file!(par_file::TempoParFile; 
                          override_params::Vector{GeneralTempoParameter},
                          nits::Union{Int64, Nothing} = nothing, 
                          gain::Union{Real, Nothing} = nothing,  
                          time_start::Union{Float64, Nothing} = nothing,
                          time_finish::Union{Float64, Nothing} = nothing)
    for param in override_params
        update_one_parameter_in_par_file!(par_file, param)
    end
    if !isnothing(nits)
        update_one_parameter_in_par_file!(par_file, TP("NITS", nits))
    end
    if !isnothing(gain)
        update_one_parameter_in_par_file!(par_file, TP("GAIN", Float64(gain)))
    end
    if !isnothing(time_start)
        update_one_parameter_in_par_file!(par_file, TP("START", time_start, flag = 1))
    end
    if !isnothing(time_finish)
        update_one_parameter_in_par_file!(par_file, TP("FINISH", time_finish, flag = 1))
    end
    return par_file
end

"""
    update_one_parameter_in_par_file!(par_file::TempoParFile, param::GeneralTempoParameter)

Adds or updates a parameter in-place:
- if the value is `nothing`, updates the `flag` only;
- otherwise, replaces or adds the parameter fully and preserves order.
"""
function update_one_parameter_in_par_file!(par_file::TempoParFile, param::GeneralTempoParameter)
    if param.value == nothing
        par_file.params[param.name_symbol].flag = param.flag
    else
        if !haskey(par_file.params, param.name_symbol)
            push!(par_file.order, param.name_symbol)
        end
        par_file.params[param.name_symbol] = param
    end
    return par_file
end


# function detect_backends(tim_file_path::String)
#     backends = Set{String}()  # Используем множество для избежания дубликатов
#     open(tim_file_path, "r") do file
#         for line in eachline(file)
#             if occursin("-be ", line)
#                 backend_name = match(r"-be (\S+)", line)
#                 if backend_name !== nothing
#                     # Извлекаем только название бэкенда, исключая "-be"
#                     push!(backends, backend_name.captures[1])
#                 end
#             end
#         end
#     end
#     return collect(backends)  # Преобразуем множество в вектор
# end

