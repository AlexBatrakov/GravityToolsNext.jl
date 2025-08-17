# TempoParFile.jl
# Representation and I/O for TEMPO/TEMPO2 .par files.

using Printf

# Assumes the following are available from TempoParameters.jl:
# - GeneralTempoParameter
# - TP = GeneralTempoParameter
# - get_par_file_representation(::GeneralTempoParameter)
# - extract_tempo_parameter_from_line(::AbstractString)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

# Ensure ".par" extension is present
_ensure_par_name(name::AbstractString) = endswith(lowercase(name), ".par") ? String(name) : String(name) * ".par"

# ------------------------------------------------------------
# TempoParFile
# ------------------------------------------------------------

"""
    TempoParFile

A representation of a TEMPO `.par` file, containing:
- `name`   : file name (e.g., `"B1937.par"`)
- `dir`    : directory path
- `path`   : full path (`joinpath(dir, name)`)
- `params` : dictionary of parameters (`Symbol => GeneralTempoParameter`)
- `order`  : ordered list of parameter symbols to preserve file layout
"""
struct TempoParFile
    name::String
    dir::String
    path::String
    params::Dict{Symbol,GeneralTempoParameter}
    order::Vector{Symbol}

    function TempoParFile(
        name::String,
        dir::String,
        params::Dict{Symbol,GeneralTempoParameter} = Dict{Symbol,GeneralTempoParameter}(),
        order::Vector{Symbol} = Symbol[]
    )
        name2 = _ensure_par_name(name)
        path  = joinpath(dir, name2)
        return new(name2, dir, path, params, order)
    end
end

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

"""
    Base.show(io::IO, par_file::TempoParFile)

Pretty-prints the `.par` file name, directory, and all parameters in order.
"""
function Base.show(io::IO, par_file::TempoParFile)
    println(io, "Tempo par file:")
    println(io, "  name: ", par_file.name)
    println(io, "  dir : ", par_file.dir)
    for (i, name_symbol) in enumerate(par_file.order)
        p = par_file.params[name_symbol]
        print(io, "    ", get_par_file_representation(p))
        i < length(par_file.order) && print(io, "\n")
    end
    return nothing
end

# ------------------------------------------------------------
# I/O
# ------------------------------------------------------------

"""
    read_par_file!(par_file::TempoParFile) -> TempoParFile

Reads and parses the `.par` file on disk.
- Ignores blank lines and comment lines beginning with `C`, `c`, or `#` (after trimming).
- Clears and fills `params` and `order` in-place.
- If a parameter appears multiple times, the **first occurrence** defines the position in `order`,
  subsequent occurrences update the value in-place.
"""
function read_par_file!(par_file::TempoParFile)
    empty!(par_file.params)
    empty!(par_file.order)

    open(par_file.path, "r") do io
        for (ln, raw) in enumerate(eachline(io))
            line = strip(raw)
            # skip blank + comment lines:
            if isempty(line) || occursin(r"^(?:#|[Cc](?:\s|$))", line)
                continue
            end

            param = try
                extract_tempo_parameter_from_line(line)
            catch err
                @warn "Failed to parse .par line $ln: $raw" exception=err
                continue
            end

            sym = param.name_symbol
            if !haskey(par_file.params, sym)
                push!(par_file.order, sym)
            end
            par_file.params[sym] = param
        end
    end
    return par_file
end

"""
    write_par_file!(par_file::TempoParFile) -> TempoParFile

Writes parameters to `par_file.path`, preserving the order defined in `order`.
If `par_file.dir` does not exist, it will be created.
"""
function write_par_file!(par_file::TempoParFile)
    isdir(par_file.dir) || mkpath(par_file.dir)
    open(par_file.path, "w") do io
        for name_symbol in par_file.order
            param = par_file.params[name_symbol]
            println(io, get_par_file_representation(param))
        end
    end
    return par_file
end

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

"""
    generate_par_file_name(base_name::String, suffix::String) -> String

Generates a new `.par` filename by removing the existing extension (if any),
appending `_suffix`, and adding `.par`.
"""
function generate_par_file_name(base_name::String, suffix::String)
    # strip extension if present
    base = replace(base_name, r"\.[Pp][Aa][Rr]$" => "")
    return "$(base)_$(suffix).par"
end

"""
    copy_par_file(par_file::TempoParFile;
                  new_name::Union{Nothing,String}=nothing,
                  suffix::Union{Nothing,String}=nothing,
                  new_dir::String = par_file.dir,
                  deep_copy::Bool = false) -> TempoParFile

Returns a copy of the `TempoParFile`, optionally changing:
- the filename (`new_name`) or by appending `suffix`,
- the directory (`new_dir`).

By default, the parameter storage is reused (no deep copy), matching the current behavior.
Set `deep_copy=true` to duplicate parameter dict and order vector.
"""
function copy_par_file(
    par_file::TempoParFile;
    new_name::Union{Nothing,String}=nothing,
    suffix::Union{Nothing,String}=nothing,
    new_dir::String = par_file.dir,
    deep_copy::Bool = false
)
    name′ = isnothing(new_name) ? par_file.name : new_name
    if isnothing(new_name) && !isnothing(suffix)
        name′ = generate_par_file_name(par_file.name, suffix)
    end
    name′ = _ensure_par_name(name′)

    params′ = deep_copy ? deepcopy(par_file.params) : par_file.params
    order′  = deep_copy ? deepcopy(par_file.order)  : par_file.order

    return TempoParFile(name′, new_dir, params′, order′)
end

# ------------------------------------------------------------
# Updates
# ------------------------------------------------------------

"""
    update_one_parameter_in_par_file!(par_file::TempoParFile, param::GeneralTempoParameter) -> TempoParFile

Adds or updates a parameter in-place:
- if `param.value === nothing`, updates only the `flag` of an existing parameter (throws if missing);
- otherwise, upserts the full parameter and preserves order (new names are appended).
"""
function update_one_parameter_in_par_file!(par_file::TempoParFile, param::GeneralTempoParameter)
    sym = param.name_symbol
    if param.value === nothing
        @assert haskey(par_file.params, sym) "Cannot set flag for a non-existent parameter: $(param.name)"
        par_file.params[sym].flag = param.flag
    else
        if !haskey(par_file.params, sym)
            push!(par_file.order, sym)
        end
        par_file.params[sym] = param
    end
    return par_file
end

"""
    update_par_file!(par_file::TempoParFile;
                     override_params::Vector{GeneralTempoParameter}=TP[],
                     nits::Union{Int,Nothing}=nothing,
                     gain::Union{Real,Nothing}=nothing,
                     time_start::Union{Float64,Nothing}=nothing,
                     time_finish::Union{Float64,Nothing}=nothing) -> TempoParFile

Applies parameter overrides and optionally sets:
- number of iterations (`NITS`),
- fitting gain factor (`GAIN`),
- time window (`START` / `FINISH`, flagged with `1`).

All changes are applied in-place to `par_file`.
"""
function update_par_file!(
    par_file::TempoParFile;
    override_params::Vector{GeneralTempoParameter}=TP[],
    nits::Union{Int,Nothing}=nothing,
    gain::Union{Real,Nothing}=nothing,
    time_start::Union{Float64,Nothing}=nothing,
    time_finish::Union{Float64,Nothing}=nothing
)
    for p in override_params
        update_one_parameter_in_par_file!(par_file, p)
    end
    if nits !== nothing
        update_one_parameter_in_par_file!(par_file, TP("NITS", Int(nits)))
    end
    if gain !== nothing
        update_one_parameter_in_par_file!(par_file, TP("GAIN", gain))
    end
    if time_start !== nothing
        update_one_parameter_in_par_file!(par_file, TP("START", time_start, flag=1))
    end
    if time_finish !== nothing
        update_one_parameter_in_par_file!(par_file, TP("FINISH", time_finish, flag=1))
    end
    return par_file
end