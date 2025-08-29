# TempoParFile.jl
# Representation and I/O for TEMPO/TEMPO2 .par files.

# Assumes the following are available from TempoParameters.jl:
# - TempoParameter
# - TP = TempoParameter
# - get_par_file_representation(::TempoParameter)
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
- `params` : dictionary of parameters (`Symbol => TempoParameter`)
- `order`  : ordered list of parameter symbols to preserve file layout
"""
struct TempoParFile
    name::String
    dir::String
    path::String
    params::Dict{Symbol,TempoParameter}
    order::Vector{Symbol}

    function TempoParFile(
        name::String,
        dir::String,
        params::Dict{Symbol,TempoParameter} = Dict{Symbol,TempoParameter}(),
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
function Base.show(io::IO, ::MIME"text/plain", pf::TempoParFile)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "TempoParFile")
    println(io, spad, "name: ", pf.name)
    println(io, spad, "dir:  ", pf.dir)
    if !isempty(pf.order)
        println(io, spad, "params:")
        pio = IOContext(io, :indent => indent + 4)
        for sym in pf.order
            haskey(pf.params, sym) || continue
            println(pio, get_par_file_representation(pf.params[sym]))
        end
    end
end

# ------------------------------------------------------------
# I/O
# ------------------------------------------------------------

# Keep only part before a trailing '#', if any (no effect if '#' absent)
_strip_inline_comment(line::AbstractString) = begin
    i = findfirst('#', line)
    i === nothing ? line : strip(first(line, i-1))
end


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
            line0 = strip(raw)
            # skip blank + comment lines at start:
            if isempty(line0) || occursin(r"^(?:#|[Cc](?:\s|$))", line0)
                continue
            end
            # optionally strip trailing inline comments:
            line = _strip_inline_comment(line0)
            isempty(line) && continue

            param = try
                extract_tempo_parameter_from_line(line)
            catch err
                @warn "Failed to parse .par line $ln" text=line0 exception=err
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
    copy_par_file(par_file; new_name=nothing, suffix=nothing, new_dir=par_file.dir, deep_copy=false)

Return a (possibly shallow) copy of `par_file` with optional file/dir changes.

Precedence for filename:
1) if `new_name` is provided -> use it,
2) else if `suffix` is provided -> append `_suffix` before ".par",
3) else -> keep original name.
Ensures the result ends with `.par`.
"""
function copy_par_file(
    par_file::TempoParFile;
    new_name::Union{Nothing,String}=nothing,
    suffix::Union{Nothing,String}=nothing,
    new_dir::String = par_file.dir,
    deep_copy::Bool = false
)
    # compute final name with clear precedence
    name_new = isnothing(new_name) ? par_file.name : new_name
    if isnothing(new_name) && !isnothing(suffix)
        name_new = generate_par_file_name(par_file.name, suffix)
    end
    name_new = _ensure_par_name(name_new)

    params_new = deep_copy ? deepcopy(par_file.params) : par_file.params
    order_new  = deep_copy ? deepcopy(par_file.order)  : par_file.order

    return TempoParFile(name_new, new_dir, params_new, order_new)
end

# ---- Utilities for TempoParFile's (params::Dict, order::Vector) ----
# Reuse _to_sym and validate_flag you already have in TempoParameters.jl

# Ensure a symbol is present in the order vector (append if missing)
function _ensure_in_order!(pf::TempoParFile, sym::Symbol)
    # keep stable order; avoid duplicates
    (sym in pf.order) || push!(pf.order, sym)
    return pf
end

"""
    has_param(pf, name) -> Bool
"""
has_param(pf::TempoParFile, name::Union{Symbol,AbstractString}) =
    haskey(pf.params, _to_sym(name))

"""
    get_param(pf, name; default=nothing) -> Union{TempoParameter,Nothing}
"""
get_param(pf::TempoParFile, name::Union{Symbol,AbstractString}; default=nothing) =
    get(pf.params, _to_sym(name), default)

"""
    upsert_param!(pf, p) -> TempoParFile

Insert or replace parameter `p` by name, and ensure it is present in `order`.
"""
function upsert_param!(pf::TempoParFile, p::TempoParameter)
    p.value === nothing &&
        throw(ArgumentError("upsert_param! expects a value; use update_one_parameter_in_par_file! for flag-only"))
    pf.params[p.name_symbol] = p
    _ensure_in_order!(pf, p.name_symbol)
    return pf
end

"""
    set_flag!(pf, name, flag) -> TempoParFile

Set the fit-flag for an existing parameter. Throws if missing.
"""
function set_flag!(pf::TempoParFile, name::Union{Symbol,AbstractString}, flag::Int)
    validate_flag(flag)
    sym = _to_sym(name)
    haskey(pf.params, sym) || throw(ArgumentError("Cannot set flag for missing parameter: $name"))
    pf.params[sym].flag = flag
    return pf
end

"""
    upsert_params!(pf, ps) -> TempoParFile

Apply `upsert_param!` for each `p in ps`. If `p.value === nothing`, treat it as a flag-only update.
"""
function upsert_params!(pf::TempoParFile, ps::AbstractVector{TempoParameter})
    @inbounds for p in ps
        if p.value === nothing
            # flag-only update; require it to exist
            f = p.flag === nothing ? pf.params[p.name_symbol].flag : p.flag
            set_flag!(pf, p.name_symbol, f)
        else
            upsert_param!(pf, p)
        end
    end
    return pf
end

"""
    delete_param!(pf, name) -> Bool

Delete by name from `params` and remove from `order`. Returns `true` if deleted.
"""
function delete_param!(pf::TempoParFile, name::Union{Symbol,AbstractString})::Bool
    sym = _to_sym(name)
    deleted = pop!(pf.params, sym, nothing) !== nothing
    if deleted
        filter!(s -> s !== sym, pf.order)
    end
    return deleted
end

"""
    delete_params!(pf, names) -> Int

Delete all parameters listed in `names`. Returns the number of deletions.
"""
function delete_params!(pf::TempoParFile, names)::Int
    cnt = 0
    for n in names
        cnt += delete_param!(pf, n) ? 1 : 0
    end
    return cnt
end

"""
    params_as_vector(pf) -> Vector{TempoParameter}

Return parameters in the stored order; skips names that may have been removed from `params`.
"""
params_as_vector(pf::TempoParFile)::Vector{TempoParameter} =
    [pf.params[sym] for sym in pf.order if haskey(pf.params, sym)]

# ------------------------------------------------------------
# Updates
# ------------------------------------------------------------

"""
    update_one_parameter_in_par_file!(par_file::TempoParFile, param::TempoParameter) -> TempoParFile

Adds or updates a parameter in-place:
- if `param.value === nothing`, updates only the `flag` of an existing parameter (throws if missing);
- otherwise, upserts the full parameter and preserves order (new names are appended).
"""

function update_one_parameter_in_par_file!(pf::TempoParFile, p::TempoParameter)
    if p.value === nothing
        # flag-only update
        set_flag!(pf, p.name_symbol, something(p.flag, pf.params[p.name_symbol].flag))
    else
        upsert_param!(pf, p)
    end
    return pf
end

# Helper: build a case-insensitive set of override names
_overrides_name_set(overrides::AbstractVector{TempoParameter}) =
    Set(uppercase.(String.(getfield.(overrides, :name_symbol))))

"""
    update_par_file!(pf::TempoParFile; ...; prefer_overrides=true) -> TempoParFile

Apply overrides and optionally set NITS/GAIN/START/FINISH with a clear precedence rule.
"""
function update_par_file!(
    pf::TempoParFile;
    override_params::Vector{TempoParameter}=TP[],
    nits::Union{Int,Nothing}=nothing,
    gain::Union{Real,Nothing}=nothing,
    time_start::Union{Float64,Nothing}=nothing,
    time_finish::Union{Float64,Nothing}=nothing,
    prefer_overrides::Bool=true,
)
    # 1) Apply explicit overrides first (flag-only and full upserts)
    for p in override_params
        update_one_parameter_in_par_file!(pf, p)
    end

    # 2) special knobs with a single uniform rule
    names_ci = _overrides_name_set(override_params)

    nits        !== nothing && (!prefer_overrides || !("NITS"   in names_ci)) && update_one_parameter_in_par_file!(pf, TP("NITS",   Int(nits)))
    gain        !== nothing && (!prefer_overrides || !("GAIN"   in names_ci)) && update_one_parameter_in_par_file!(pf, TP("GAIN",   gain))
    time_start  !== nothing && (!prefer_overrides || !("START"  in names_ci)) && update_one_parameter_in_par_file!(pf, TP("START",  time_start;  flag=1))
    time_finish !== nothing && (!prefer_overrides || !("FINISH" in names_ci)) && update_one_parameter_in_par_file!(pf, TP("FINISH", time_finish; flag=1))

    return pf
end