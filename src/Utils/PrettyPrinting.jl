# src/Utils/PrettyPrinting.jl
# internal pretty-print helpers (not exported)


"""
    rpad_display(s, w) -> String

Pad string `s` on the right (using display width via `Unicode.textwidth`) to width `w`.
Returns a new `String`. If `s` already fits, returns `String(s)` unchanged.
"""
rpad_display(s::AbstractString, w::Integer) =
    (pad = w - Unicode.textwidth(s); pad <= 0 ? String(s) : String(s) * repeat(" ", pad))

"""
    lpad_display(s, w) -> String

Pad string `s` on the left (using display width via `Unicode.textwidth`) to width `w`.
Returns a new `String`. If `s` already fits, returns `String(s)` unchanged.
"""
lpad_display(s::AbstractString, w::Integer) =
    (pad = w - Unicode.textwidth(s); pad <= 0 ? String(s) : repeat(" ", pad) * String(s))

"""
    cpad_display(s, w) -> String

Center-pad string `s` to width `w` (using display width). Left pad gets the floor,
right pad the remainder. Returns a new `String`.
"""
function cpad_display(s::AbstractString, w::Integer)
    pad = w - Unicode.textwidth(s)
    pad <= 0 && return String(s)
    left = pad ÷ 2
    right = pad - left
    return repeat(" ", left) * String(s) * repeat(" ", right)
end

# Convert a "path-like" column specification into `Vector{Symbol}`:
# - `:a`, `"a.b"`, `(:a, :b)`, `[:a, :b]` → `[:a, :b]`

"""
    pathlike_to_syms(x) -> Vector{Symbol}

Normalize a column/field path specification to a vector of symbols.
Supported inputs:
- `Symbol`
- `AbstractString` with dots (e.g. `"stats.std"`)
- `Tuple`/`Vector{Symbol}`
Throws `ArgumentError` otherwise.
"""
function pathlike_to_syms(x)
    if x isa Symbol
        return [x]
    elseif x isa AbstractString
        return Symbol.(split(x, '.'))
    elseif x isa Tuple
        return collect(Symbol.(x))
    elseif x isa Vector{Symbol}
        return x
    else
        throw(ArgumentError("Unsupported pathlike type: $(typeof(x))"))
    end
end

# Traverse nested `getproperty` by a symbol path.

"""
    get_in(obj, path::AbstractVector{Symbol})

Follow `path` through nested `getproperty` calls: e.g. `get_in(x, [:stats, :std])`
is equivalent to `x.stats.std`.
"""
function get_in(obj, path::AbstractVector{Symbol})
    for p in path
        obj = getproperty(obj, p)
    end
    return obj
end

# Build an accessor triple from a column spec:
#   (label::String, getter::Function, key_for_dicts)
# `key` is used to index `formats`/`renderers`/`aligns`.

"""
    make_accessor(spec) -> (label::String, getter::Function, key)

Create an accessor triple from a column spec. Supported forms:

- `:field` — direct field access, label `"field"`, key `:field`
- `"a.b"` — nested path by dots, label `"a.b"`, key `"a.b"`
- `(:a, :b)` / `Vector{Symbol}` — nested path, label `"a.b"`, key `"a.b"`
- `label => path_or_function` — custom header; right side can be a function `(row)->val`
  or a path-like spec as above. `key` is the `label`.
- `Function` — use function both as getter and as `key`, label `"(f)"`.

Throws `ArgumentError` for unsupported types.
"""
function make_accessor(spec)
    if spec isa Symbol
        label = string(spec)
        getter = r -> getproperty(r, spec)
        key = spec
        return (label, getter, key)
    elseif spec isa AbstractString
        label = spec
        path = pathlike_to_syms(spec)
        getter = r -> get_in(r, path)
        key = spec
        return (label, getter, key)
    elseif spec isa Tuple || spec isa Vector{Symbol}
        path = pathlike_to_syms(spec)
        label = join(string.(path), ".")
        getter = r -> get_in(r, path)
        key = label
        return (label, getter, key)
    elseif spec isa Pair
        # Pair{<:AbstractString,<:Any} or Pair{Symbol,<:Any}
        label = string(spec.first)
        val = spec.second
        if val isa Function
            getter = val
        else
            # val is path-like
            path = pathlike_to_syms(val)
            getter = r -> get_in(r, path)
        end
        key = spec.first
        return (label, getter, key)
    elseif spec isa Function
        label = "(f)"
        getter = spec
        key = spec
        return (label, getter, key)
    else
        throw(ArgumentError("Unsupported column spec type: $(typeof(spec))"))
    end
end

"""
    print_aligned_table(io, rows;
        namecol, cols;
        indent=0, header=true,
        formats=Dict(), renderers=Dict(), aligns=Dict(),
        header_align=:auto,
        sort_by=nothing, order=:asc,
        names_override=nothing, name_label_override=nothing)

Flexible, allocation-friendly table printer for a **vector** of rows (typically structs).

**Parameters**
- `namecol` — how to obtain the left-most "name" column from each row. One of:
  `Symbol`, dotted `String` (e.g. `"stats.std"`), tuple `(:stats,:std)`, or a function `(row)->name`.
- `cols` — collection of column specs (right-hand side of the table). Each element can be:
  - `Symbol` — field name
  - dotted `String` — nested path (e.g. `"stats.std"`)
  - tuple / `Vector{Symbol}` — nested path
  - `label => getter` — custom label with either a function or a path-like spec
  - `Function` — getter; label becomes `"(f)"`
- `formats::Dict` — optional number formats (Printf) per-column, keyed by the column `key`
  (as returned by `make_accessor`). Defaults to `"%10.6f"`.
- `renderers::Dict` — optional custom cell renderers per-column: `v -> AbstractString`.
  Takes precedence over `formats`.
- `aligns::Dict` — optional alignment per-column: `:left | :right | :center`.
  Defaults to `:right` for numeric columns, `:center` for Bool, otherwise `:left`.
- `header_align` — alignment for header labels:
  - `:auto` (default): center for right/center columns, left for text columns
  - `:left`, `:center`, or `:match` (same as data alignment)
- Sorting: `sort_by` can be `nothing` (no sorting), `:__name__` (sort by resolved names),
  any path-like spec, or a function. `order` is `:asc` or `:desc`.
- Name overrides: `names_override::Vector{String}` can be used to directly supply the row names
  (useful when the name is not obtainable from the row itself). `name_label_override` sets the
  header text for the name column.

**Examples**
```julia
print_aligned_table(stdout, rows;
    namecol = :backend,
    cols = (
        :efac, :equad, "stats.std", "chisqr" => (:stats, :chisqr),
        :success, :converged,
    ),
    formats = Dict(
        :efac => "%10.6f", :equad => "%10.4f", "stats.std" => "%10.6f",
        "chisqr" => "%10.4f"
    ),
    renderers = Dict(
        :success => x -> (x ? "✓" : "✗"),
        :converged => x -> (x ? "✓" : "✗"),
    ),
    aligns = Dict(:success=>:center, :converged=>:center),
    header_align = :auto,
    sort_by = :ad_objective, order = :asc,
)
```
"""
function print_aligned_table(
    io::IO,
    rows::AbstractVector;
    namecol,
    cols,
    indent::Union{Int,Nothing} = nothing,
    header::Bool = true,
    formats::AbstractDict{<:Any,<:AbstractString} = Dict{Any,String}(),
    renderers::AbstractDict{<:Any,<:Function} = Dict{Any,Function}(),
    aligns::AbstractDict{<:Any,Symbol} = Dict{Any,Symbol}(),
    header_align::Symbol = :auto,
    # extras:
    sort_by::Union{Nothing,Symbol,AbstractString,Tuple,Function} = nothing,
    order::Symbol = :asc,
    names_override::Union{Nothing,Vector{String}} = nothing,
    name_label_override::Union{Nothing,String} = nothing,
)
    indent = indent === nothing ? get(io, :indent, 0) : indent
    isempty(rows) && (println(io, repeat(" ", indent), "(empty)"); return)

    # 1) resolve name column (and header)
    if names_override === nothing
        (name_label, name_get, _) = make_accessor(namecol)
        names = string.(name_get.(rows))
    else
        name_label = something(name_label_override, string(namecol))
        names = names_override
    end
    wname = maximum(Unicode.textwidth.(names))

    # 1.5) optional sort
    if sort_by !== nothing
        perm = if sort_by === :__name__
            sortperm(names; rev = (order === :desc))
        else
            (_, getter, _) = make_accessor(sort_by)
            sortperm(getter.(rows); rev = (order === :desc))
        end
        rows  = rows[perm]
        names = names[perm]
    end

    # 2) prepare columns
    col_specs   = collect(cols)
    ncols       = length(col_specs)
    header_lbls = Vector{String}(undef, ncols)
    col_cells   = Vector{Vector{String}}(undef, ncols)
    col_widths  = Vector{Int}(undef, ncols)
    col_aligns  = Vector{Symbol}(undef, ncols)
    col_keys    = Vector{Any}(undef, ncols)

    for (i, spec) in enumerate(col_specs)
        label, getter, key = make_accessor(spec)
        col_keys[i]     = key
        header_lbls[i]  = label
        vals            = getter.(rows)
        cells           = Vector{String}(undef, length(vals))

        r    = get(renderers, key, nothing)
        ffmt = Printf.Format(get(formats, key, "%10.6f"))

        for j in eachindex(vals)
            v = vals[j]
            s = if r !== nothing
                r(v)
            elseif v isa Real
                io_tmp = IOBuffer()
                Printf.format(io_tmp, ffmt, float(v))
                String(take!(io_tmp))
            elseif v isa Bool
                v ? "✓" : "✗"
            else
                string(v)
            end
            cells[j] = s
        end

        col_cells[i] = cells
        w = maximum(Unicode.textwidth.(cells))
        w = max(w, Unicode.textwidth(label))
        col_widths[i] = w

        alg = get(aligns, key,
                  any(v -> v isa Real, vals) ? :right  :
                  any(v -> v isa Bool, vals) ? :center : :left)
        col_aligns[i] = alg
    end

    pad = repeat(" ", indent)
    sep = "  "

    # header alignment helper
    _align_hdr(s, w, alg) = begin
        mode = header_align === :match  ? alg :
               header_align === :left   ? :left :
               header_align === :center ? :center :
               # :auto — center for right/center columns, left for text columns
               (alg === :right || alg === :center ? :center : :left)
        mode === :right  ? lpad_display(s, w) :
        mode === :center ? cpad_display(s, w) :
                           rpad_display(s, w)
    end

    # 3) header
    if header
        print(io, pad, rpad_display(name_label, wname))
        for i in 1:ncols
            print(io, sep, _align_hdr(header_lbls[i], col_widths[i], col_aligns[i]))
        end
        println(io)
    end

    # 4) rows
    for (i, name) in enumerate(names)
        print(io, pad, rpad_display(name, wname))
        for j in 1:ncols
            s   = col_cells[j][i]
            w   = col_widths[j]
            alg = col_aligns[j]
            print(io, sep,
                  alg === :right  ? lpad_display(s, w) :
                  alg === :center ? cpad_display(s, w) :
                                    rpad_display(s, w))
        end
        println(io)
    end
end

"""
    print_aligned_table(io, d::AbstractDict; kwargs...)

Dictionary-friendly overload. Treats **values** as rows, and uses **keys** as the
left-most "name" column. Keys are converted to strings. Supports the same options
as the vector overload, plus:

- `name_label`: header text for the name column (default: `"key"`).
- `sort_by = :__key__` sorts by dictionary keys (name column). Use `nothing` to keep
  insertion order, or any other valid `sort_by` (applied to values).

All other options (`cols`, `formats`, `renderers`, `aligns`, `header_align`, `order`) are forwarded.
"""
function print_aligned_table(
    io::IO,
    d::AbstractDict;
    namecol = :__key__,   # logical "label" — vector version does not read it,
                          # but we set label to "key" (or what you provide).
    cols,
    indent::Union{Int,Nothing} = nothing,
    header::Bool = true,
    formats::AbstractDict{<:Any,<:AbstractString} = Dict{Any,String}(),
    renderers::AbstractDict{<:Any,<:Function} = Dict{Any,Function}(),
    aligns::AbstractDict{<:Any,Symbol} = Dict{Any,Symbol}(),
    header_align::Symbol = :auto,
    sort_by::Union{Nothing,Symbol,AbstractString,Tuple,Function} = :__key__,
    order::Symbol = :asc,
    name_label::AbstractString = "key",
)
    indent = indent === nothing ? get(io, :indent, 0) : indent

    isempty(d) && (println(io, repeat(" ", indent), "(empty)"); return)

    ks = collect(keys(d))
    vs = collect(values(d))
    names = string.(ks)

    # sorting
    if sort_by === :__key__
        perm = sortperm(names; rev = (order === :desc))
    elseif sort_by === nothing
        perm = eachindex(vs)
    else
        (_, getter, _) = make_accessor(sort_by)
        perm = sortperm(getter.(vs); rev = (order === :desc))
    end

    vs    = vs[perm]
    names = names[perm]

    # forward to the vector version, overriding names
    print_aligned_table(io, vs;
        namecol = namecol,                  # not used to compute names, but appears in header if name_label not changed
        cols    = cols,
        indent  = indent,
        header  = header,
        formats = formats,
        renderers = renderers,
        aligns  = aligns,
        header_align = header_align,
        sort_by = nothing,                  # already sorted
        names_override = names,
        name_label_override = name_label,
    )
end

"""
    print_aligned_table(data; kwargs...)

Convenience overload that prints to `stdout`.
"""
print_aligned_table(data; kwargs...) = print_aligned_table(stdout, data; kwargs...)

# const INDENT_STEP = 2

# "Current indent level from IOContext."
# indent_level(io)::Int = get(io, :indent, 0)

# "Pad string for current indent + extra."
# pad(io, add::Int=0) = repeat(" ", indent_level(io) + add)

# "Child IOContext with increased indent; inherits :compact/:limit."
# with_indent(io; add::Int=INDENT_STEP,
#                  compact=get(io, :compact, false),
#                  limit=get(io, :limit, true),
#                  color=get(io, :color, false)) =
#     IOContext(io,
#         :indent  => indent_level(io) + add,
#         :compact => compact,
#         :limit   => limit,
#         :color   => color)

# "Print `key: value` with indent."
# print_kv(io, key, value; add::Int=INDENT_STEP) =
#     println(io, pad(io, add), key, ": ", value)

# "Print a titled section and show nested object."
# function print_section(io, title::AbstractString, obj;
#                        add::Int=INDENT_STEP, subadd::Int=INDENT_STEP+2,
#                        mime::MIME=MIME"text/plain"())
#     println(io, pad(io, add), title, ":")
#     show(with_indent(io; add=subadd), mime, obj)
#     println(io)  # trailing newline
# end

# "Compact numeric pretty-print (handy inside kv)."
# print_short_real(x::Real; sig::Int=6) = @sprintf("%.*g", sig, float(x))