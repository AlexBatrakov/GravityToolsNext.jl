# =============================================================================
#  AdaptiveGridFramework / RefinementSettings.jl
#
#  Purpose
#  -------
#  Defines refinement units and a container (`RefinementSettings`) used by the
#  adaptive 2D grid engine. A *refinement unit* describes *where* and *how* to
#  refine along a named variable (axis). The container aggregates one or more
#  units together with execution options (desired refinement depth, parallel
#  preference) and a list of metrics to extract at each grid point.
#
#  Notes
#  -----
#  • This file is agnostic to any particular solver; it only encodes policies
#    for selecting candidate points during refinement.
#  • Pretty `show` methods for `text/plain` are provided for human‑readable
#    summaries in the REPL and logs.
# =============================================================================

# == Refinement Units ==========================================================

abstract type AbstractRefinementUnit end

"""
    LocalMinimaUnit(name; min=-Inf, max=Inf, from_min=true)

Refine until a local minimum is detected along variable `name` within the
interval `[min, max]`.

Arguments
- `name::Symbol` : the variable/axis identifier.
- `min`, `max`   : bounds of the search interval.
- `from_min`     : if `true`, sweep from the lower bound (otherwise from upper).
"""
struct LocalMinimaUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    from_min::Bool
end

LocalMinimaUnit(name; min = -Inf, max = Inf, from_min = true) =
    LocalMinimaUnit(name, min, max, from_min)

function Base.show(io::IO, ::MIME"text/plain", ru::LocalMinimaUnit)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Local minima unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value:    ", ru.min)
    println(io, spad, "Maximum value:    ", ru.max)
    println(io, spad, "From min value:   ", ru.from_min)
end

"""
    FullUnit(name; min=-Inf, max=Inf)

Always refine over the whole range `[min, max]` of variable `name`.
"""
struct FullUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
end

FullUnit(name; min = -Inf, max = Inf) = FullUnit(name, min, max)

function Base.show(io::IO, ::MIME"text/plain", ru::FullUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Full refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value:    ", ru.min)
    println(io, spad, "Maximum value:    ", ru.max)
end

"""
    DiffUnit(name; min=-Inf, max=Inf, diff, from_min=true)

Refine where the absolute change Δ(metric) exceeds `diff` within `[min, max]`.
If `from_min` is `true`, scanning proceeds from the lower bound.
"""
struct DiffUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    diff::Float64
    from_min::Bool
end

DiffUnit(name; min = -Inf, max = Inf, diff, from_min = true) =
    DiffUnit(name, min, max, diff, from_min)

function Base.show(io::IO, ::MIME"text/plain", ru::DiffUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Difference refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value:    ", ru.min)
    println(io, spad, "Maximum value:    ", ru.max)
    println(io, spad, "Maximal difference: ", ru.diff)
    println(io, spad, "From min value:   ", ru.from_min)
end

"""
    RelDiffUnit(name; min=-Inf, max=Inf, rel_diff, from_min=true)

Refine where the relative change exceeds `rel_diff` within `[min, max]`.
"""
struct RelDiffUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    rel_diff::Float64
    from_min::Bool
end

RelDiffUnit(name; min = -Inf, max = Inf, rel_diff, from_min = true) =
    RelDiffUnit(name, min, max, rel_diff, from_min)

function Base.show(io::IO, ::MIME"text/plain", ru::RelDiffUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Relative difference refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value:    ", ru.min)
    println(io, spad, "Maximum value:    ", ru.max)
    println(io, spad, "Maximal relative difference: ", ru.rel_diff)
    println(io, spad, "From min value:   ", ru.from_min)
end

"""
    ContourUnit(name; min=-Inf, max=Inf, contours, from_min=true)

Refine near specified contour levels (`contours`) of the target metric within
`[min, max]`. `contours` can be a number, a tuple, or a vector of levels.
"""
struct ContourUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    contours::Vector{Float64}
    from_min::Bool
end

# primary keyword constructor (vector of levels)
ContourUnit(name::Symbol; min = -Inf, max = Inf, contours::AbstractVector{<:Real}, from_min::Bool = true) =
    ContourUnit(name, min, max, Float64.(contours), from_min)

# positional conveniences that *do* dispatch:
ContourUnit(name::Symbol, contours::AbstractVector{<:Real}; min = -Inf, max = Inf, from_min::Bool = true) =
    ContourUnit(name, min, max, Float64.(contours), from_min)
ContourUnit(name::Symbol, contour::Real; min = -Inf, max = Inf, from_min::Bool = true) =
    ContourUnit(name, min, max, [Float64(contour)], from_min)
ContourUnit(name::Symbol, contours::Tuple{Vararg{Real}}; min = -Inf, max = Inf, from_min::Bool = true) =
    ContourUnit(name, min, max, Float64.(collect(contours)), from_min)

function Base.show(io::IO, ::MIME"text/plain", ru::ContourUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Contour refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value:    ", ru.min)
    println(io, spad, "Maximum value:    ", ru.max)
    println(io, spad, "Contour levels:   ", ru.contours)
    println(io, spad, "From min value:   ", ru.from_min)
end

"""
    DiffContourUnit(name; min=-Inf, max=Inf, diffs, contours, from_min=true)

Refine near combinations of absolute difference thresholds (`diffs`) and
contour levels (`contours`) within `[min, max]`.
`diffs`/`contours` accept a number, tuple, or vector.
"""
struct DiffContourUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    diffs::Vector{Float64}
    contours::Vector{Float64}
    from_min::Bool
end

# primary keyword constructor (vectors)
DiffContourUnit(name::Symbol; min = -Inf, max = Inf,
                diffs::AbstractVector{<:Real}, contours::AbstractVector{<:Real},
                from_min::Bool = true) =
    DiffContourUnit(name, min, max, Float64.(diffs), Float64.(contours), from_min)

# positional conveniences with dispatch
DiffContourUnit(name::Symbol, diffs::AbstractVector{<:Real}, contours::AbstractVector{<:Real};
                min = -Inf, max = Inf, from_min::Bool = true) =
    DiffContourUnit(name, min, max, Float64.(diffs), Float64.(contours), from_min)

DiffContourUnit(name::Symbol, diff::Real, contours::AbstractVector{<:Real};
                min = -Inf, max = Inf, from_min::Bool = true) =
    DiffContourUnit(name, min, max, [Float64(diff)], Float64.(contours), from_min)

DiffContourUnit(name::Symbol, diffs::AbstractVector{<:Real}, contour::Real;
                min = -Inf, max = Inf, from_min::Bool = true) =
    DiffContourUnit(name, min, max, Float64.(diffs), [Float64(contour)], from_min)

DiffContourUnit(name::Symbol, diffs::Tuple{Vararg{Real}}, contours::AbstractVector{<:Real};
                min = -Inf, max = Inf, from_min::Bool = true) =
    DiffContourUnit(name, min, max, Float64.(collect(diffs)), Float64.(contours), from_min)

DiffContourUnit(name::Symbol, diffs::AbstractVector{<:Real}, contours::Tuple{Vararg{Real}};
                min = -Inf, max = Inf, from_min::Bool = true) =
    DiffContourUnit(name, min, max, Float64.(diffs), Float64.(collect(contours)), from_min)

DiffContourUnit(name::Symbol, diff::Real, contour::Real;
                min = -Inf, max = Inf, from_min::Bool = true) =
    DiffContourUnit(name, min, max, [Float64(diff)], [Float64(contour)], from_min)

function Base.show(io::IO, ::MIME"text/plain", ru::DiffContourUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Difference and contour refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value:    ", ru.min)
    println(io, spad, "Maximum value:    ", ru.max)
    println(io, spad, "Maximal differences: ", ru.diffs)
    println(io, spad, "Contour levels:   ", ru.contours)
    println(io, spad, "From min value:   ", ru.from_min)
end

# == Refinement Settings Container ============================================

"""
    RefinementSettings(units...; desired_refinement_level, parallel=false, params_to_save=())

Container for adaptive‑grid refinement configuration.

Arguments
- `units...`                 : one or more refinement units (order matters).
- `desired_refinement_level` : how many refinement rounds to perform (≥ 0).
- `parallel`                 : whether the grid engine may run points in parallel.
- `params_to_save`           : tuple of metric names to extract at each point.
"""
struct RefinementSettings{T<:Tuple{Vararg{AbstractRefinementUnit}}}
    params_to_save::Tuple{Vararg{Symbol}}
    desired_refinement_level::Int64
    parallel::Bool
    units::T
end

function RefinementSettings(units...;
    desired_refinement_level::Int64,
    parallel::Bool=false,
    params_to_save::Tuple{Vararg{Symbol}}=())
    return RefinementSettings(params_to_save, desired_refinement_level, parallel, units)
end

# == Pretty Printing ===========================================================

function Base.show(io::IO, ::MIME"text/plain", ref_sets::RefinementSettings)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    iop = IOContext(io, :indent => indent + 4)

    println(io, pad,  "Grid Refinement settings:")
    println(io, spad, "Parameters to save:        ", ref_sets.params_to_save)
    println(io, spad, "Desired refinement level:  ", ref_sets.desired_refinement_level)
    println(io, spad, "Parallel computation:      ", ref_sets.parallel)

    if length(ref_sets.units) == 0
        println(io, spad, "Units: (none)")
    else
        println(io, spad, "Units:")
        for (k, u) in pairs(ref_sets.units)
            show(iop, MIME"text/plain"(), u)
        end
    end
end