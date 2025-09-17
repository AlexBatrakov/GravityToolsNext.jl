# src/AdaptiveGridFramework/GridAxes.jl
# -----------------------------------------------------------------------------
# 1D Grid axis abstraction used by adaptive grid tasks
#
# This file defines:
#   * `AbstractGridRule` — a discretization rule (linear/log/explicit)
#   * `GridAxis`        — a typed axis with name, range, resolution and rule
#   * linspace(ax)      — linear parameterization (rule space)
#   * values(ax)        — physical coordinates corresponding to linspace(ax)
#
# Design notes
# - We keep a single lightweight type per axis. Concrete rule is encoded in the
#   type parameter `R`, so dispatch stays fast and allocations are minimal.
# - For `ExplicitRule` the axis stores user-provided values verbatim; `min/max/N`
#   are derived at construction time and used only for pretty-printing.
# - `LogRule()` (no args) means "auto sign": the sign is inferred from `min` at
#   the time of sampling. This is convenient for negative-only logarithmic axes.
# -----------------------------------------------------------------------------

# == Rules ====================================================================

"""
    AbstractGridRule

A discretization rule for a 1D grid axis. Implementations:
- [`LinRule`](@ref) — uniform spacing in the physical domain.
- [`LogRule`](@ref) — uniform spacing in log10(|x|); sign can be auto-detected.
- [`ExplicitRule`](@ref) — use a user-provided set of values verbatim.
"""
abstract type AbstractGridRule end

"""
    LinRule <: AbstractGridRule

Linear spacing rule in the *physical* coordinate: `values(ax)` is an evenly
spaced `LinRange(min, max, N)`.
"""
struct LinRule    <: AbstractGridRule end

"""
    LogRule(sign)
    LogRule()

Logarithmic spacing rule: uniformly spaced in `log10(|x|)`.

Arguments
- `sign::Int` : `+1` or `-1`. Controls the sign of returned physical values.
                Use `LogRule()` (no args) to enable *auto sign* — inferred from
                the axis `min` at sampling time.

Notes
- Physical values are computed as sign * 10 .^ linspace(ax), where linspace(ax) is uniform in log10(|x|).
"""
struct LogRule    <: AbstractGridRule
    sign::Int  # +1 / -1; 0 => auto from axis min
end

# Convenience constructor: auto-sign (0 means choose sign from `ax.min` later)
LogRule() = LogRule(0)

"""
    ExplicitRule(vals)

Explicit set of coordinates for the axis. The rule carries the values; the axis
keeps `min/max/N` for convenience.
"""
struct ExplicitRule <: AbstractGridRule
    vals::Vector{Float64}
end

# == Rule-specific sampling via multiple dispatch =============================

"""
    _rule_linvals(rule::LinRule, min::Float64, max::Float64, N::Int)

Return linear parameterization for a linear axis.
"""
_rule_linvals(::LinRule, min::Float64, max::Float64, N::Int) = collect(LinRange(min, max, N))

"""
    _rule_values(rule::LinRule, lin::Vector{Float64}; min_hint::Float64=NaN)

For a linear axis, physical values coincide with the linear parameterization.
`min_hint` is accepted for a uniform API and ignored.
"""
_rule_values(::LinRule, lin::Vector{Float64}; min_hint::Float64=NaN) = lin

"""
    _rule_linvals(rule::LogRule, min::Float64, max::Float64, N::Int)

Return uniformly spaced samples in log10(|x|) for a logarithmic axis.
Throws if min/max are zero or have different signs.
"""
function _rule_linvals(rule::LogRule, min::Float64, max::Float64, N::Int)
    (min == 0 || max == 0) && error("LogRule: min/max must be nonzero")
    sign(min) != sign(max) && error("LogRule: min and max must have the same sign")
    lo = log10(abs(min)); hi = log10(abs(max))
    return collect(LinRange(lo, hi, N))
end

"""
    _rule_values(rule::LogRule, lin::Vector{Float64}; min_hint::Float64=1.0)

Map logarithmic linear samples back to physical space, using the rule sign if
provided or the sign inferred from `min_hint` (axis min) for auto-sign.
"""
function _rule_values(rule::LogRule, lin::Vector{Float64}; min_hint::Float64=1.0)
    s = rule.sign == 0 ? sign(min_hint) : rule.sign
    return s .* (10.0 .^ lin)
end

"""
    _rule_linvals(rule::ExplicitRule, min::Float64, max::Float64, N::Int)

For an explicit axis, the linear parameterization equals the stored values.
"""
_rule_linvals(rule::ExplicitRule, ::Float64, ::Float64, ::Int) = copy(rule.vals)

"""
    _rule_values(rule::ExplicitRule, lin::Vector{Float64}; min_hint::Float64=NaN)

For an explicit axis, physical values equal the stored values.
`min_hint` is accepted for a uniform API and ignored.
"""
_rule_values(::ExplicitRule, lin::Vector{Float64}; min_hint::Float64=NaN) = lin

# == Axis =====================================================================

"""
    GridAxis(name; min, max, N, rule)
    GridAxis(name, values::AbstractVector)

A typed grid axis with a chosen discretization rule. Call [`values`](@ref) to
obtain physical coordinates and [`linspace`](@ref) to obtain the rule's linear
parameterization.

Arguments
- `name::Symbol` : semantic axis name (e.g. `:PBDOT`, `:beta0`).
- `min, max`     : axis bounds (ignored for `ExplicitRule`).
- `N::Integer`   : number of grid points (must be ≥ 1).
- `rule`         : one of `LinRule()`, `LogRule([sign])`, or `ExplicitRule`.
- values       : alternative constructor; makes an ExplicitRule from a supplied vector of nodes.

Examples
```julia
ax1 = GridAxis(:x; min=-1.0, max=1.0, N=21, rule=LinRule())
ax2 = GridAxis(:m; min=-1e-4, max=-1e-2, N=16, rule=LogRule())  # auto-sign (−)
ax3 = GridAxis(:θ, [-1.0, -0.5, 0.0, 0.5, 1.0])                 # explicit nodes
```
"""
struct GridAxis{R<:AbstractGridRule}
    name::Symbol
    min::Float64
    max::Float64
    N::Int
    rule::R
    lin_values::Vector{Float64}
    values::Vector{Float64}
    function GridAxis{R}(name::Symbol, min::Real, max::Real, N::Integer, rule::R) where {R<:AbstractGridRule}
        N ≤ 0 && error("GridAxis: N must be ≥ 1, got $N")
        minf = float(min); maxf = float(max); nint = Int(N)

        # Compute rule-space samples and physical values via dispatch
        lin_vals = _rule_linvals(rule, minf, maxf, nint)
        vals     = _rule_values(rule, lin_vals; min_hint=minf)

        return new{R}(name, minf, maxf, nint, rule, lin_vals, vals)
    end
end

# -- Convenience constructors --------------------------------------------------

"""
    GridAxis(name::Symbol; min, max, N, rule::AbstractGridRule=LinRule())
Construct an axis on `[min, max]` with `N` nodes, using the provided `rule`.
Supported rules: `LinRule()`, `LogRule([sign])`, or `ExplicitRule(vals)`
(via the positional constructor `GridAxis(name, vals)`).
"""
GridAxis(name::Symbol; min, max, N, rule::AbstractGridRule=LinRule()) =
    GridAxis{typeof(rule)}(name, min, max, N, rule)

"""
    GridAxis(name::Symbol, vals::AbstractVector)
Construct an explicit axis that uses the provided `vals` verbatim.
"""
GridAxis(name::Symbol, vals::AbstractVector{<:Real}) =
    GridAxis{ExplicitRule}(name, first(vals), last(vals), length(vals), ExplicitRule(collect(float.(vals))))

# -- Aliases for readability ---------------------------------------------------
const LinAxis      = GridAxis{LinRule}
const LogAxis      = GridAxis{LogRule}
const ExplicitAxis = GridAxis{ExplicitRule}

# == Sampling helpers ==========================================================

"""
    linspace(ax)

Return the *linear* parameterization for an axis:
- For `LinAxis`, this equals the physical coordinates.
- For `LogAxis`, this is uniform in `log10(|x|)`.
- For `ExplicitAxis`, returns a copy of the stored values.

The result is a copy (safe to mutate by the caller).
"""
linspace(ax::GridAxis) = copy(ax.lin_values)

"""
    values(ax)

Return the *physical* coordinates for an axis. For `LogAxis`, uses the rule
sign (`+1/−1`) or infers it from `ax.min` when the rule is constructed with
`LogRule()`.

The result is a copy (safe to mutate by the caller).
"""
values(ax::GridAxis) = copy(ax.values)

"""
    axisvalues(ax)
Alias for [`values(ax)`](@ref). Provided to avoid accidental name shadowing
in user code that defines local variables named `values`.
"""
axisvalues(ax::GridAxis) = values(ax)

# == Printing =================================================================

function Base.show(io::IO, ax::GridAxis)
    print(io, string(ax.name), " axis: ",
          "N=", ax.N, ", range=(", ax.min, ", ", ax.max, "), rule=", nameof(typeof(ax.rule)),
          ", cached samples")
end

# == Refinement ================================================================

"""
    refine(ax::GridAxis) -> GridAxis

Refine a 1D axis by inserting midpoints to reach `2N-1` nodes.

- `LinAxis`  : midpoints in physical space on [min, max].
- `LogAxis`  : midpoints uniform in log10(|x|) (geometric means in physical space).
- `ExplicitAxis` : midpoints in physical space between consecutive user nodes.

Returns a **new** axis; the original is unchanged.
"""
function refine(ax::LinAxis)
    GridAxis(ax.name; min=ax.min, max=ax.max, N=2*ax.N - 1, rule=LinRule())
end

function refine(ax::LogAxis)
    # Preserve the rule’s sign/auto-sign choice
    GridAxis(ax.name; min=ax.min, max=ax.max, N=2*ax.N - 1, rule=ax.rule)
end

function refine(ax::ExplicitAxis)
    v = values(ax)
    n = length(v)
    n <= 1 && return ax  # nothing to refine

    # Ensure monotonicity to avoid self-intersections
    if !(issorted(v; lt = <) || issorted(v; lt = >))
        error("ExplicitAxis refine: values must be strictly monotonic to insert midpoints")
    end

    w = Vector{Float64}(undef, 2n - 1)
    for i in 1:n-1
        w[2i - 1] = v[i]
        w[2i]     = 0.5 * (v[i] + v[i+1])
    end
    w[end] = v[end]

    # Rebuild as explicit axis; min/max/N will update automatically.
    return GridAxis(ax.name, w)
end

"""
    refine(ax::GridAxis, levels::Integer) -> GridAxis

Apply `refine` repeatedly `levels` times.
"""
function refine(ax::GridAxis, levels::Integer)
    levels < 0 && error("levels must be ≥ 0")
    out = ax
    for _ in 1:levels
        out = refine(out)
    end
    return out
end
