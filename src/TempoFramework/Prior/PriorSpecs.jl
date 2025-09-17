# PriorSpecs.jl
# Core types and helpers for prior-based marginalization over a single parameter.

# ---------------------------
# Prior interface
# ---------------------------

abstract type AbstractPriorSpec end

# Friendly error helper
_notimpl(pr, fn) = error(string(nameof(fn)), " not implemented for ", typeof(pr))

"""
    prior_invcdf(pr, u::Real) -> Float64

Inverse-CDF mapping u ∈ (0,1) to θ. Concrete prior types must implement this scalar method.
The caller may clamp u into (0,1) if appropriate.
"""
prior_invcdf(pr::AbstractPriorSpec, u::Real) = _notimpl(pr, prior_invcdf)

"""
    prior_pdf(pr, x::Real) -> Float64

Prior density at x. Optional but recommended for diagnostics/plots.
"""
prior_pdf(pr::AbstractPriorSpec, x::Real) = _notimpl(pr, prior_pdf)

"""
    prior_logpdf(pr, x::Real) -> Float64

Log prior density at x. Optional but recommended for diagnostics/plots.
"""
prior_logpdf(pr::AbstractPriorSpec, x::Real) = _notimpl(pr, prior_logpdf)

"""
    prior_support(pr) -> (xmin::Float64, xmax::Float64)

Optional domain hint used for plotting/clamping; override when available.
"""
prior_support(pr::AbstractPriorSpec) = (NaN, NaN)

"""
    prior_cdf(pr, x::Real) -> Float64

Prior cumulative distribution function at x. Optional but recommended for diagnostics/plots.
"""
prior_cdf(pr::AbstractPriorSpec, x::Real) = _notimpl(pr, prior_cdf)

# Vectorized conveniences (built on the scalar methods)
prior_invcdf(pr::AbstractPriorSpec, us::AbstractArray{<:Real}) =
    prior_invcdf.(Ref(pr), us)

prior_pdf(pr::AbstractPriorSpec, xs::AbstractArray{<:Real}) =
    prior_pdf.(Ref(pr), xs)

prior_logpdf(pr::AbstractPriorSpec, xs::AbstractArray{<:Real}) =
    prior_logpdf.(Ref(pr), xs)

prior_cdf(pr::AbstractPriorSpec, xs::AbstractArray{<:Real}) =
    prior_cdf.(Ref(pr), xs)

"""
    prior_tail_mass(pr, lo, hi) -> Float64
Return prior mass outside [lo, hi]: CDF(lo) + (1 - CDF(hi)).
"""
prior_tail_mass(pr::AbstractPriorSpec, lo::Real, hi::Real) = begin
    clo = prior_cdf(pr, lo)
    chi = prior_cdf(pr, hi)
    return max(clo, 0.0) + max(1.0 - chi, 0.0)
end

# Niceties
prior_median(pr::AbstractPriorSpec) = prior_invcdf(pr, 0.5)
prior_quantile(pr::AbstractPriorSpec, q::Real) = prior_invcdf(pr, q)
prior_quantiles(pr::AbstractPriorSpec, qs::AbstractVector{<:Real}) = prior_invcdf(pr, qs)

# ---------------------------
# Analytic Prior
# ---------------------------

struct AnalyticPrior{D<:UnivariateDistribution} <: AbstractPriorSpec
    dist::D
end

function Base.show(io::IO, ::MIME"text/plain", pr::AnalyticPrior)
    lo, hi = prior_support(pr)
    println(io, "AnalyticPrior(", pr.dist, ")")
    println(io, "  support≈[", @sprintf("%.6g", lo), ", ", @sprintf("%.6g", hi), "]")
end

const PRIOR_U_EPS       = 1e-12
const PRIOR_SUPPORT_EPS = 1e-6

# inverse CDF (core for node mapping u -> θ)
prior_invcdf(pr::AnalyticPrior, u::Real) =
    quantile(pr.dist, clamp(float(u), PRIOR_U_EPS, 1 - PRIOR_U_EPS))

prior_pdf(   pr::AnalyticPrior, x::Real) = pdf(pr.dist, x)
prior_logpdf(pr::AnalyticPrior, x::Real) = logpdf(pr.dist, x)
prior_cdf(   pr::AnalyticPrior, x::Real) = cdf(pr.dist, x)

prior_support(pr::AnalyticPrior) = begin
    lo = try minimum(pr.dist) catch; -Inf end
    hi = try maximum(pr.dist) catch;  Inf end
    if isfinite(lo) && isfinite(hi)
        (lo, hi)
    else
        (quantile(pr.dist, PRIOR_SUPPORT_EPS),
         quantile(pr.dist, 1 - PRIOR_SUPPORT_EPS))
    end
end

# ---------------------------
# Grid Prior
# ---------------------------

"""
    GridPrior

Tabulated prior on a single parameter.

Fields
- `x`   : strictly increasing parameter grid (e.g., θ-values)
- `pdf` : prior density on `x` (non-negative; integrates to 1 w.r.t. x)
- `cdf` : cumulative prior on `x` (strictly increasing, ~[0,1])
"""
struct GridPrior <: AbstractPriorSpec
    x::Vector{Float64}
    pdf::Vector{Float64}
    cdf::Vector{Float64}
end

"""
    validate_prior(pr::GridPrior; atol_mass=1e-3) -> Bool

Validate monotonicity and normalization of a `GridPrior`.
Throws an error with a diagnostic message if invalid.
"""
function validate_prior(pr::GridPrior; atol_mass::Float64 = 1e-3)
    nx = length(pr.x)
    nx == length(pr.pdf) == length(pr.cdf) || error("GridPrior: x/pdf/cdf length mismatch")
    nx >= 2 || error("GridPrior: grid must have at least 2 points")

    # x strictly increasing
    @inbounds for i in 2:nx
        pr.x[i] > pr.x[i-1] || error("GridPrior: x must be strictly increasing")
    end

    # pdf non-negative
    @inbounds for v in pr.pdf
        v >= 0.0 || error("GridPrior: pdf must be non-negative")
    end

    # cdf strictly increasing and within [0,1] (allow small tolerances)
    @inbounds for i in 2:nx
        pr.cdf[i] > pr.cdf[i-1] || error("GridPrior: cdf must be strictly increasing")
    end
    c0, c1 = pr.cdf[1], pr.cdf[end]
    (c0 >= -1e-8 && c0 <= 1e-3)      || error("GridPrior: cdf[1]≈0 expected, got $c0")
    (c1 <= 1 + 1e-8 && c1 >= 1 - 1e-3) || error("GridPrior: cdf[end]≈1 expected, got $c1")

    # total mass
    mass = c1 - c0
    abs(mass - 1.0) <= atol_mass || error("GridPrior: total mass $(mass) deviates from 1 by > $(atol_mass)")

    return true
end

function Base.show(io::IO, ::MIME"text/plain", pr::GridPrior)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    n   = length(pr.x)
    xmin, xmax = pr.x[1], pr.x[end]
    c0, c1 = pr.cdf[1], pr.cdf[end]
    mass   = c1 - c0
    pdfmin = minimum(pr.pdf)
    pdfmax = maximum(pr.pdf)

    println(io, pad, "GridPrior")
    println(io, spad, "n:        ", n)
    println(io, spad, @sprintf("x[min]:  %.6g", xmin))
    println(io, spad, @sprintf("x[max]:  %.6g", xmax))
    println(io, spad, @sprintf("cdf[1]:  %.3g    cdf[end]: %.3g    mass: %.6f", c0, c1, mass))
    println(io, spad, @sprintf("pdf[min]: %.3g    pdf[max]: %.3g", pdfmin, pdfmax))
end

# ---------- small linear interpolation helpers (endpoint-clamped) ----------

@inline function _lininterp(xg::AbstractVector{<:Real}, yg::AbstractVector{<:Real}, x::Real)
    # assumes xg strictly increasing; returns endpoint values outside support
    j = searchsortedlast(xg, x)
    if j <= 0
        return yg[1]
    elseif j >= length(xg)
        return yg[end]
    else
        x0 = xg[j];  x1 = xg[j+1]
        y0 = yg[j];  y1 = yg[j+1]
        t = (float(x) - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)
    end
end

_lininterp(xg, yg, xs::AbstractArray{<:Real}) = _lininterp.(Ref(xg), Ref(yg), xs)

# ---------- interface implementations for GridPrior ----------

"""
    prior_invcdf(pr::GridPrior, u::Real) -> Float64

Invert the tabulated CDF with linear interpolation.
`u` is clamped into `[cdf[1], cdf[end]]` before inversion.
"""
function prior_invcdf(pr::GridPrior, u::Real)::Float64
    ui = clamp(float(u), pr.cdf[1], pr.cdf[end])
    # find first j with cdf[j] ≥ ui
    j = searchsortedfirst(pr.cdf, ui)
    if j <= 1
        return pr.x[1]
    elseif j > length(pr.x)
        return pr.x[end]
    else
        c0, c1 = pr.cdf[j-1], pr.cdf[j]
        x0, x1 = pr.x[j-1],  pr.x[j]
        return x0 + (ui - c0) * (x1 - x0) / (c1 - c0)
    end
end

function prior_pdf(pr::GridPrior, x::Real)
    xv = pr.x
    if x < xv[1] || x > xv[end]
        return 0.0
    else
        return _lininterp(pr.x, pr.pdf, x)
    end
end

prior_logpdf(pr::GridPrior, x::Real) =
    (x < pr.x[1] || x > pr.x[end]) ? -Inf : log(prior_pdf(pr, x))

"""
    prior_cdf(pr::GridPrior, x::Real) -> Float64

Evaluate the tabulated CDF with linear interpolation.
Returns `cdf[1]` if `x < x[1]` and `cdf[end]` if `x > x[end]`.
"""
function prior_cdf(pr::GridPrior, x::Real)
    if x < pr.x[1]
        return pr.cdf[1]
    elseif x > pr.x[end]
        return pr.cdf[end]
    else
        return _lininterp(pr.x, pr.cdf, x)
    end
end

prior_support(pr::GridPrior) = (pr.x[1], pr.x[end])

# ---------- convenience constructor: GridPrior(x, pdf) ----------

"""
    GridPrior(x::AbstractVector{<:Real}, pdf::AbstractVector{<:Real})

Build a normalized `GridPrior` from tabulated `x` and (possibly unnormalized) `pdf`.
- Sorts by `x` if needed;
- Clamps negative `pdf` to 0;
- Computes CDF via cumulative trapezoid;
- Normalizes both `pdf` and `cdf` to unit mass;
- Validates the result.
"""
function GridPrior(x::AbstractVector{<:Real}, pdf::AbstractVector{<:Real})
    xv = collect(float.(x))
    pv = collect(float.(pdf))
    length(xv) == length(pv) || error("GridPrior: x/pdf length mismatch")
    length(xv) >= 2 || error("GridPrior: need at least 2 grid points")

    # sort by x if unsorted
    if !issorted(xv)
        p = sortperm(xv)
        xv = xv[p]; pv = pv[p]
    end

    # enforce non-negativity
    @. pv = max(pv, 0.0)

    # cumulative trapezoid for cdf
    cdf = similar(xv)
    cdf[1] = 0.0
    s = 0.0
    @inbounds for i in 2:length(xv)
        dx = xv[i] - xv[i-1]
        dx > 0 || error("GridPrior: x must be strictly increasing")
        s += 0.5 * dx * (pv[i] + pv[i-1])
        cdf[i] = s
    end
    s > 0 || error("GridPrior: pdf integrates to zero")

    # normalize mass to 1
    @. pv  = pv  / s
    @. cdf = cdf / s

    gp = GridPrior(xv, pv, cdf)
    validate_prior(gp)
    return gp
end

# ---------------------------
# Sampled Prior (KDE-backed → GridPrior with lazy cache)
# ---------------------------

struct SampledPrior <: AbstractPriorSpec
    samples::Vector{Float64}
    bandwidth::Union{Nothing,Float64}
    _grid::Base.RefValue{Union{Nothing,GridPrior}}
    function SampledPrior(samples::AbstractVector{<:Real};
                          bandwidth::Union{Nothing,Real}=nothing)
        s = Float64[si for si in samples if isfinite(si)]
        length(s) >= 2 || error("SampledPrior: need at least 2 finite samples")
        return new(
            s,
            bandwidth === nothing ? nothing : Float64(bandwidth),
            Ref{Union{Nothing,GridPrior}}(nothing)  # ← typed Ref that can later hold a GridPrior
        )
    end
end

# Pretty-print SampledPrior with indentation and optional GridPrior details
function Base.show(io::IO, ::MIME"text/plain", pr::SampledPrior)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "SampledPrior")
    println(io, spad, "n:         ", length(pr.samples))
    if pr.bandwidth === nothing
        println(io, spad, "bandwidth: (auto)")
    else
        println(io, spad, "bandwidth: ", pr.bandwidth)
    end

    if pr._grid[] === nothing
        println(io, spad, "materialized: false")
    else
        println(io, spad, "materialized: true")
        # Show the cached GridPrior with deeper indentation
        show(IOContext(io, :indent => indent + 4), MIME"text/plain"(), pr._grid[]::GridPrior)
    end
end

"""
    materialize(pr::SampledPrior; n_max=4096) -> GridPrior

Build (and cache) a `GridPrior` via KDE of the samples. Subsequent calls reuse the cache.
`n_max` is an upper bound on the KDE output grid size (the actual size is controlled by KernelDensity.jl).
"""
function materialize(pr::SampledPrior; n_max::Int=4096)::GridPrior
    g = pr._grid[]
    if g !== nothing
        return g
    end

    # --- KDE ---
    # Requires KernelDensity.jl in your project:
    # using KernelDensity
    kd = pr.bandwidth === nothing ? KernelDensity.kde(pr.samples) :
                                    KernelDensity.kde(pr.samples; bandwidth=pr.bandwidth)

    x   = collect(float.(kd.x))
    pdf = collect(float.(kd.density))

    # Optionally subsample very long tails to keep the grid modest
    if length(x) > n_max
        idxf = range(1, length(x), length=n_max)              # Float indices
        p    = clamp.(round.(Int, idxf), 1, length(x))        # safe Ints
        unique!(p)                                            # just in case
        x   = x[p]
        pdf = pdf[p]
    end

    # Build tabulated prior (normalization and cdf are handled inside)
    gp = GridPrior(x, pdf)

    # Cache
    pr._grid[] = gp
    return gp
end

# Delegate interface to the materialized grid
prior_invcdf(pr::SampledPrior, u::Real)    = prior_invcdf(materialize(pr), u)
prior_pdf(   pr::SampledPrior, x::Real)    = prior_pdf(materialize(pr), x)
prior_logpdf(pr::SampledPrior, x::Real)    = prior_logpdf(materialize(pr), x)
prior_cdf(   pr::SampledPrior, x::Real)    = prior_cdf(materialize(pr), x)
prior_support(pr::SampledPrior)            = prior_support(materialize(pr))


# ---------------------------
# I/O helpers for priors
# ---------------------------

"""
    read_prior_samples(path; column=1, comment="#", delim=nothing, transform=identity) -> Vector{Float64}

Read numeric samples from a text file to build an empirical prior.

- Skips empty lines and lines starting with `comment`.
- Splits each line by `delim` (default: regex for commas/whitespace).
- `column` can be an `Int` (1-based) or `:auto` to pick the first parseable number on each line.
- `transform` is applied to each parsed value (e.g., to convert units or apply a sign flip).

Returns a vector of finite Float64 values. Errors if fewer than 2 values are found.
"""
function read_prior_samples(path::AbstractString;
    column::Union{Int,Symbol}=1,
    comment::AbstractString="#",
    delim=nothing,
    transform::Function=identity,
)::Vector{Float64}
    isfile(path) || error("read_prior_samples: no such file: $path")
    vals = Float64[]
    open(path, "r") do io
        for ln in eachline(io)
            s = strip(ln)
            isempty(s) && continue
            startswith(s, comment) && continue
            toks = split(s, delim === nothing ? r"[,\s]+" : delim; keepempty=false)
            isempty(toks) && continue

            v::Union{Nothing,Float64} = nothing
            if column === :auto
                idx = findfirst(t -> tryparse(Float64, t) !== nothing, toks)
                idx === nothing && continue
                v = parse(Float64, toks[idx])
            else
                1 <= column <= length(toks) || continue
                p = tryparse(Float64, toks[column])
                p === nothing && continue
                v = p
            end

            v = transform(v)
            isfinite(v) && push!(vals, float(v))
        end
    end
    length(vals) > 1 || error("read_prior_samples: collected $(length(vals)) finite values from $path")
    return vals
end

"""
    SampledPrior(path::AbstractString;
                 column::Union{Int,Symbol}=1,
                 comment::AbstractString="#",
                 delim=nothing,
                 transform::Function=identity,
                 bandwidth::Union{Nothing,Real}=nothing) -> SampledPrior

Convenience constructor: read numeric samples from `path` and build a `SampledPrior`.
- `column`: 1-based column index or `:auto` (first parseable number per line).
- `comment`: lines starting with this string are skipped.
- `delim`: custom delimiter for `split` (default: regex for commas/whitespace).
- `transform`: applied to each parsed value (e.g. unit conversion).
- `bandwidth`: optional KDE bandwidth override used at materialization time.
"""
function SampledPrior(path::AbstractString;
    column::Union{Int,Symbol}=1,
    comment::AbstractString="#",
    delim=nothing,
    transform::Function=identity,
    bandwidth::Union{Nothing,Real}=nothing,
)
    s = read_prior_samples(path; column=column, comment=comment, delim=delim, transform=transform)
    return SampledPrior(s; bandwidth=bandwidth)
end