# PriorMarginalizedTempoTask.jl
# This file depends on:


"""
    _resolve_ref_theta(s, theta, chi2; chi2_spline=nothing) -> Float64

Choose the reference θ according to `s.ref_strategy`.

- `:custom_value`   → the provided `ref_value` (clamped to prior support if finite).
- `:prior_median`   → prior median (does not have to coincide with a node).
- `:grid_argmin`    → θ at the discrete argmin of `chi2` on the provided nodes.
- `:spline_argmin`  → θ minimizing a smooth χ²(θ) proxy on [min(theta), max(theta)].
                      If `chi2_spline` is provided (callable f(θ)), it is used.
                      Otherwise a light fallback is used (dense linear interp sampling).

Notes:
- The function returns the θ itself. Оценку χ²(θ_ref) ты потом получаешь либо
  на узле (если совпало), либо интерполяцией/сплайном — это уже следующий шаг.
"""
function _resolve_ref_theta(
    s::PriorMarginalizationSettings,
    theta::AbstractVector{<:Real},
    chi2::AbstractVector{<:Real};
    chi2_spline::Union{Nothing,Function}=nothing,
)::Float64
    length(theta) == length(chi2) || error("theta/chi2 length mismatch")

    # ensure ascending θ (keep chi2 in sync)
    if !issorted(theta)
        p = sortperm(theta)
        θ = collect(float.(theta[p]))
        c = collect(float.(chi2[p]))
    else
        θ = collect(float.(theta))
        c = collect(float.(chi2))
    end

    rs = s.ref_strategy
    if rs === :custom_value
        v = s.ref_value
        v === nothing && error("ref_strategy=:custom_value requires settings.ref_value")
        lo, hi = prior_support(s.prior)
        return (isfinite(lo) && isfinite(hi)) ? clamp(v::Float64, lo, hi) : (v::Float64)

    elseif rs === :prior_median
        return prior_median(s.prior)

    elseif rs === :grid_argmin
        i = argmin(c)
        return θ[i]

    elseif rs === :spline_argmin
        # 1) get a smooth evaluator f(θ)
        f = if chi2_spline !== nothing
            chi2_spline
        else
            # fallback: linear interpolation with dense sampling
            # (без внешних зависимостей и довольно устойчиво)
            (x -> begin
                # piecewise-linear χ² at x
                j = searchsortedlast(θ, x)
                if j <= 0
                    c[1]
                elseif j >= length(θ)
                    c[end]
                else
                    x0=θ[j]; x1=θ[j+1]; y0=c[j]; y1=c[j+1]
                    t = (x - x0) / (x1 - x0)
                    y0 + t*(y1 - y0)
                end
            end)
        end

        # 2) crude 1D minimization on [θmin, θmax] by dense sampling
        a = θ[1]; b = θ[end]
        nprobe = max(64, 5*length(θ))
        xs = range(a, b; length=nprobe)
        ys = similar(xs, Float64)
        @inbounds for i in eachindex(xs)
            ys[i] = f(xs[i])
        end
        return xs[argmin(ys)]

    else
        error("Unknown ref_strategy: $rs")
    end
end


# ---------------------------
# Task wrapper
# ---------------------------

"""
    PriorMarginalizedTempoTask{T,P,N} <: SingleTempoTask

Wraps a base task (e.g., `BasicTempoTask`) and runs it across a set of prior nodes,
then marginalizes the chosen likelihood metric over the prior.
"""
struct PriorMarginalizedTempoTask{T<:SingleTempoTask,P<:AbstractPriorSpec,N<:AbstractNodeRule} <: SingleTempoTask
    unit_task::T
    settings::PriorMarginalizationSettings{P,N}
end

function Base.show(io::IO, ::MIME"text/plain", t::PriorMarginalizedTempoTask)
    println(io, "PriorMarginalizedTempoTask")
    println(io, "  unit_task:")
    show(IOContext(io, :indent => 4), MIME"text/plain"(), t.unit_task)
    println(io, "  settings:")
    show(IOContext(io, :indent => 4), MIME"text/plain"(), t.settings)
end

# ---------------------------
# Small utilities
# ---------------------------

# nonuniform trapezoid ∫ f(θ) dθ on sorted θ
@inline function _trapz(theta::AbstractVector{<:Real}, f::AbstractVector{<:Real})
    n = length(theta); n == length(f) || error("trapz: length mismatch")
    n >= 2 || return 0.0
    acc = 0.0
    @inbounds for i in 1:n-1
        dθ = theta[i+1] - theta[i]
        acc += 0.5 * dθ * (f[i] + f[i+1])
    end
    return acc
end

@inline function _nearest_index(x::Real, grid::AbstractVector{<:Real})
    j = searchsortedfirst(grid, x)
    if j <= 1
        return 1
    elseif j > length(grid)
        return length(grid)
    else
        # choose closer of j and j-1
        abs(grid[j] - x) < abs(x - grid[j-1]) ? j : (j-1)
    end
end

# ---------------------------
# Core runner
# ---------------------------

"""
    run_task(task::PriorMarginalizedTempoTask) -> GeneralTempoResult

Algorithm (high-level):
1) Build θ-nodes from `(settings.prior, settings.nodes)`.
2) For each node, call `settings.make_task(θ, i, unit_task)` and `run_task` it.
3) Collect metric values (e.g., χ²) from each node result.
4) Choose a reference χ² by `settings.ref_strategy` and form Δχ² = χ² - χ²_ref.
5) Compute the marginal likelihood integral ∫ prior(θ) * exp(-Δχ²/2) dθ (trapezoid on θ-nodes).
6) Report the corrected χ²: χ²_ref - 2 log(integral).
7) Return a top-level `GeneralTempoResult` whose `final` is the representative node's final
   iteration; include all node results as `subresults` if requested; push key metrics.
"""
function run_task(task::PriorMarginalizedTempoTask)::GeneralTempoResult
    s = task.settings

    # 1) nodes in parameter space
    θ_nodes = eval_nodes(s.prior, s.nodes)
    length(θ_nodes) >= 2 || error("Prior marginalization needs at least 2 nodes")
    # ensure sorted (eval_nodes already sorts, but be defensive)
    sort!(θ_nodes)
    unique!(θ_nodes)

    # 2) run per-node tasks
    node_results = Vector{GeneralTempoResult}(undef, length(θ_nodes))
    metrics_vals = Vector{Float64}(undef, length(θ_nodes))

    for (i, θ) in pairs(θ_nodes)
        node_task = s.make_task(θ, i, task.unit_task)
        res_i = run_task(node_task)
        node_results[i] = res_i
        metrics_vals[i] = get(res_i.metrics, s.likelihood_source, NaN)
    end

    # Drop nodes with NaN metric (failed or missing); keep them in subresults but exclude from integral
    valid = findall(isfinite, metrics_vals)
    if length(valid) < 2
        error("Not enough valid node metrics to marginalize (have $(length(valid))).")
    end
    θv = θ_nodes[valid]
    χ²v = metrics_vals[valid]

    # 3) choose reference χ²
    ref_idx = begin
        if s.ref_strategy == :prior_median
            _nearest_index(prior_median(s.prior), θv)
        elseif s.ref_strategy == :node_min
            argmin(χ²v)
        elseif s.ref_strategy == :node_median
            _nearest_index(median(θv), θv)
        else
            error("Unknown ref_strategy=$(s.ref_strategy)")
        end
    end
    χ²_ref = χ²v[ref_idx]
    θ_ref  = θv[ref_idx]

    # 4) marginalization integral over θ (nonuniform trapezoid on nodes)
    Δχ² = χ²v .- χ²_ref
    L   = @. exp(-0.5 * Δχ²)                     # likelihood up to the reference constant
    pθ  = prior_pdf(s.prior, θv)                 # prior density at nodes
    integrand = @. pθ * L

    Z = _trapz(θv, integrand)                    # ∫ p(θ) L(θ) dθ  (units of θ)
    if !(Z > 0 && isfinite(Z))
        error("Marginalization integral is non-positive or invalid: Z=$Z")
    end

    χ²_marg = χ²_ref - 2 * log(Z)                # corrected chi^2

    # 5) choose representative node for the top-level "final"
    rep_idx = begin
        if s.representative == :prior_median
            _nearest_index(prior_median(s.prior), θ_nodes)
        elseif s.representative == :node_min
            argmin(metrics_vals)
        elseif s.representative == :first
            1
        elseif s.representative == :last
            length(θ_nodes)
        elseif s.representative == :index
            clamp(s.rep_index, 1, length(θ_nodes))
        else
            error("Unknown representative=$(s.representative)")
        end
    end

    rep_res = node_results[rep_idx]
    rep_it  = rep_res.final  # InternalIterationResult

    # 6) build top-level result (use representative final iteration as the "final")
    meta = Dict{Symbol,Any}(
        :parameter          => s.parameter,
        :theta_nodes        => θ_nodes,
        :theta_valid        => θv,
        :metric_name        => s.likelihood_source,
        :metric_all         => metrics_vals,
        :ref_strategy       => s.ref_strategy,
        :ref_idx_valid      => ref_idx,
        :ref_theta          => θ_ref,
        :integral_method    => :trapezoid,
        :Z                  => Z,
        :representative     => s.representative,
        :rep_idx            => rep_idx,
        :keep_node_results  => s.keep_node_results,
    )

    subresults = s.keep_node_results ? node_results : GeneralTempoResult[]
    top = GeneralTempoResult(
        rep_it;
        subresults      = subresults,
        subresult_type  = :prior_node,
        metadata        = meta,
    )

    # 7) attach compact scalar metrics at the top level
    top.metrics[:chi2_ref]          = χ²_ref
    top.metrics[:chi2_marginalized] = χ²_marg
    top.metrics[:theta_ref]         = θ_ref             # note: Float64, OK in metrics
    top.metrics[:nodes_valid]       = length(valid)

    if s.metrics_hook !== nothing
        # allow user to add custom metrics (must be Pair{Symbol,Float64})
        for (k,v) in s.metrics_hook(top, θ_ref)
            top.metrics[k] = v
        end
    end

    return top
end