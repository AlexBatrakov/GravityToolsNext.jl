# ---------------------------
# Node rules
# ---------------------------

abstract type AbstractNodeRule end

# Fallback: make the interface explicit
eval_nodes(prior::AbstractPriorSpec, rule::AbstractNodeRule)::Vector{Float64} = _notimpl(rule, eval_nodes)

# Small helper for clamping & tidying u in (0,1)
@inline function _clamp01_sort_unique!(u::Vector{Float64}; lo::Float64=PRIOR_U_EPS, hi::Float64=1-PRIOR_U_EPS)
    @inbounds for i in eachindex(u)
        u[i] = clamp(u[i], lo, hi)
    end
    sort!(u)
    unique!(u)
    return u
end

"""
    ClenshawCurtisNodes(n; eps=1e-6)

Chebyshev–Clenshaw–Curtis nodes on u ∈ (0,1) with optional trimming `eps`
to avoid exact endpoints {0,1}. These u-nodes are later mapped to parameter
values via `prior_invcdf(prior, u)`.
"""
struct ClenshawCurtisNodes <: AbstractNodeRule
    n::Int
    eps::Float64
    function ClenshawCurtisNodes(n::Integer; eps::Real=1e-6)
        n < 1 && error("ClenshawCurtisNodes: n must be ≥ 1")
        new(Int(n), Float64(eps))
    end
end

"""
    eval_nodes(prior, rule::ClenshawCurtisNodes) -> θ::Vector{Float64}

Chebyshev–Clenshaw–Curtis nodes on u∈(0,1), trimmed at `max(rule.eps, PRIOR_U_EPS)`,
then mapped to parameter space via `prior_invcdf`.
"""
function eval_nodes(prior::AbstractPriorSpec, rule::ClenshawCurtisNodes)::Vector{Float64}
    n   = rule.n
    eps = rule.eps

    # u_k = (1 + cos(π k / n)) / 2, k = 0..n
    u = Vector{Float64}(undef, n + 1)
    @inbounds for k in 0:n
        u[k + 1] = (1 + cospi(k / n)) / 2
    end

    # Clamp only the exact endpoints to avoid collapsing interior CC nodes
    lo = max(eps, PRIOR_U_EPS); hi = 1 - lo
    @inbounds for i in eachindex(u)
        ui = u[i]
        if ui == 0.0 || ui == 1.0
            u[i] = clamp(ui, lo, hi)
        end
    end
    sort!(u); unique!(u)

    return prior_invcdf(prior, u)
end

"""
    QuantileNodes(q)

Explicit quantiles `q ⊂ (0,1)`; they will be clamped into (ε, 1-ε) and
mapped by `prior_invcdf`.
"""
struct QuantileNodes <: AbstractNodeRule
    q::Vector{Float64}
    function QuantileNodes(q::AbstractVector{<:Real})
        isempty(q) && error("QuantileNodes: q is empty")
        new(collect(float.(q)))
    end
end

"""
    eval_nodes(prior, rule::QuantileNodes) -> θ::Vector{Float64}

Clamp provided quantiles to (PRIOR_U_EPS, 1-PRIOR_U_EPS) and map via `prior_invcdf`.
"""
function eval_nodes(prior::AbstractPriorSpec, rule::QuantileNodes)::Vector{Float64}
    u = copy(rule.q)
    _clamp01_sort_unique!(u)
    return prior_invcdf(prior, u)
end

"""
    ExplicitThetaNodes(theta)

Explicit parameter values; used as-is (no inverse-CDF mapping).
"""
struct ExplicitThetaNodes <: AbstractNodeRule
    theta::Vector{Float64}
    function ExplicitThetaNodes(theta::AbstractVector{<:Real})
        isempty(theta) && error("ExplicitThetaNodes: theta is empty")
        new(collect(float.(theta)))
    end
end

"""
    eval_nodes(prior, rule::ExplicitThetaNodes) -> θ::Vector{Float64}

Return a sorted, deduplicated copy of the provided θ-values (no inv-CDF mapping).
"""
function eval_nodes(::AbstractPriorSpec, rule::ExplicitThetaNodes)::Vector{Float64}
    theta = copy(rule.theta)
    sort!(theta)
    unique!(theta)
    return theta
end