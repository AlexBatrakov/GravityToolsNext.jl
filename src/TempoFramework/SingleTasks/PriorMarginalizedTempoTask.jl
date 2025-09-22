# src/TempoFramework/SingleTasks/PriorMarginalizedTempoTask.jl
# A prior-marginalized wrapper around a `SingleTempoTask`.
# It runs the base task across a grid of θ values and aggregates results using a given prior.

# --------------------------------------------------------------------------------------------------------------
# Task wrapper
# --------------------------------------------------------------------------------------------------------------

"""
    PriorMarginalizedTempoTask{T,P,N} <: SingleTempoTask

A task wrapper that evaluates a base single-tempo task over a grid of prior values `θ`
and aggregates node results into one posterior-aware result.

Type parameters
- `T` — a concrete `SingleTempoTask` (e.g. `BasicTempoTask`)
- `P` — a prior specification implementing the `AbstractPriorSpec` API
- `N` — a node rule implementing the `AbstractNodeRule` API

Fields
- `base_task`  — the underlying single task to run at each grid node
- `settings`   — prior definition, node grid rule, execution options, etc.
"""
struct PriorMarginalizedTempoTask{T<:SingleTempoTask,P<:AbstractPriorSpec,N<:AbstractNodeRule} <: SingleTempoTask
    base_task::T
    settings::PriorMarginalizationSettings{P,N}
end

function Base.show(io::IO, ::MIME"text/plain", task::PriorMarginalizedTempoTask)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    iop    = IOContext(io, :indent => indent + 4)

    println(io, pad,  "PriorMarginalizedTempoTask")
    println(io, spad, "base_task:")
    show(iop, MIME"text/plain"(), task.base_task)
    println(io, spad, "settings:")
    show(iop, MIME"text/plain"(), task.settings)
end

# --------------------------------------------------------------------------------------------------------------
# Internal node result (helper)
# --------------------------------------------------------------------------------------------------------------

"""
    _PriorNodeResult

Internal container holding data for a single grid node.

Fields
- `index`         — 1-based node position on the grid
- `theta`         — θ value at the node
- `workdir`       — absolute path to the node's working directory
- `like_metric`   — scalar metric used as a likelihood proxy (e.g. `metrics[:chi2_fit]`)
- `prior_density` — prior pdf value π(θ) at this node
- `result`        — node's `GeneralTempoResult`
"""
struct _PriorNodeResult
    index::Int
    theta::Float64
    workdir::String
    like_metric::Float64
    prior_density::Float64
    result::GeneralTempoResult
end


"Convenience constructor that computes `like_metric` and `prior_density`."
_PriorNodeResult(index::Int, theta::Float64, workdir::AbstractString,
                 res::GeneralTempoResult, like_key::Symbol, prior) =
    _PriorNodeResult(
        index,
        theta,
        String(workdir),
        result_metric(res, like_key),
        prior_pdf(prior, theta),
        res,
    )

# ---------------------------
# Pure helpers (no I/O)
# ---------------------------
"""
    _build_thetas(s) -> Vector{Float64}

Construct the θ grid according to `s.nodes` using the given prior.
"""
_build_thetas(s::PriorMarginalizationSettings) = eval_nodes(s.prior, s.nodes)

# mref: baseline used to shift χ²-like sources (nodes-only baseline; no spline/interp)
_resolve_mref(::PriorMarginalizationSettings, thetas::Vector{Float64}, metrics::Vector{Float64}) =
    minimum(metrics)  # nodes-only baseline; no spline/interp

"""
    _resolve_reference(s, thetas, metrics) -> (θref, mref)

Choose a reference θ and its metric value used for χ² shifting and reporting.
Supported strategies:
- `:grid_argmin`     — θ at argmin of node metrics
- `:spline_argmin`   — currently uses node argmin (spline argmin may be added later)
- `:prior_median`    — θ at prior median (nearest node)
- `:custom_value`    — user-provided `s.ref_value`
"""
function _resolve_reference(s::PriorMarginalizationSettings,
                            thetas::Vector{Float64},
                            metrics::Vector{Float64})
    if s.ref_strategy === :grid_argmin || s.ref_strategy === :spline_argmin
        i = argmin(metrics)
        return (thetas[i], metrics[i])
    elseif s.ref_strategy === :prior_median
        θref = prior_median(s.prior)
        i = findmin(abs.(thetas .- θref))[2]    # nearest node to median
        return (θref, metrics[i])
    elseif s.ref_strategy === :custom_value
        θref = s.ref_value::Float64
        i = findmin(abs.(thetas .- θref))[2]
        return (θref, metrics[i])
    else
        error("unknown ref_strategy $(s.ref_strategy)")
    end
end

# Convert a scalar metric into a log-likelihood according to `source`.
# For χ²-like sources we use: logL = -0.5*(χ² - mref).
function _as_loglike(source::Symbol, val::Float64, mref::Float64=0.0)::Float64
    if startswith(String(source), "chi2")
        return -0.5 * (val - mref)
    elseif source === :loglike
        return val
    elseif source === :negloglike
        return -val
    else
        error("unsupported likelihood_source = $source")
    end
end

"""
    _posterior_weights(thetas, loglikes, prior) -> Vector{Float64}

Compute normalized posterior weights on the node grid:
w_i ∝ exp(logL_i) * π(θ_i), normalized with log-sum-exp stabilization.
"""
function _posterior_weights(thetas::Vector{Float64},
                            loglikes::Vector{Float64},
                            prior::AbstractPriorSpec)
    @assert length(thetas) == length(loglikes)
    n = length(thetas)
    w = similar(loglikes)
    @inbounds for i in 1:n
        # log posterior ∝ logL + log prior
        w[i] = loglikes[i] + log(max(prior_pdf(prior, thetas[i]), eps()))
    end
    # log-sum-exp normalize
    m = maximum(w)
    @inbounds for i in 1:n
        w[i] = exp(w[i] - m)
    end
    s = sum(w); s > 0 || error("posterior weights sum to zero")
    @. w = w / s
    return w
end

"""
    _posterior_summary(thetas, w) -> Dict{Symbol,Float64}

Return basic discrete summaries on the node grid: posterior mean, median, mode, and a normalization check.
"""
function _posterior_summary(thetas::Vector{Float64}, w::Vector{Float64})
    p = cumsum(w)
    i_med  = searchsortedfirst(p, 0.5)
    i_mode = argmax(w)
    return Dict{Symbol,Float64}(
        :post_mean      => sum(@. thetas * w),
        :post_median    => thetas[clamp(i_med, 1, length(thetas))],
        :post_mode_grid => thetas[i_mode],
        :post_normcheck => p[end],
    )
end

"""
    _choose_representative(s, thetas, metrics) -> θ

Select a representative θ for reporting/subresult extraction (e.g. prior median or grid argmin).
"""
function _choose_representative(s::PriorMarginalizationSettings,
                                thetas::Vector{Float64},
                                metrics::Vector{Float64})
    if s.representative === :prior_median
        return prior_median(s.prior)
    elseif s.representative === :grid_argmin
        return thetas[argmin(metrics)]
    else
        error("unsupported representative $(s.representative)")
    end
end

# ---------------------------
# I/O helpers (dirs, node task, extract metric)
# ---------------------------

"""
    _prepare_node_workdir(base_dir, exec, idx, param, theta) -> String

Create (if needed) and return the node working directory under `base_dir`.
Also writes a small `node_meta.txt` with index/parameter/θ for traceability.
"""
function _prepare_node_workdir(base_dir::AbstractString,
                               exec::PriorExecutionOptions,
                               idx::Int, param::Symbol, theta::Float64)
    dirname = build_node_dirname(exec.node_dir_prefix, idx;
                                 pad=exec.index_pad, param=param, theta=theta,
                                 mode=exec.dir_name_mode, sig=exec.value_sig)
    path = joinpath(base_dir, dirname)
    isdir(path) || mkpath(path)
    write_node_metadata(joinpath(path, "node_meta.txt");
                        index=idx, parameter=String(param), theta=theta)
    return path
end

# Build a per-node task by pinning parameter => theta with the pin-mode's flag.
# function _make_node_task(base_task::SingleTempoTask,
#                          s::PriorMarginalizationSettings,
#                          theta::Float64, workdir::AbstractString)
#     # Expect a project-level helper: tempo_flag(pin_mode)::Int
#     # (If you prefer, replace with your own dispatcher.)
#     flag = tempo_flag(s.pin_mode)
#     return task_with_param(base_task, s.parameter, theta, flag; work_dir=workdir)
# end

# ---------------------------
# Derivers (par_output / temp_dir)
# ---------------------------

"""
    _make_node_task(base_task, s, i, θ) -> (node_task, node_tag, node_dir, par_out)

Construct a per-node task by:
- deriving a leaf `node_dir` name (no leading "nodes/"),
- deriving `node_tag` = basename(node_dir),
- deriving per-node `par_output` name,
- cloning `base_task` with updated `temp_dir`, `par_output`, `io_mirror`,
  and an upserted parameter `(s.parameter => θ)` with pin-mode flag.
"""
function _make_node_task(base_task::SingleTempoTask,
                         s::PriorMarginalizationSettings,
                         i::Int, θ::Float64)

    # 1) build leaf name (no leading "nodes/"), then derive tag and temp_dir
    node_dir = build_node_dirname(s.exec_options.node_dir_prefix, i;
                                   pad  = s.exec_options.index_pad,
                                   param= s.parameter,
                                   theta= θ,
                                   mode = s.exec_options.dir_name_mode,
                                   sig  = s.exec_options.value_sig)

    # tag used in filenames is the last segment only
    node_tag = basename(normpath(String(node_dir)))

    @info "Preparing node $i: node_tag=$(node_tag), θ=$(round(θ, sigdigits=6))"
    @info "  node_dir: $(node_dir)"

    # 3) per-node par_output name (from base par_output stem)
    par_out = task_derive_par_output(base_task, node_tag)

    @info "  par_output: $par_out"

    # 4) build the cloned task
    node_task = task_copy_with(base_task;
        temp_dir   = node_dir,
        par_output = par_out,
        io_mirror  = (:depth_minus, 2),                 # can take from s.exec_options
        override_params_upsert = [TP(String(s.parameter), θ, flag = tempo_flag(s.pin_mode))],
        # optionally: snapshot_par / overwrite / layout / keep_tmp_* etc.
    )

    return node_task, node_tag, node_dir, par_out
end

"""
    _extract_metric(res::GeneralTempoResult, key::Symbol) -> Float64

Fetch a scalar metric from `res.metrics`, with a friendly error if missing.
"""
function _extract_metric(res::GeneralTempoResult, key::Symbol)::Float64
    if haskey(res.metrics, key)
        return res.metrics[key]
    end
    avail = join(sort!(collect(keys(res.metrics))), ", ")
    error("metric :$key not found in result.metrics. Available: [$avail]")
end

# ---------------------------
# Schedulers
# ---------------------------

"""
    _run_nodes_serial!(task, thetas) -> Vector{_PriorNodeResult}

Run node tasks sequentially (optionally in chained mode) and collect node results.
Handles per-node task derivation, execution, metric extraction, and chaining source updates.
"""
function _run_nodes_serial!(task::PriorMarginalizedTempoTask,
                            thetas::Vector{Float64})::Vector{_PriorNodeResult}
    s        = task.settings
    base_dir = task_workdir(task.base_task)
    n        = length(thetas)
    out      = Vector{_PriorNodeResult}(undef, n)

    # traversal order
    idxs = (s.exec_options.chain_direction === :backward) ? collect(n:-1:1) : collect(1:n)

    # relative par_input for the next node (within job_root)
    prev_par_rel::Union{Nothing,String} = nothing

    for i in idxs
        θ = thetas[i]
        node_task, node_tag, node_dir, _ = _make_node_task(task.base_task, s, i, θ)

        # chaining: pass previous node's final par as input
        if s.exec_options.mode === :chained && prev_par_rel !== nothing
            node_task = task_copy_with(node_task;
                par_input    = prev_par_rel,
                snapshot_par = s.exec_options.chain_snapshot_par ? true : nothing,
            )
        end

        res = run_task(node_task)
        m   = _extract_metric(res, s.likelihood_source)

        # working directory (from metadata if available -- should be available)
        run_cwd = get(res.metadata, :run_cwd, joinpath(base_dir, node_dir))

        out[i] = _PriorNodeResult(i, θ, run_cwd, m, prior_pdf(s.prior, θ), res)

        # update chaining source for the next node (make path relative to base_dir)
        if s.exec_options.mode === :chained
            prev_par_abs = get(res.metadata, :par_out_path, nothing)
            if prev_par_abs !== nothing
                prev_par_rel = relpath(prev_par_abs, base_dir)
            else
                prev_par_rel = nothing
            end
        end
    end

    return out
end

function _run_nodes_distributed!(task::PriorMarginalizedTempoTask,
                                 thetas::Vector{Float64})::Vector{_PriorNodeResult}
    error("distributed scheduler: TODO")
end

# ---------------------------
# Assemble top-level result
# ---------------------------

"""
    _integrate_posterior(thetas, metrics, prior; kwargs...) -> NamedTuple

Spline-approximate χ²(θ) from node metrics, combine with the prior to form a continuous
posterior p(θ) ∝ exp(logL(θ))·π(θ), and evaluate evidence, posterior moments, and χ² summaries.

Returns a named tuple with fields:
- `Z` (evidence), `mean`, `var`, `std`, `median`, `mode`, `quad_error`,
- `θgrid`, `pgrid` (normalized density on the returned grid),
- `chi2_post_mean_spline`, `chi2_at_post_mean_spline`, `chi2_marginalized`.
"""
function _integrate_posterior(
    thetas::Vector{Float64},
    metrics::Vector{Float64},
    prior::AbstractPriorSpec;
    likelihood_source::Symbol = :chi2_fit,
    mref::Float64,
    k::Int=3, s::Float64=0.0,          # spline χ²
    refine::Int=801,                   # grid density for quantiles
    rtol::Float64=1e-6
)
    # 1) spline of the metric
    Sχ2 = Spline1D(thetas, metrics; k=k, s=s)
    θmin, θmax = minimum(thetas), maximum(thetas)

    # Prior support; restrict integration to overlap with node range
    θpmin, θpmax = prior_support(prior)
    θlo, θhi = max(θmin, θpmin), min(θmax, θpmax)

    tail_mass = prior_cdf(prior, θlo) + (1.0 - prior_cdf(prior, θhi))
    if tail_mass > 1e-3
        @warn "Non-negligible prior mass outside node interval" tail_mass θlo θhi θmin θmax
    end

    # 2) loglike
    loglike(θ) = likelihood_source === :chi2_fit ? -0.5*(Sχ2(θ) - mref) :
                 error("only :chi2_fit demoed here")

    # 3) logprior
    logprior(θ) = log(max(prior_pdf(prior, θ), eps()))

    # 4) log-posterior
    L(θ) = loglike(θ) + logprior(θ)

    # 5) Stabilization/grid: prefer the prior's native grid (e.g., KDE grid),
    #    otherwise fall back to a uniform grid on [θlo, θhi].
    function _try_prior_grid(pr)
        try
            g = getfield(pr, :_grid)
            if g !== nothing
                gx = getfield(g[], :x)
                return gx
            end
        catch
            return nothing
        end
        return nothing
    end
    θgrid_raw = _try_prior_grid(prior)
    if θgrid_raw isa AbstractVector
        θgrid_v = Float64.(θgrid_raw)
        θgrid_v = filter(x -> x ≥ θlo && x ≤ θhi, θgrid_v)
        issorted(θgrid_v) || sort!(θgrid_v)
        length(θgrid_v) ≥ 3 || (θgrid_v = collect(range(θlo, θhi; length=refine)))
        θgrid = θgrid_v
    else
        θgrid = collect(range(θlo, θhi; length=refine))
    end
    Lgrid = [L(θ) for θ in θgrid]
    c = maximum(Lgrid)

    # 6) integrals
    f(θ) = exp(L(θ) - c)
    Z0, errZ = quadgk(f, θlo, θhi; rtol=rtol)
    Z = Z0 * exp(c)

    m1_0, _ = quadgk(θ -> θ * f(θ), θlo, θhi; rtol=rtol)
    meanθ = (m1_0 * exp(c)) / Z

    # Compute second moment, variance, and std of the posterior
    m2_0, _ = quadgk(θ -> θ^2 * f(θ), θlo, θhi; rtol=rtol)
    varθ = (m2_0 * exp(c)) / Z - meanθ^2
    stdθ = sqrt(max(varθ, 0.0))

    # 7) grid-based CDF/quantiles/mode
    pgrid0 = exp.(Lgrid .- c)
    # normalize with a trapezoidal rule
    function trapz(x, y)
        s = 0.0
        @inbounds for i in 1:length(x)-1
            s += 0.5*(y[i] + y[i+1])*(x[i+1]-x[i])
        end
        s
    end
    Z0_grid = trapz(θgrid, pgrid0)
    pgrid = pgrid0 ./ Z0_grid              # normalized density (∫ pgrid dθ = 1)
    # CDF
    cdf = similar(pgrid); cdf[1] = 0.0
    @inbounds for i in 1:length(θgrid)-1
        cdf[i+1] = cdf[i] + 0.5*(pgrid[i] + pgrid[i+1])*(θgrid[i+1]-θgrid[i])
    end
    # quantiles
    idx_med = searchsortedfirst(cdf, 0.5)
    θ_median = θgrid[clamp(idx_med, 1, length(θgrid))]
    # mode
    i_mode = argmax(pgrid)
    θ_mode = θgrid[i_mode]

    # χ² summaries under the continuous posterior
    # E[χ²(θ)] = ∫ χ²(θ) p(θ) dθ (trapezoidal on the returned grid)
    chi2_post_mean_spline = begin
        acc = 0.0
        @inbounds for i in 1:length(θgrid)-1
            f1 = Sχ2(θgrid[i])   * pgrid[i]
            f2 = Sχ2(θgrid[i+1]) * pgrid[i+1]
            acc += 0.5 * (f1 + f2) * (θgrid[i+1] - θgrid[i])
        end
        acc
    end
    chi2_at_post_mean_spline = Sχ2(meanθ)
    chi2_marginalized = -2 * log(Z) + mref
    return (
        Z = Z,
        mean = meanθ,
        var = varθ,
        std = stdθ,
        median = θ_median,
        mode = θ_mode,
        quad_error = errZ,
        θgrid = θgrid,
        pgrid = pgrid,
        chi2_post_mean_spline = chi2_post_mean_spline,
        chi2_at_post_mean_spline = chi2_at_post_mean_spline,
        chi2_marginalized = chi2_marginalized,
    )
end

"""
    _marginalize_param_estimates(nodes, w) -> Dict{Symbol, ParamEstimate}

Compute posterior-marginalized parameter estimates across nodes using weights `w`.
For each parameter name seen in any node, compute the pooled mean and the total
uncertainty via the law of total variance.
"""
function _marginalize_param_estimates(nodes::Vector{_PriorNodeResult},
                                      w::Vector{Float64})::Dict{Symbol, ParamEstimate}
    # collect all parameter names that ever appeared
    all_names = Set{Symbol}()
    for n in nodes
        r = n.result
        r === nothing && continue
        for k in keys(r.param_estimates)
            push!(all_names, k)
        end
    end

    out = Dict{Symbol, ParamEstimate}()
    for name in all_names
        μs  = Float64[]
        σ2s = Float64[]
        ws  = Float64[]

        for (j, n) in enumerate(nodes)
            r = n.result
            r === nothing && continue
            pe = get(r.param_estimates, name, nothing)
            pe === nothing && continue
            μ = pe.value
            isfinite(μ) || continue
            σ2 = pe.uncertainty
            σ2 = isfinite(σ2) ? σ2^2 : 0.0
            push!(μs, μ); push!(σ2s, σ2); push!(ws, w[j])
        end

        isempty(ws) && continue
        s = sum(ws); s <= 0 && continue
        @. ws = ws / s
        μbar = sum(ws .* μs)
        m2   = sum(ws .* (σ2s .+ μs.^2)) - μbar^2
        σbar = sqrt(max(m2, 0.0))
        out[name] = (value=μbar, uncertainty=σbar)
    end
    return out
end


"""
    _assemble_prior_result(task, thetas, nodes) -> GeneralTempoResult

Assemble the top-level prior-marginalized result:
- choose a reference (θref, mref) and convert node metrics into log-likelihoods;
- compute discrete posterior weights on the node grid and summarize them;
- integrate a continuous (spline-based) posterior for moments/evidence/χ² summaries;
- select a representative node result for the structural fields;
- marginalize model parameter estimates across nodes using posterior weights;
- merge metrics/extras and metadata; optionally attach node subresults.
"""
function _assemble_prior_result(task::PriorMarginalizedTempoTask,
                                thetas::Vector{Float64},
                                nodes::Vector{_PriorNodeResult})
    s = task.settings

    # node metrics and reference point
    node_metrics = getfield.(nodes, :like_metric)
    θref, mref   = _resolve_reference(s, thetas, node_metrics)

    # log-likelihoods relative to the reference
    loglikes = map(node_metrics) do m
        _as_loglike(s.likelihood_source, m, mref)
    end

    # posterior weights and grid summary
    w    = _posterior_weights(thetas, loglikes, s.prior)
    post = _posterior_summary(thetas, w)

    # continuous (spline-based) posterior integration over θ
    integ = _integrate_posterior(
        thetas,
        node_metrics,
        s.prior;
        likelihood_source = s.likelihood_source,
        mref = mref,
    )

    # representative θ and its index
    θrep = _choose_representative(s, thetas, node_metrics)
    irep = findmin(abs.(thetas .- θrep))[2]

    # representative node result
    rep_res = nodes[irep].result
    rep_res === nothing && error("Representative node result is missing")

    # posterior-marginalized model parameter estimates
    param_est_marg = _marginalize_param_estimates(nodes, w)

    # extras: combine discrete (grid) summary and continuous (spline) summary + evidence
    extras = Dict{Symbol,Float64}(
        :ref_theta           => θref,
        :ref_metric          => mref,
        :rep_theta           => θrep,

        # grid (discrete) summaries from node weights
        :post_mean           => post[:post_mean],
        :post_median         => post[:post_median],
        :post_mode_grid      => post[:post_mode_grid],
        :post_normcheck      => post[:post_normcheck],
        :post_ess            => 1.0 / sum(w.^2),

        # spline-based continuous summaries and evidence
        :post_mean_spline    => integ.mean,
        :post_std_spline     => integ.std,
        :post_var_spline     => integ.var,
        :post_median_spline  => integ.median,
        :post_mode_spline    => integ.mode,
        :evidence            => integ.Z,
        :evidence_quad_error => integ.quad_error,
        :chi2_marginalized        => integ.chi2_marginalized,
        :chi2_post_mean_spline    => integ.chi2_post_mean_spline,
        :chi2_at_post_mean_spline => integ.chi2_at_post_mean_spline,
    )
    # if the source is χ², also report grid-average χ²
    if startswith(String(s.likelihood_source), "chi2")
        extras[:chi2_post_mean] = sum(w .* node_metrics)
    end

    # No longer need _GTN__χ2extras_tmp cleanup.

    # prepare subresults if requested
    subres = s.save_node_results ? [n.result for n in nodes] : GeneralTempoResult[]
    subtyp = s.save_node_results ? :prior_nodes : nothing

    # merge metrics safely
    new_metrics = copy(rep_res.metrics)
    merge!(new_metrics, extras)

    # assemble extra metadata
    meta = copy(rep_res.metadata)
    merge!(meta, Dict{Symbol,Any}(
        :prior_parameter   => s.parameter,
        :likelihood_source => s.likelihood_source,
        :ref_strategy      => s.ref_strategy,
        :representative    => s.representative,
        :thetas            => thetas,
        :node_metrics      => node_metrics,
        :node_loglikes     => loglikes,
        :node_weights      => w,
        :node_prior_pdf    => getfield.(nodes, :prior_density),
        :node_workdirs     => getfield.(nodes, :workdir),
        :rep_index         => irep,
        :post_theta_grid   => integ.θgrid,
        :post_pdf_grid     => integ.pgrid,
    ))

    # final aggregate: reuse representative internals, override estimates/metrics, attach extras
    return GeneralTempoResult(
        rep_res.iterations,
        rep_res.final_index,
        rep_res.final,
        rep_res.convergence,
        rep_res.par_file_final,
        param_est_marg,
        rep_res.residuals,
        rep_res.white_noise,
        new_metrics,
        subres,
        subtyp,
        meta,
    )
end

# ---------------------------
# Public API
# ---------------------------
"""
    run_task(task::PriorMarginalizedTempoTask) -> GeneralTempoResult

Execute the prior-marginalized task: build θ grid, run/supply node results,
and assemble the aggregated `GeneralTempoResult`.
"""
function run_task(task::PriorMarginalizedTempoTask)::GeneralTempoResult
    s      = task.settings
    thetas = _build_thetas(s)
    node_results = s.exec_options.scheduler === :serial ?
        _run_nodes_serial!(task, thetas) : _run_nodes_distributed!(task, thetas)
    # node_results = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR/node_results.jld2", "node_results")
    # @info "Loaded $(length(node_results)) node results from JLD2."
    return _assemble_prior_result(task, thetas, node_results)
end

function save_result_jld2(res::GeneralTempoResult; filename::AbstractString="final_result.jld2")
    job_root = res.metadata[:job_root]
    file_path = joinpath(job_root, "results", filename)
    mkpath(dirname(file_path))
    @save file_path res
end

# --------------------------------------------------------------------------------------------------------------