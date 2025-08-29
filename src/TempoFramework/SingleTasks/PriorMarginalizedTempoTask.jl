# src/TempoFramework/SingleTasks/PriorMarginalizedTempoTask.jl
# Prior-marginalized wrapper around a SingleTempoTask:
# runs per-node tasks over a theta grid and aggregates results under a given prior.

# --------------------------------------------------------------------------------------------------------------
# Task wrapper
# --------------------------------------------------------------------------------------------------------------

"""
    PriorMarginalizedTempoTask{T,P,N} <: SingleTempoTask

Wrapper that runs a base single-tempo task on a grid of `theta` values and aggregates
results according to the prior `P` and node rule `N`.
- `base_task`  : the underlying single task to run at each grid node
- `settings`   : prior definition, node grid rule, execution options, etc.
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

Internal container for a single grid node.

Fields
- `index`         — 1-based node position on the grid
- `theta`         — theta value at the node
- `workdir`       — absolute path to the node's working directory
- `like_metric`   — scalar metric used as likelihood proxy (e.g. `metrics[:chi2_fit]`)
- `prior_density` — prior pdf value π(theta) at this node
- `result`        — full `GeneralTempoResult` produced by the node
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
_build_thetas(s::PriorMarginalizationSettings) = eval_nodes(s.prior, s.nodes)

# Return only mref used to shift chi2-like sources.
_resolve_mref(::PriorMarginalizationSettings, thetas::Vector{Float64}, metrics::Vector{Float64}) =
    minimum(metrics)  # nodes-only baseline; no spline/interp

# Return (θref, mref) where mref is the metric at the chosen reference
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

# Use mref only for chi2-like; ignore otherwise.
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
function _make_node_task(base_task::SingleTempoTask,
                         s::PriorMarginalizationSettings,
                         theta::Float64, workdir::AbstractString)
    # Expect a project-level helper: tempo_flag(pin_mode)::Int
    # (If you prefer, replace with your own dispatcher.)
    flag = tempo_flag(s.pin_mode)
    return task_with_param(base_task, s.parameter, theta, flag; work_dir=workdir)
end

# Delegate metric extraction to a small adapter you implement for your GeneralTempoResult
# e.g. result_metric(res::GeneralTempoResult, ::Val{:chi2_fit}) = ...
function _extract_metric(res::GeneralTempoResult, key::Symbol)::Float64
    if haskey(res.metrics, key)
        return res.metrics[key]
    end
    # дружелюбная ошибка со списком доступных ключей
    avail = join(sort!(collect(keys(res.metrics))), ", ")
    error("metric :$key not found in result.metrics. Available: [$avail]")
end

# ---------------------------
# Schedulers
# ---------------------------

function _run_nodes_serial!(task::PriorMarginalizedTempoTask,
                            thetas::Vector{Float64})::Vector{_PriorNodeResult}
    s = task.settings
    base_dir = task_workdir(task.base_task)
    out = Vector{_PriorNodeResult}(undef, length(thetas))

    for (i, θ) in enumerate(thetas)
        wdir = _prepare_node_workdir(base_dir, s.exec_options, i, s.parameter, θ)

        # build node-task and stage inputs (name-only scheme)
        node_task = _make_node_task(task.base_task, s, θ, wdir)
        task_stage_inputs!(node_task, wdir)

        res = run_task(node_task)
        m   = _extract_metric(res, s.likelihood_source)

        out[i] = _PriorNodeResult(
            i, θ, wdir, m,
            prior_pdf(s.prior, θ),
            res,         
        )
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

function _integrate_posterior(
    thetas::Vector{Float64},
    metrics::Vector{Float64},
    prior::AbstractPriorSpec;
    likelihood_source::Symbol = :chi2_fit,
    mref::Float64,
    k::Int=3, s::Float64=0.0,          # сплайн χ²
    refine::Int=801,                   # плотность сетки для квантилей
    rtol::Float64=1e-6
)
    # 1) сплайн метрики
    Sχ2 = Spline1D(thetas, metrics; k=k, s=s)
    θmin, θmax = minimum(thetas), maximum(thetas)

    # 2) loglike
    loglike(θ) = likelihood_source === :chi2_fit ? -0.5*(Sχ2(θ) - mref) :
                 error("only :chi2_fit demoed here")

    # 3) logprior
    logprior(θ) = log(max(prior_pdf(prior, θ), eps()))

    # 4) лог-постериор
    ℓ(θ) = loglike(θ) + logprior(θ)

    # 5) стабилизация: возьмём максимум на плотной сетке
    θgrid = range(θmin, θmax; length=refine)
    ℓgrid = [ℓ(θ) for θ in θgrid]
    c = maximum(ℓgrid)

    # 6) интегралы
    f(θ) = exp(ℓ(θ) - c)
    Z0, errZ = quadgk(f, θmin, θmax; rtol=rtol)
    Z = Z0 * exp(c)

    m1_0, _ = quadgk(θ -> θ * f(θ), θmin, θmax; rtol=rtol)
    meanθ = (m1_0 * exp(c)) / Z

    # 7) сетка для CDF/квантилей/моды
    pgrid0 = exp.(ℓgrid .- c)
    # нормировка трапецией
    function trapz(x, y)
        s = 0.0
        @inbounds for i in 1:length(x)-1
            s += 0.5*(y[i] + y[i+1])*(x[i+1]-x[i])
        end
        s
    end
    Z0_grid = trapz(θgrid, pgrid0)
    pgrid = (pgrid0 ./ Z0_grid) .* exp(c)   # нормированная плотность
    # CDF
    cdf = similar(pgrid); cdf[1] = 0.0
    @inbounds for i in 1:length(θgrid)-1
        cdf[i+1] = cdf[i] + 0.5*(pgrid[i] + pgrid[i+1])*(θgrid[i+1]-θgrid[i])
    end
    # квантили
    idx_med = searchsortedfirst(cdf, 0.5)
    θ_median = θgrid[clamp(idx_med, 1, length(θgrid))]
    # мода
    i_mode = argmax(pgrid)
    θ_mode = θgrid[i_mode]

    return (
        Z = Z,
        mean = meanθ,
        median = θ_median,
        mode = θ_mode,
        quad_error = errZ,
        θgrid = θgrid,
        pgrid = pgrid,
    )
end

# Маргинализация оценок параметров по весам постериора
function _marginalize_param_estimates(nodes::Vector{_PriorNodeResult},
                                      w::Vector{Float64})::Dict{Symbol, ParamEstimate}
    # набор всех имён, которые вообще появлялись
    all_names = Set{Symbol}()
    for n in nodes
        r = n.base_result
        r === nothing && continue
        for k in keys(r.param_estimates)
            push!(all_names, k)
        end
    end

    out = Dict{Symbol, ParamEstimate}()
    for name in all_names
        μs  = Float64[]   # значения
        σ2s = Float64[]   # дисперсии
        ws  = Float64[]   # веса тех узлов, где параметр есть

        for (j, n) in enumerate(nodes)
            r = n.base_result
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


function _assemble_prior_result(task::PriorMarginalizedTempoTask,
                                thetas::Vector{Float64},
                                nodes::Vector{_PriorNodeResult})
    s = task.settings

    # узловые метрики и опорная точка
    node_metrics = getfield.(nodes, :metric)
    θref, mref   = _resolve_reference(s, thetas, node_metrics)

    # лог-правдоподобия относительно опоры
    loglikes = map(node_metrics) do m
        _as_loglike(s.likelihood_source, m, mref)
    end

    # постериорные веса и краткая сводка
    w    = _posterior_weights(thetas, loglikes, s.prior)
    post = _posterior_summary(thetas, w)

    # представитель и его индекс
    θrep = _choose_representative(s, thetas, node_metrics)
    irep = findmin(abs.(thetas .- θrep))[2]

    # результат представителя (должен быть, мы его всегда сохраняем в раннере)
    rep_res = nodes[irep].base_result
    rep_res === nothing && error("Representative node result is missing")

    # маргинализованные оценки параметров модели
    param_est_marg = _marginalize_param_estimates(nodes, w)

    # дополнительные метрики поверх core-метрик representative-узла
    extras = Dict{Symbol,Float64}(
        :ref_theta      => θref,
        :ref_metric     => mref,
        :rep_theta      => θrep,
        :post_mean      => post[:post_mean],
        :post_median    => post[:post_median],
        :post_mode_grid => post[:post_mode_grid],
        :post_normcheck => post[:post_normcheck],
        :post_ess       => 1.0 / sum(w.^2),
    )
    # если источник — χ², добавим среднее χ² по постериору
    if startswith(String(s.likelihood_source), "chi2")
        extras[:chi2_post_mean] = sum(w .* node_metrics)
    end

    # подготовим subresults по флагу
    subres = s.exec_options.save_node_results ? [n.base_result for n in nodes] : GeneralTempoResult[]
    subtyp = s.exec_options.save_node_results ? :prior_nodes : nothing

    # аккуратно смёржим метрики
    new_metrics = copy(rep_res.metrics)
    merge!(new_metrics, extras)

    # соберём полезную мета-информацию
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
        :node_prior_pdf    => getfield.(nodes, :prior_pdf),
        :node_workdirs     => getfield.(nodes, :workdir),
        :rep_index         => irep,
    ))

    # финальный агрегат: берём «внутренности» representative-результата,
    # заменяем param_estimates и metrics, добавляем (опционально) subresults и метаданные
    return GeneralTempoResult(
        rep_res.iterations,
        rep_res.final_index,
        rep_res.final,
        rep_res.convergence,
        rep_res.par_file_final,
        param_est_marg,                 # ← маргинализованные оценки параметров
        rep_res.residuals,
        rep_res.white_noise,
        new_metrics,                    # ← метрики + extras
        subres,
        subtyp,
        meta,
    )
end

# ---------------------------
# Public API
# ---------------------------
function run_task(task::PriorMarginalizedTempoTask)::GeneralTempoResult
    s      = task.settings
    thetas = _build_thetas(s)
    node_results = s.exec_options.scheduler === :serial ?
        _run_nodes_serial!(task, thetas) : _run_nodes_distributed!(task, thetas)
    return _assemble_prior_result(task, thetas, node_results)
end