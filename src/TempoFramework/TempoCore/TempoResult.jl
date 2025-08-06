#--------------------------------------------------------------------------------------------------------------

struct BasicResidualStats
    mean::Float64
    wmean::Float64
    median::Float64
    rms::Float64
    wrms::Float64
end

function Base.show(io::IO, stats::BasicResidualStats)
    println(io, "Basic residual stats:")
    println(io, "   mean   = ", @sprintf("%.5f", stats.mean))
    println(io, "   wmean  = ", @sprintf("%.5f", stats.wmean))
    println(io, "   median = ", @sprintf("%.5f", stats.median))
    println(io, "   rms    = ", @sprintf("%.5f", stats.rms))
    println(io, "   wrms   = ", @sprintf("%.5f", stats.wrms))
end

function build_basic_residual_statistics(x::Vector{Float64}, weights::Vector{Float64})::BasicResidualStats
    return BasicResidualStats(
        mean(x),
        sum(x .* weights) / sum(weights),
        median(x),
        sqrt(mean(x .^ 2)),
        sqrt(sum((x .^ 2) .* weights) / sum(weights))
    )
end

#--------------------------------------------------------------------------------------------------------------

struct NormalizedResidualStats
    mean::Float64
    median::Float64
    std::Float64
    skewness::Float64
    kurtosis::Float64
    chisqr::Float64
    ad_statistics::Float64
    ad_p_value::Float64
    ks_statistic::Float64
    ks_p_value::Float64
    jb_statistic::Float64
    jb_p_value::Float64
    min::Float64
    max::Float64
    n_points::Int
end

function Base.show(io::IO, stats::NormalizedResidualStats)
    println(io, "Normalized residual stats:")
    println(io, "   minumum    = ", @sprintf("%.5f", stats.min))
    println(io, "   maximum    = ", @sprintf("%.5f", stats.max))
    println(io, "   mean       = ", @sprintf("%.5f", stats.mean))
    println(io, "   median     = ", @sprintf("%.5f", stats.median))
    println(io, "   std        = ", @sprintf("%.5f", stats.std))
    println(io, "   skewness   = ", @sprintf("%.5f", stats.skewness))
    println(io, "   kurtosis   = ", @sprintf("%.5f", stats.kurtosis))
    println(io, "   chisqr     = ", @sprintf("%.5f", stats.chisqr))
    println(io, "   Anderson-Darling statistic   = ", @sprintf("%.5f", stats.ad_statistics))
    println(io, "   Kolmogorov-Smirnov statistic = ", @sprintf("%.5f", stats.ks_statistics))
    println(io, "   Jarque-Bera statistic        = ", @sprintf("%.5f", stats.jb_statistics))
    println(io, "   number of points = ", stats.n_points)
end

function build_normalized_residual_statistics(x::Vector{Float64})::NormalizedResidualStats
    ad_test = OneSampleADTest(x, Normal(0, 1))
    ks_test = ApproximateOneSampleKSTest(x, Normal(0, 1))
    jb_test = JarqueBeraTest(x)
    return NormalizedResidualStats(
        mean(x),
        median(x),
        std(x),
        skewness(x),
        kurtosis(x),
        sum(x .^ 2),
        ad_test.A²,
        pvalue(ad_test),
        ks_test.δ,
        pvalue(ks_test),
        jb_test.JB,
        pvalue(jb_test),
        minimum(x),
        maximum(x),
        length(x)
    )
end

#--------------------------------------------------------------------------------------------------------------

struct ResidualStatistics
    raw::BasicResidualStats
    tn::BasicResidualStats
    norm_global::NormalizedResidualStats
    norm_local::NormalizedResidualStats
end

function Base.show(io::IO, stats::ResidualStatistics)
    println(io, "Residual Statistics:\n")

    print(io, " → Raw (residual):\n")
    show(IOContext(io, :compact => true), stats.raw)
    println(io)

    print(io, "\n → Whitened (residual_tn):\n")
    show(IOContext(io, :compact => true), stats.tn)
    println(io)

    print(io, "\n → Globally centered normalized whitened ((residual_tn - global_wmean) / uncertainty):\n")
    show(IOContext(io, :compact => true), stats.norm_global)
    println(io)

    print(io, "\n → Locally centered normalized whitened ((residual_tn - local_wmean) / uncertainty):\n")
    show(IOContext(io, :compact => true), stats.norm_local)
end

function build_residual_statistics(
    entries::Vector{DataEntry},
    global_wmean::Float64
    )::ResidualStatistics

    residuals     = [e.residual for e in entries]
    residuals_tn  = [e.residual_tn for e in entries]
    uncertainties = [e.uncertainty for e in entries]
    weights       = [e.weight for e in entries]

    raw = build_basic_residual_statistics(residuals, weights)
    tn  = build_basic_residual_statistics(residuals_tn, weights)

    norm_global = (residuals_tn .- global_wmean) ./ uncertainties
    norm_local  = (residuals_tn .- tn.wmean) ./ uncertainties

    norm_global_stats = build_normalized_residual_statistics(norm_global)
    norm_local_stats  = build_normalized_residual_statistics(norm_local)

    return ResidualStatistics(raw, tn, norm_global_stats, norm_local_stats)
end

#--------------------------------------------------------------------------------------------------------------

struct ResidualStatisticsEntry
    all::ResidualStatistics
    by_backend::Dict{String, ResidualStatistics}
end

function Base.show(io::IO, entry::ResidualStatisticsEntry)
    indent = get(io, :indent, 0)
    println(io, " " ^ indent, "All TOAs:")
    show(IOContext(io, :indent => indent + 4), entry.all)

    println(io, "\n" ^ indent, "Backends:")
    # for (backend, stats) in entry.by_backend
    #     println(io, " " ^ (indent + 2), "[$backend]:")
    #     show(IOContext(io, :indent => indent + 6), stats)
    # end
end

function build_residual_statistics_entry(
    entries::Vector{DataEntry},
    global_wmean::Float64
    )::ResidualStatisticsEntry
    # --- Общая статистика по всем записям
    all_stats = build_residual_statistics(entries, global_wmean)

    # --- Группировка по backend
    grouped = Dict{String, Vector{DataEntry}}()
    for e in entries
        push!(get!(grouped, e.backend, DataEntry[]), e)
    end

    # --- Статистика по каждому backend с тем же global_wmean
    backend_stats = Dict{String, ResidualStatistics}()
    for (backend, group_entries) in grouped
        backend_stats[backend] = build_residual_statistics(group_entries, global_wmean)
    end

    return ResidualStatisticsEntry(
        all_stats,
        backend_stats
    )
end

#--------------------------------------------------------------------------------------------------------------

struct ResidualStatisticsGroup
    in_fit::ResidualStatisticsEntry
    in_tim::ResidualStatisticsEntry
end

function Base.show(io::IO, group::ResidualStatisticsGroup)
    println(io, "Residual Statistics Group:")
    
    println(io, "\n→ In-Fit:")
    show(IOContext(io, :indent => 4), group.in_fit)

    println(io, "\n→ In-TIM:")
    show(IOContext(io, :indent => 4), group.in_tim)
end

function compute_global_wmean(entries::Vector{DataEntry})
    residuals_tn = [e.residual_tn for e in entries]
    weights      = [e.weight for e in entries]
    return sum(residuals_tn .* weights) / sum(weights)
end

function build_residual_statistics_group(entries::Vector{DataEntry})::ResidualStatisticsGroup
    entries_in_fit = filter(e -> e.in_fit, entries)
    global_wmean = compute_global_wmean(entries_in_fit)

    return ResidualStatisticsGroup(
        build_residual_statistics_entry(entries_in_fit, global_wmean),
        build_residual_statistics_entry(entries, global_wmean),
    )
end

#--------------------------------------------------------------------------------------------------------------

struct WhiteNoiseBackendFitResult
    backend::Union{String, Symbol}
    EFAC::Float64
    EQUAD::Float64
    offset::Float64
    AD_objective::Float64
    stats::NormalizedResidualStats
end

struct WhiteNoiseFitResult
    by_backend::Dict{Symbol, WhiteNoiseBackendFitResult}
    global_stats::NormalizedResidualStats
end

function groupby_backend(entries::Vector{CombinedTOAEntry}; in_fit::Bool = true)
    groups = Dict{Symbol, Vector{CombinedTOAEntry}}()
    for entry in entries
        if entry.in_fit != in_fit
            continue
        end
        backend = Symbol(entry.backend)
        if !haskey(groups, backend)
            groups[backend] = CombinedTOAEntry[]
        end
        push!(groups[backend], entry)
    end
    return groups
end

function build_white_noise_fit(
    combined_entries::Vector{CombinedTOAEntry}
    )::WhiteNoiseFitResult
    by_backend = Dict{Symbol, WhiteNoiseBackendFitResult}()

    # Сбор нормализованных резидуалов по всем бэкендам
    global_normalized_residuals = Float64[]
    global_uncertainties_orig = Float64[]

    grouped_entries = groupby_backend(combined_entries, in_fit = true)

    for (backend, entries) in grouped_entries
        residuals = [e.residual_tn for e in entries]
        uncertainties_orig = [e.uncertainty_orig for e in entries]

        # Выполняем оценку параметров белого шума
        EFAC, EQUAD, offset, AD_obj = estimate_WhiteNoise_AD_with_offset(
            residuals,
            uncertainties_orig;
            backend = string(backend),
            print_results = false,
            plot_results = false,
        )

        # Строим трансформированные ошибки и нормализуем
        uncertainties_transformed = similar(uncertainties_orig)
        transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)

        residuals_norm = (residuals .- offset) ./ uncertainties_transformed

        # Статистика по нормализованным остаткам
        stats = build_normalized_residual_statistics(residuals_norm)

        result = WhiteNoiseBackendFitResult(
            backend,
            EFAC,
            EQUAD,
            offset,
            AD_obj,
            length(residuals),
            stats
        )

        by_backend[backend] = result

        # Для глобальной статистики
        append!(global_normalized_residuals, residuals_norm)
        append!(global_uncertainties_orig, uncertainties_orig)  # может пригодиться
    end

    # Глобальная статистика
    global_stats = build_normalized_residual_statistics(global_normalized_residuals)

    return WhiteNoiseFitResult(by_backend, global_stats)
end


#--------------------------------------------------------------------------------------------------------------

struct InternalIterationResult
    output::InternalIterationOutput
    residuals::Union{Vector{RawResidualEntry}, Nothing}
    stats::Union{ResidualStatisticsGroup, Nothing}
    white_noise_fit::Union{WhiteNoiseFitResult, Nothing}
    metadata::Dict{Symbol, Any}
end

function Base.show(io::IO, result::InternalIterationResult)
    println(io, "Internal Iteration Result:")
    
    # Output summary
    println(io, "  ▸ Basic TEMPO output:")
    show(io, result.output.basic)

    if iserror(result.output.error)
        println(io, "  ▸ TEMPO error:")
        show(io, result.output.error)
    else
        println(io, "  ▸ TEMPO error: none")
    end

    # Residual stats
    if isnothing(result.stats)
        println(io, "  ▸ Residual statistics: none")
    else
        println(io, "  ▸ Residual statistics (in-fit only):")
        show(io, result.stats.in_fit)
    end

    # Residuals presence
    if isnothing(result.residuals)
        println(io, "  ▸ Residuals: not stored")
    else
        println(io, "  ▸ Residuals: $(length(result.residuals)) entries loaded")
    end

    # Metadata summary
    if isempty(result.metadata)
        println(io, "  ▸ Metadata: none")
    else
        println(io, "  ▸ Metadata: $(join(keys(result.metadata), ", "))")
    end
end

function build_internal_iteration_result(
    output::InternalIterationOutput,
    residual_path::Union{Nothing, String},
    tim_entries::Vector{TimTOAEntry};
    save_residuals::Bool = false,
    analyze_white_noise::Bool = false,
    time_start::Union{Nothing, Float64} = nothing,
    time_finish::Union{Nothing, Float64} = nothing
)::InternalIterationResult

    # Обработка случая с ошибкой
    if output.error ≠ TempoOutputError()
        return InternalIterationResult(
            output,
            nothing,
            nothing,
            nothing,  # white_noise_fit
            Dict(:error => output.error),
        )
    end

    # Чтение резидуалов, если файл существует
    residuals = isnothing(residual_path) ? nothing : read_residual_file(residual_path)

    # Комбинация TOA и резидуалов
    combined_entries = !isnothing(residuals) ? 
        combine_tim_and_residuals(
            tim_entries, residuals;
            time_start = time_start,
            time_finish = time_finish
        ) :
        nothing

    # Статистика по резидуалам
    stats = isnothing(combined_entries) ? nothing :
            build_residual_statistics_group(combined_entries)

    # Анализ белого шума (по объединённым данным)
    white_noise_fit = (analyze_white_noise && !isnothing(combined_entries)) ?
            build_white_noise_fit(combined_entries) :
            nothing

    # Возврат результата
    return InternalIterationResult(
        output,
        save_residuals ? residuals : nothing,
        stats,
        white_noise_fit,
        Dict{Symbol, Any}()  # для логов или отладочной информации
    )
end

#--------------------------------------------------------------------------------------------------------------

struct ConvergenceSeries
    values::Vector{Float64}                      # Значения метрики по итерациям
    abs_deltas::Vector{Float64}                  # abs(Δₙ) = abs(xₙ - xₙ₋₁)
    rel_deltas::Vector{Float64}                  # rel(Δₙ) = abs(Δₙ / xₙ₋₁)
    final_abs_delta::Union{Float64, Nothing}     # Последнее абсолютное изменение
    final_rel_delta::Union{Float64, Nothing}     # Последнее относительное изменение
end

function build_convergence_series(values::Vector{Float64})::ConvergenceSeries
    N = length(values)

    if N < 2
        return ConvergenceSeries(values, Float64[], Float64[], nothing, nothing)
    end

    abs_deltas = [abs(values[i+1] - values[i]) for i in 1:N-1]
    rel_deltas = [abs_deltas[i] / (abs(values[i]) + eps()) for i in 1:N-1]  # стабильно при делении

    return ConvergenceSeries(
        values,
        abs_deltas,
        rel_deltas,
        abs_deltas[end],
        rel_deltas[end]
    )
end

#--------------------------------------------------------------------------------------------------------------

function find_worst_parameter(fit_params::Vector{FitParameter})
    worst_ratio = 0.0
    worst_parameter = nothing

    for p in fit_params
        if p.fit_flag && isfinite(p.uncertainty) && p.uncertainty ≠ 0
            ratio = p.difference / p.uncertainty
            if abs(ratio) > abs(worst_ratio)
                worst_ratio = ratio
                worst_parameter = (
                    name = p.name_symbol,
                    delta = p.difference,
                    uncertainty = p.uncertainty,
                    ratio = ratio,
                )
            end
        end
    end

    return worst_parameter
end

#--------------------------------------------------------------------------------------------------------------

struct ConvergenceInfo
    wrms::ConvergenceSeries
    wrms_tn::ConvergenceSeries
    chisqr::ConvergenceSeries
    worst_parameter::Union{NamedTuple, Nothing}  # ← новое поле

    final_pre_post::Float64
    threshold::Dict{Symbol, Float64}
    converged::Bool
end

function Base.show(io::IO, info::ConvergenceInfo)
    println(io, "ConvergenceInfo:")
    println(io, "  ✓ Converged: ", info.converged ? "✅ Yes" : "❌ No")
    println(io)

    metrics = [
        (:wrms,        info.wrms),
        (:wrms_tn,     info.wrms_tn),
        (:chisqr,      info.chisqr),
    ]

    for (name, conv) in metrics
        println(io, "  ── $name ───────────────────────────")
        println(io, "     last value        : ", @sprintf("%.6f", conv.values[end]))
        println(io, "     abs delta         : ", isnothing(conv.final_abs_delta) ? "-" : @sprintf("%.6e", conv.final_abs_delta))
        println(io, "     rel delta         : ", isnothing(conv.final_rel_delta) ? "-" : @sprintf("%.6e", conv.final_rel_delta))

        thresh_abs = get(info.threshold, Symbol("abs_$name"), nothing)
        thresh_rel = get(info.threshold, Symbol("rel_$name"), nothing)

        if thresh_abs !== nothing
            println(io, "     threshold abs     : ", @sprintf("%.1e", thresh_abs))
        end
        if thresh_rel !== nothing
            println(io, "     threshold rel     : ", @sprintf("%.1e", thresh_rel))
        end
        println(io)
    end

    # final_pre_post
    println(io, "  ── final_pre_post ───────────────")
    println(io, "     value             : ", @sprintf("%.8f", info.final_pre_post))
    if haskey(info.threshold, :final_pre_post)
        println(io, "     threshold         : ", @sprintf("%.1e", info.threshold[:final_pre_post]))
    end
    println(io)

    # worst_param (если есть)
    if info.worst_parameter !== nothing
        wp = info.worst_parameter
        println(io, "  ── worst_parameter ─────────────────────")
        println(io, "     name              : ", wp.name)
        println(io, "     Δ (difference)    : ", @sprintf("%.6g", wp.delta))
        println(io, "     σ (uncertainty)   : ", @sprintf("%.6g", wp.uncertainty))
        println(io, "     Δ / σ             : ", @sprintf("%.3g", wp.ratio))
    end
end

function get_convergence_metric(info::ConvergenceInfo, key::Symbol)::Float64
    if key == :abs_wrms
        return info.wrms.final_abs_delta
    elseif key == :rel_wrms
        return info.wrms.final_rel_delta
    elseif key == :abs_wrms_tn
        return info.wrms_tn.final_abs_delta
    elseif key == :rel_wrms_tn
        return info.wrms_tn.final_rel_delta
    elseif key == :abs_chisqr
        return info.chisqr.final_abs_delta
    elseif key == :rel_chisqr
        return info.chisqr.final_rel_delta
    elseif key == :final_pre_post
        return abs(info.final_pre_post - 1.0)
    else
        throw(ArgumentError("Unknown convergence key: $key"))
    end
end

function default_convergence_thresholds()::Dict{Symbol, Float64}
    Dict(
        :abs_wrms_tn     => 1e-2,
        :abs_chisqr      => 1e-2,
        :final_pre_post  => 1e-6
    )
end

function is_converged(info::ConvergenceInfo)::Bool
    for (key, limit) in info.threshold
        actual = get_convergence_metric(info, key)
        if isnan(actual) || actual > limit
            return false
        end
    end
    return true
end

function is_converged_by(info::ConvergenceInfo, keys::Symbol...)
    all(key -> get_convergence_metric(info, key) ≤ info.threshold[key], keys)
end


function build_convergence_info(
    iterations::Vector{InternalIterationResult};
    threshold::Dict{Symbol, Float64} = default_convergence_thresholds()
    )::ConvergenceInfo

    wrms_vals     = Float64[]
    wrms_tn_vals  = Float64[]
    chisqr_vals   = Float64[]

    for iter in iterations
        if isnothing(iter.stats)
            continue  # пропускаем итерации с ошибками
        end
        push!(wrms_vals, iter.stats.in_fit.all.raw.wrms)
        push!(wrms_tn_vals, iter.stats.in_fit.all.tn.wrms)
        push!(chisqr_vals, iter.stats.in_fit.all.norm_global.chisqr)
    end

    wrms     = build_convergence_series(wrms_vals)
    wrms_tn  = build_convergence_series(wrms_tn_vals)
    chisqr   = build_convergence_series(chisqr_vals)

    # --- Итоговая итерация (последняя успешная) ---
    final_success = findlast(x -> x.output.error == TempoOutputError(), iterations)

    if isnothing(final_success)
        # Всё сломалось — возвращаем пустой объект с converged = false
        return ConvergenceInfo(
            wrms,
            wrms_tn,
            chisqr,
            nothing,        # worst_parameter
            Inf,            # final_pre_post
            threshold,
            false
        )
    end

    last_ok = iterations[final_success]
    final_pre_post = last_ok.output.basic.pre_post
    worst_parameter = find_worst_parameter(last_ok.output.fit_parameters)

    info = ConvergenceInfo(
        wrms,
        wrms_tn,
        chisqr,
        worst_parameter,
        final_pre_post,
        threshold,
        false
    )

    return ConvergenceInfo(
        wrms,
        wrms_tn,
        chisqr,
        worst_parameter,
        final_pre_post,
        threshold,
        is_converged(info)
    )
end

#--------------------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------------------

struct GeneralTempoResult
    all_internal_iterations::Union{Vector{InternalIterationResult}, Nothing}
    final_result::InternalIterationResult

    convergence_info::ConvergenceInfo
    final_parfile::Union{TempoParFile, Nothing}

    model_parameters::Dict{Symbol, @NamedTuple{value::Float64, uncertainty::Float64}}
    residual_stats::Dict{Symbol, Any}
    white_noise_fit::Union{Dict{String, Tuple{Float64, Float64}}, Nothing}
    
    subresults::Union{Vector{GeneralTempoResult}, Nothing}
    subresult_type::Union{Symbol, Nothing}

    metadata::Dict{Symbol, Any}
end
