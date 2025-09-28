# Results and Diagnostics

How to read, compare, and diagnose results produced by Tempo tasks.

- Top-level: `GeneralTempoResult` — the container you get from `run_task`.
- Iteration-level: `InternalIterationResult` — per TEMPO iteration outputs and stats.
- Residual statistics: raw/whitened/normalized, grouped by backend and by in-fit vs in-tim.
- White-noise fit: EFAC/EQUAD/offset per backend with normality diagnostics.
- Compact metrics: `res.metrics` for quick scoring and ranking.

This page uses the same types your tests and examples return, so you can copy/paste the snippets below into a REPL.

## GeneralTempoResult

The top-level container with everything you usually want after a run.

Key properties:
- `res.iterations::Vector{InternalIterationResult}` — all iterations, in order
- `res.final` — alias of the last iteration (`res.iterations[end]`)
- `res.last_successful` — last iteration with no engine/parse error and computed stats (or `nothing`)
- `res.success::Bool` and `res.status::Symbol` — quick verdict and reason
- `res.convergence` — summary across iterations (wrms_tn and chi2 deltas)
- `res.metrics::Dict{Symbol,Float64}` — small set of scalar numbers for ranking/plots
- `res.param_estimates::Dict{Symbol, (value, uncertainty)}` — final fit parameters
- `res.par_file_final::Union{TempoParFile,Nothing}` — the output par-file if available
- `res.subresults` — nested results (e.g., nodes, grid cells) when a task runs many jobs
- `res.metadata::Dict{Symbol,Any}` — timings, paths, seeds, policies, etc.

Convenience fields for quick access:
- `res.residual_stats` — alias to `res.final.stats`
- `res.white_noise_fit` — alias to `res.final.white_noise_fit`

Example checks:
```julia
res.success          # true/false
res.status           # :ok | :engine_failed | :parse_failed | :files_missing | :unknown
res.final.output     # parsed TEMPO basic block and fit table
res.final.stats      # residual statistics (see below)
res.metrics[:wrms_fit]  # weighted RMS (fit window, whitened)
```

## Compact metrics (res.metrics)

We compute a compact, task-agnostic set of scalars from the final iteration and convergence:
- `:chi2_fit_basic` — Fit chi-square reported by TEMPO (basic block)
- `:wrms_fit`, `:wrms_tn_fit`, `:chi2_fit`, `:chi2r_fit` — from residual stats (in-fit set)
- `:pre_post_final` — TEMPO pre/post RMS ratio for the final iteration
- `:delta_wrms_tn`, `:delta_chi2` — absolute deltas between last two iterations (convergence)
- `:ad_white_fit` — global Anderson–Darling A² after a white-noise fit, if available

When stats include a separate in-TIM set and you enable it, these may also appear:
- `:wrms_tim`, `:wrms_tn_tim`, `:chi2_tim`, `:chi2r_tim`

Safe lookup helper:
```julia
using GravityToolsNext: result_metric
wrms = result_metric(res, :wrms_fit)        # NaN if missing
chi2 = result_metric(res, :chi2_fit)
```

Tip: Use `isfinite` when ranking many results and fall back to NaN-aware sorting.

## Residual statistics

`InternalIterationResult.stats` holds a `ResidualStatisticsGroup` with two entries:
- `in_fit` — statistics restricted to the fit window
- `in_tim` — statistics over all TOAs (may equal `in_fit` if there is no time window)

Each entry has:
- `all::ResidualStatistics` — overall (not split by backend)
- `by_backend::Dict{Symbol,ResidualStatistics}` — per-backend stats

And each `ResidualStatistics` includes four views:
- `raw` — basic stats of residuals
- `tn` — basic stats of whitened residuals (after TN plugin)
- `norm_global` — normalized whitened residuals centered by a global weighted mean
- `norm_local` — normalized whitened residuals centered by a local weighted mean

Quick peek:
```julia
stats = res.residual_stats  # == res.final.stats
show(stats)                 # pretty, multi-section text/plain view

# Pull a few numbers
wrms_fit = stats.in_fit.all.tn.wrms
rchi2    = stats.in_fit.all.norm_global.red_chisqr
n_by_be  = Dict(k => v.norm_global.n for (k,v) in stats.in_fit.by_backend)
```

Helpers:
- `in_fit_equals_in_tim(stats_group)` — tells if the same entry is reused (no time window)

## White-noise fit (EFAC, EQUAD, offset)

When requested via settings, the final iteration may include a per-backend white-noise calibration:
```julia
wn = res.white_noise_fit   # ::Union{WhiteNoiseFitResult,Nothing}
wn === nothing && @info "No white-noise fit performed"
```

If present:
- `wn.by_backend::Dict{Symbol,WhiteNoiseBackendFitResult}` — parameters and diagnostics per backend
- `wn.failed_backends::Vector{Symbol}` — those that failed the solver or produced non-finite params
- `wn.global_stats::NormalizedResidualStats` — summary over concatenated normalized residuals from successful backends

Each backend entry stores:
- `efac::Float64`, `equad::Float64` (µs), and `offset::Float64` (µs)
- `ad_objective::Float64` — minimized A²
- `stats::NormalizedResidualStats` — detailed normality metrics for the calibrated residuals
- `success::Bool`, `converged::Bool`

Printing a rich, aligned report:
```julia
using GravityToolsNext: print_white_noise_fit_report
for (be, bres) in wn.by_backend
    # Reconstruct the two vectors the report expects (residual_tn and original uncertainties)
    # If you kept CombinedTOAEntry vectors per backend, pass them directly; otherwise skip.
    # print_white_noise_fit_report(be, residuals_tn, uncertainties_orig, bres.efac, bres.equad, bres.offset, bres.ad_objective)
end

# Or quickly inspect the global stats
show(wn.global_stats)
```

## Parameter estimates

The final fit table is summarized into `res.param_estimates::Dict{Symbol,NamedTuple}` with
(value, uncertainty) for each parameter that has finite values:
```julia
est = res.param_estimates
get(est, :F0, (value=NaN, uncertainty=NaN))
```

## Success and convergence

- `res.success` — by default, mirrors final iteration error: `!iserror(res.final.output.error)`
- `res.convergence` — contains fields with deltas and a `converged::Bool` flag summarizing the path

You can build your own criteria by combining `res.metrics` and `res.convergence`.

## Working with nested results

Tasks that execute multiple jobs (prior marginalization, adaptive grid) return a parent `GeneralTempoResult`
with `res.subresults::Vector{GeneralTempoResult}` and a tag `res.subresult_type` (e.g., `:node` or `:grid`).
Traverse or aggregate using the same fields and helpers described above:
```julia
parent = res
for child in parent.subresults
    @info "child wrms" result_metric(child, :wrms_fit)
end
```

## Troubleshooting

- Missing metrics: use `result_metric(res, key)` which returns `NaN` if absent.
- No time window: expect `res.residual_stats.in_fit === res.residual_stats.in_tim`.
- White-noise failures: check `wn.failed_backends` and per-backend `success=false` entries.
- Engine crashes: inspect `res.final.output.error` and `res.metadata[:status]` if provided.
