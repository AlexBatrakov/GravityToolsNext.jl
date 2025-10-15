# TempoFramework/Prior/PriorSettings.jl

# --------------------------------------------------------------------------------------------------------------
# Execution options for running prior-marginalization nodes
# --------------------------------------------------------------------------------------------------------------

_choices_str(allowed) = join((":" .* string.(allowed)), " | ")

@inline function validate_choice(name::AbstractString, val::Symbol, allowed)::Symbol
    val in allowed || throw(ArgumentError("$name must be $(_choices_str(allowed)); got :$val"))
    return val
end

const _MODES            = (:independent, :chained)
const _SCHEDULERS       = (:serial, :distributed)
const _ON_ERROR         = (:stop, :skip, :collect)
const _WORKDIR_LAYOUT   = (:per_node,)
const _DIRNAME_MODE     = (:index_only, :with_value)
const _CHAIN_DIRECTIONS = (:forward, :backward)

"""
    PriorExecutionOptions

Execution policy for running prior-marginalization nodes.

Fields
- `mode`              : node chaining policy (`:independent` | `:chained`)
- `chain_direction`   : traversal order when chaining (`:forward` | `:backward`)
- `chain_snapshot_par`: force snapshot of per-node input par when chaining (true/false)
- `scheduler`         : how to schedule nodes (`:serial` | `:distributed`)
- `max_workers`       : cap for distributed workers (ignored for `:serial`; `0` means "all available")
- `workdir_layout`    : directory layout policy (currently only `:per_node`)
- `node_dir_prefix`   : subdirectory prefix under the base task's work dir (e.g. `"nodes/node_"`)
- `keep_node_dirs`    : keep (true) or delete (false) per-node directories after run
- `on_error`          : error handling (`:stop` | `:skip` | `:collect`)

Directory naming:
- `dir_name_mode`     : directory name style (`:index_only` | `:with_value`)
- `index_pad`         : zero-padding width for node index
- `value_sig`         : significant digits used when embedding the parameter value into the name
"""
struct PriorExecutionOptions
    mode::Symbol             # :independent | :chained
    chain_direction::Symbol  # :forward | :backward
    chain_snapshot_par::Bool # force snapshot per-node input when chaining
    scheduler::Symbol        # :serial | :distributed
    max_workers::Int
    workdir_layout::Symbol   # currently only :per_node
    node_dir_prefix::String  # e.g. "nodes/node_"
    keep_node_dirs::Bool
    on_error::Symbol         # :stop | :skip | :collect

    dir_name_mode::Symbol    # :index_only | :with_value
    index_pad::Int           # zero-padding width for node index
    value_sig::Int           # significant digits for value in name

    function PriorExecutionOptions(;
        mode::Symbol             = :independent,
        chain_direction::Symbol  = :forward,
        chain_snapshot_par::Bool = false,
        scheduler::Symbol        = :serial,
        max_workers::Int         = 0,
        workdir_layout::Symbol   = :per_node,
        node_dir_prefix::String  = "nodes/node_",
        keep_node_dirs::Bool     = true,
        on_error::Symbol         = :stop,
        dir_name_mode::Symbol    = :index_only,   # or :with_value
        index_pad::Int           = 3,
        value_sig::Int           = 6,
    )
        validate_choice("mode",             mode,             _MODES)
        validate_choice("chain_direction",  chain_direction,  _CHAIN_DIRECTIONS)
        validate_choice("scheduler",        scheduler,        _SCHEDULERS)
        validate_choice("workdir_layout",   workdir_layout,   _WORKDIR_LAYOUT)
        validate_choice("on_error",         on_error,         _ON_ERROR)
        validate_choice("dir_name_mode",    dir_name_mode,    _DIRNAME_MODE)
        index_pad >= 0 || error("index_pad must be ≥ 0")
        value_sig >= 1 || error("value_sig must be ≥ 1")
        !isempty(node_dir_prefix) || error("node_dir_prefix must not be empty")
        new(mode, chain_direction, chain_snapshot_par, scheduler, max_workers, workdir_layout, node_dir_prefix,
            keep_node_dirs, on_error, dir_name_mode, index_pad, value_sig)
    end
end

function Base.show(io::IO, ::MIME"text/plain", e::PriorExecutionOptions)
    indent = get(io, :indent, 0)
    pad, spad = repeat(" ", indent), repeat(" ", indent + 2)
    println(io, pad, "PriorExecutionOptions")
    println(io, spad, "mode:              ", e.mode)
    println(io, spad, "chain_direction:   ", e.chain_direction)
    println(io, spad, "chain_snapshot_par:", e.chain_snapshot_par)
    println(io, spad, "scheduler:         ", e.scheduler)
    println(io, spad, "max_workers:       ", e.max_workers)
    println(io, spad, "workdir_layout:    ", e.workdir_layout)
    println(io, spad, "node_dir_prefix:   ", e.node_dir_prefix)
    println(io, spad, "keep_node_dirs:    ", e.keep_node_dirs)
    println(io, spad, "on_error:          ", e.on_error)
    println(io, spad, "dir_name_mode:     ", e.dir_name_mode)
    println(io, spad, "index_pad:         ", e.index_pad)
    println(io, spad, "value_sig:         ", e.value_sig)
end

"""
    copy_with(e::PriorExecutionOptions; kwargs...) -> PriorExecutionOptions

Return a copy of `e` with any provided keyword overrides applied. Valid keys:
- `mode`, `chain_direction`, `chain_snapshot_par`
- `scheduler`, `max_workers`
- `workdir_layout`, `node_dir_prefix`, `keep_node_dirs`
- `on_error`
- `dir_name_mode`, `index_pad`, `value_sig`
"""
function copy_with(e::PriorExecutionOptions; kwargs...)
    mode              = get(kwargs, :mode,              e.mode)
    chain_direction   = get(kwargs, :chain_direction,   e.chain_direction)
    chain_snapshot_par= get(kwargs, :chain_snapshot_par,e.chain_snapshot_par)
    scheduler         = get(kwargs, :scheduler,         e.scheduler)
    max_workers       = get(kwargs, :max_workers,       e.max_workers)
    workdir_layout    = get(kwargs, :workdir_layout,    e.workdir_layout)
    node_dir_prefix   = get(kwargs, :node_dir_prefix,   e.node_dir_prefix)
    keep_node_dirs    = get(kwargs, :keep_node_dirs,    e.keep_node_dirs)
    on_error          = get(kwargs, :on_error,          e.on_error)
    dir_name_mode     = get(kwargs, :dir_name_mode,     e.dir_name_mode)
    index_pad         = get(kwargs, :index_pad,         e.index_pad)
    value_sig         = get(kwargs, :value_sig,         e.value_sig)

    return PriorExecutionOptions(; mode, chain_direction, chain_snapshot_par,
        scheduler, max_workers, workdir_layout, node_dir_prefix, keep_node_dirs,
        on_error, dir_name_mode, index_pad, value_sig)
end

Base.copy(e::PriorExecutionOptions) = copy_with(e)

# --------------------------------------------------------------------------------------------------------------
# Naming & metadata helpers
# --------------------------------------------------------------------------------------------------------------

"""
    format_short(x; sig=6) -> String

Compact scientific formatting for real numbers. Preserves `NaN`/`Inf` in a readable form.
"""
format_short(x::Real; sig::Int=6) = isnan(x)  ? "NaN" :
                                    isinf(x)  ? (x > 0 ? "Inf" : "-Inf") :
                                    @sprintf("%.*g", sig, float(x))

"""
    sanitize_name(s) -> String

Make a string safe for filesystem by keeping only alphanumerics plus `._=-+`,
replacing the rest with `_`.
"""
sanitize_name(s::AbstractString) = replace(s, r"[^A-Za-z0-9._=\-\+]" => "_")

"""
    build_node_dirname(prefix, idx; pad=3, param=:θ, theta=nothing, mode=:index_only, sig=6) -> String

Build a per-node directory name according to the selected naming mode.

- `:index_only` → uses only an index with zero-padding.
- `:with_value` → additionally appends `__<PARAM>=<VALUE>` where `<VALUE>` is compact-formatted.
"""
function build_node_dirname(prefix::AbstractString, idx::Integer;
    pad::Int=3, param::Symbol=:theta, theta::Union{Nothing,Real}=nothing,
    mode::Symbol=:index_only, sig::Int=6)

    base = string(prefix, lpad(idx, pad, '0'))
    if mode === :index_only || theta === nothing
        return base
    elseif mode === :with_value
        val = sanitize_name(format_short(theta; sig=sig))
        par = sanitize_name(String(param))
        return string(base, "__", par, "=", val)
    else
        error("build_node_dirname: unknown dir_name_mode = $mode (expected :index_only or :with_value)")
    end
end

"""
    write_node_metadata(path; meta...) -> Nothing

Write a simple `key: value` metadata file at `path` (e.g., `"node_meta.txt"`).
- `Real` values are formatted via `format_short`.
- `AbstractVector{<:Real}` values are comma-joined with each element formatted via `format_short`.
- Other vectors are comma-joined via `string`.
- Other scalars are stringified.
"""
function write_node_metadata(path::AbstractString; meta...)
    open(path, "w") do io
        for (k, v) in pairs(meta)
            key = string(k)
            if v isa AbstractVector{<:Real}
                println(io, key, ": ", join(format_short.(v; sig=12), ","))
            elseif v isa AbstractVector
                println(io, key, ": ", join(string.(v), ","))
            elseif v isa Real
                println(io, key, ": ", format_short(v; sig=12))
            else
                println(io, key, ": ", string(v))
            end
        end
    end
    nothing
end

# --------------------------------------------------------------------------------------------------------------
# Per-node seed specification (starting par-files for prior nodes)
# --------------------------------------------------------------------------------------------------------------

abstract type AbstractNodeSeedSpec end

"""
    NoSeeds()

Disable per-node seeding (default).
"""
struct NoSeeds <: AbstractNodeSeedSpec end

"""
    SeedPaths(paths::Vector{String})

Explicit per-node par-file paths mapped by **node index**. The length must
match the number of nodes produced by the `nodes` rule.
"""
struct SeedPaths <: AbstractNodeSeedSpec
    paths::Vector{String}
end

# Pretty printer for SeedPaths
function Base.show(io::IO, ::MIME"text/plain", s::SeedPaths)
    indent = get(io, :indent, 0)
    pad, spad = repeat(" ", indent), repeat(" ", indent + 2)
    n = length(s.paths)
    println(io, pad, "SeedPaths")
    println(io, spad, "count: ", n)
    # Preview a few paths for convenience
    preview = min(n, 5)
    for i in 1:preview
        println(io, spad, lpad(string(i), 3), ": ", s.paths[i])
    end
    if n > preview
        println(io, spad, "… (", n - preview, " more)")
    end
end

"""
    SeedPaths(res::GeneralTempoResult) -> SeedPaths

Convenience: extract per-node `par_out_path` from a prior run's result and build
`SeedPaths` in the same order as `res.subresults`.
"""
function SeedPaths(res::GeneralTempoResult)
    n = length(res.subresults)
    n == 0 && error("SeedPaths(result): result has no subresults")
    paths = Vector{String}(undef, n)
    for (i, r) in enumerate(res.subresults)
        p = get(r.metadata, :par_out_path, nothing)
        p isa AbstractString || error("SeedPaths(result): subresult $i has no :par_out_path in metadata")
        paths[i] = String(p)
    end
    return SeedPaths(paths)
end

# --------------------------------------------------------------------------------------------------------------
# Settings that fully describe a prior-marginalized single-task
# --------------------------------------------------------------------------------------------------------------

const _PIN_MODES       = (:force, :fixed, :fit)
const _REF_STRATS      = (:prior_median, :grid_argmin, :spline_argmin, :custom_value)
const _REPRESENTATIVES = (:prior_median, :grid_argmin)

"""
    tempo_flag(pin_mode::Symbol) -> Int

Map a high-level pin mode to the TEMPO flag:
- `:force` → `-1`  (force this value; do not recompute even if derived)
- `:fixed` → `0`   (set, but do not fit)
- `:fit`   → `1`   (fit this parameter)
"""
@inline tempo_flag(pin_mode::Symbol) =
    pin_mode === :force ? -1 :
    pin_mode === :fixed ?  0 :
    pin_mode === :fit   ?  1 :
    error("Unknown pin_mode: $pin_mode")

"""
    PriorMarginalizationSettings{P,N}

Settings describing prior-based marginalization over a single parameter.

Fields
- `parameter`          : TEMPO parameter name (e.g., `:GAMMA`, `:XPBDOT`)
- `pin_mode`           : how to set the parameter per node (`:force` | `:fixed` | `:fit`)
- `prior`              : prior specification (analytic / grid / sampled)
- `nodes`              : rule to generate θ-nodes (mapped via prior inverse CDF, if applicable)
- `likelihood_source`  : metric key used to build likelihood (a key in `result.metrics`)
- `ref_strategy`       : reference choice for Δ (e.g., `:prior_median` | `:grid_argmin` | `:spline_argmin` | `:custom_value`)
- `ref_value`          : custom reference value (required iff `ref_strategy == :custom_value`)
- `representative`     : which θ to report at top-level (e.g., `:prior_median`)
- `save_node_results`  : whether to keep `subresults` (per-node `GeneralTempoResult`)
- `exec_options`       : execution policy (chaining, scheduler, dirs, error handling, naming)
- `metrics_hook`       : optional callback to extend the metrics dictionary
- `seed_spec`          : per-node seed specification (par-file inputs for nodes)

Notes
- `metrics_hook` is expected to have signature like:
  `(final_it::InternalIterationResult, conv::ConvergenceInfo, metrics::Dict{Symbol,Float64}) -> Nothing`
  and can mutate `metrics` in-place to add custom scalars.
"""
struct PriorMarginalizationSettings{P<:AbstractPriorSpec, N<:AbstractNodeRule}
    parameter::Symbol
    pin_mode::Symbol
    prior::P
    nodes::N

    likelihood_source::Symbol
    ref_strategy::Symbol
    ref_value::Union{Nothing,Float64}

    representative::Symbol
    save_node_results::Bool

    exec_options::PriorExecutionOptions
    metrics_hook::Union{Nothing,Function}
    seed_spec::AbstractNodeSeedSpec

    function PriorMarginalizationSettings(;
        parameter::Symbol,
        pin_mode::Symbol,
        prior::P,
        nodes::N,
        likelihood_source::Symbol = :chi2_fit,
        ref_strategy::Symbol      = :prior_median,  # also supports :grid_argmin | :spline_argmin | :custom_value
        ref_value::Union{Nothing,Real} = nothing,
        representative::Symbol    = :prior_median,  # e.g. :prior_median | :grid_argmin
        save_node_results::Bool   = true,
        exec_options::PriorExecutionOptions = PriorExecutionOptions(),
        metrics_hook::Union{Nothing,Function} = nothing,
        seed_spec::AbstractNodeSeedSpec = NoSeeds(),
    ) where {P<:AbstractPriorSpec, N<:AbstractNodeRule}

        validate_choice("pin_mode",       pin_mode,       _PIN_MODES)
        validate_choice("ref_strategy",   ref_strategy,   _REF_STRATS)
        validate_choice("representative", representative, _REPRESENTATIVES)

        # Custom reference: require finite value and sanity-check support
        local ref_val64::Union{Nothing,Float64} =
            ref_value === nothing ? nothing : Float64(ref_value)

        if ref_strategy === :custom_value
            (ref_val64 !== nothing && isfinite(ref_val64)) ||
                error("ref_strategy=:custom_value requires a finite ref_value")
            # Check against prior support (best-effort)
            lo, hi = prior_support(prior)
            if isfinite(lo) && isfinite(hi)
                (ref_val64 >= lo && ref_val64 <= hi) ||
                    error("ref_value=$ref_val64 is outside prior support [$lo, $hi]")
            end
        end

        # Ensure node list is non-empty (materializes sampled prior if needed)
        theta = eval_nodes(prior, nodes)
        isempty(theta) && error("nodes rule produced an empty set of theta values")

        new{P,N}(
            parameter,
            pin_mode,
            prior,
            nodes,
            likelihood_source,
            ref_strategy,
            ref_val64,
            representative,
            save_node_results,
            exec_options,
            metrics_hook,
            seed_spec,
        )
    end
end


# Pretty-print with indentation (compact, avoids recomputing nodes)
function Base.show(io::IO, ::MIME"text/plain", s::PriorMarginalizationSettings)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad, "PriorMarginalizationSettings")
    println(io, spad, "parameter:          ", s.parameter)
    println(io, spad, "pin_mode:           ", s.pin_mode)
    println(io, spad, "likelihood_source:  ", s.likelihood_source)
    if s.ref_strategy === :custom_value
        println(io, spad, "ref_strategy:       ", s.ref_strategy, " (ref_value=", s.ref_value, ")")
    else
        println(io, spad, "ref_strategy:       ", s.ref_strategy)
    end
    println(io, spad, "representative:     ", s.representative)
    println(io, spad, "save_node_results:  ", s.save_node_results)

    println(io, spad, "prior:")
    show(IOContext(io, :indent => indent + 4), MIME"text/plain"(), s.prior)

    println(io, spad, "nodes:")
    show(IOContext(io, :indent => indent + 4), MIME"text/plain"(), s.nodes)

    println(io, spad, "exec_options:")
    show(IOContext(io, :indent => indent + 4), MIME"text/plain"(), s.exec_options)

    if s.metrics_hook !== nothing
        println(io, spad, "metrics_hook:       provided")
    end

    println(io, spad, "seeds:")
    if s.seed_spec isa NoSeeds
        println(io, repeat(" ", indent + 4), "none")
    elseif s.seed_spec isa SeedPaths
        println(io, repeat(" ", indent + 4), "paths (", length((s.seed_spec::SeedPaths).paths), ")")
    else
        println(io, repeat(" ", indent + 4), string(typeof(s.seed_spec)))
    end
end

"""
    copy_with(s::PriorMarginalizationSettings; kwargs...) -> PriorMarginalizationSettings

Create a modified copy of `s`, overriding any subset of fields via keyword arguments.
All validations are performed by the regular constructor.

Keyword arguments mirror the settings fields:
- `parameter`, `pin_mode`, `prior`, `nodes`
- `likelihood_source`, `ref_strategy`, `ref_value`
- `representative`, `save_node_results`
- `exec_options`, `metrics_hook`, `seed_spec`
"""
function copy_with(s::PriorMarginalizationSettings; kwargs...)
    # pull values with defaults from `s`
    p_parameter        = get(kwargs, :parameter,        s.parameter)
    p_pin_mode         = get(kwargs, :pin_mode,         s.pin_mode)
    p_prior            = get(kwargs, :prior,            s.prior)
    p_nodes            = get(kwargs, :nodes,            s.nodes)
    p_like_src         = get(kwargs, :likelihood_source, s.likelihood_source)
    p_ref_strategy     = get(kwargs, :ref_strategy,     s.ref_strategy)
    p_ref_value        = get(kwargs, :ref_value,        s.ref_value)
    p_representative   = get(kwargs, :representative,   s.representative)
    p_save_nodes       = get(kwargs, :save_node_results, s.save_node_results)
    p_exec_opts        = get(kwargs, :exec_options,     s.exec_options)
    p_metrics_hook     = get(kwargs, :metrics_hook,     s.metrics_hook)
    p_seed_spec        = get(kwargs, :seed_spec,        s.seed_spec)

    return PriorMarginalizationSettings(
        parameter        = p_parameter,
        pin_mode         = p_pin_mode,
        prior            = p_prior,
        nodes            = p_nodes,
        likelihood_source= p_like_src,
        ref_strategy     = p_ref_strategy,
        ref_value        = p_ref_value,
        representative   = p_representative,
        save_node_results= p_save_nodes,
        exec_options     = p_exec_opts,
        metrics_hook     = p_metrics_hook,
        seed_spec        = p_seed_spec,
    )
end

"""
    Base.copy(s::PriorMarginalizationSettings)

Return an identical copy of `s`.
"""
Base.copy(s::PriorMarginalizationSettings) = copy_with(s)

# Pick a reference θ that depends only on the prior (no node info required).
# Supported here: :prior_median, :custom_value
function _ref_theta_prior_only(s::PriorMarginalizationSettings)::Float64
    pr = s.prior
    if s.ref_strategy === :custom_value
        v = s.ref_value
        v === nothing && error("ref_strategy=:custom_value requires settings.ref_value")
        lo, hi = prior_support(pr)
        return (isfinite(lo) && isfinite(hi)) ? clamp(v::Float64, lo, hi) : (v::Float64)
    elseif s.ref_strategy === :prior_median
        return prior_median(pr)
    else
        error("ref_strategy=$(s.ref_strategy) is not a prior-only choice")
    end
end

"""
    resolve_node_seed_paths(settings, thetas) -> Vector{Union{Nothing,String}}

Materialize `settings.seed_spec` into a per-node list of par paths (or `nothing`).
For `NoSeeds()` returns a vector of `nothing`s. For `SeedPaths`, validates length.
"""

function resolve_node_seed_paths(s::PriorMarginalizationSettings,
                                 thetas::Vector{Float64})::Vector{Union{Nothing,String}}
    N = length(thetas)
    spec = getfield(s, :seed_spec)
    if spec isa NoSeeds
        return fill(nothing, N)
    elseif spec isa SeedPaths
        paths = (spec::SeedPaths).paths
        length(paths) == N || error("SeedPaths: expected $N paths, got $(length(paths))")
        # return a fresh copy of strings to avoid accidental mutation
        return String.(paths)
    else
        error("Unknown seed spec $(typeof(spec))")
    end
end
