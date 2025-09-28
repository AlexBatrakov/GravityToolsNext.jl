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
end

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