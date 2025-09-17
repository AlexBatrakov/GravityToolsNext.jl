# src/TempoFramework/MultiPointTasks/Adaptive2DGridTask.jl
# ----------------------------------------------------------------------------
# Adaptive2DGridTask
# A thin adapter that runs an arbitrary SingleTempoTask on an adaptive 2D grid.
# It wires your grid engine (AdaptiveRefinement2DGrid) to a `target_function`
# that:
#   • clones the base SingleTempoTask for a point (x, y) via `task_copy_with`,
#   • runs it,
#   • returns a NamedTuple with the requested scalars for the grid.
#
# Notes
# -----
# * We keep isolation per grid point by running each point in its own job directory (`work_mode=:jobdir`, distinct `job_name`) and using a per-point par_output stem ...
# * We do not assume any particular SingleTempoTask internals beyond the public
#   helpers (task_copy_with, task_derive_par_output) and GeneralTempoResult API.
# ----------------------------------------------------------------------------

# -- Options ------------------------------------------------------------------

Base.@kwdef struct GridWorkspaceOptions
    tag_mode::Symbol = :with_value           # :index_only | :with_value
    value_sig::Int = 6                       # significant digits for values in tags
    point_job_prefix::String = "grid_points/"       # job_name prefix per point
    point_overwrite::Symbol = :reuse         # :error | :reuse | :unique | :clean
    point_layout::Symbol = :split            # :flat | :split (inside each job dir)
    io_mirror::Union{Symbol,Int,Tuple{Symbol,Int}} = :none
    snapshot_par::Bool = true
    link_tim::Bool = false
    keep_tmp_on_success::Bool = false
    save_result_jld2::Bool = false
    results_dir::String = "results"
end

# -- Task ---------------------------------------------------------------------

struct Adaptive2DGridTask{T1<:SingleTempoTask, RX<:AbstractGridRule, RY<:AbstractGridRule, RU<:Tuple{Vararg{AbstractRefinementUnit}}} <: MultiPointTask
    base_task::T1
    x::GridAxis{RX}
    y::GridAxis{RY}
    ref_settings::RefinementSettings{RU}
    opts::GridWorkspaceOptions
end

function Adaptive2DGridTask(;
    base_task::T1,
    x::GridAxis{RX},
    y::GridAxis{RY},
    ref_settings::RefinementSettings{RU},
    opts::GridWorkspaceOptions=GridWorkspaceOptions(),
) where {T1<:SingleTempoTask, RX<:AbstractGridRule, RY<:AbstractGridRule, RU<:Tuple{Vararg{AbstractRefinementUnit}}}
    return Adaptive2DGridTask{T1,RX,RY,RU}(base_task, x, y, ref_settings, opts)
end

# -- Helpers ------------------------------------------------------------------

# compact numeric to string for tags
_format_val(x::Real; sig::Int=6) = format_short(x; sig=sig)

# Build a stable, filesystem-safe tag for a point
function _point_tag(axX::GridAxis, axY::GridAxis, xv::Real, yv::Real, opts::GridWorkspaceOptions)
    if opts.tag_mode === :index_only
        ix = findfirst(==(xv), axis_values(axX))::Union{Int,Nothing}
        iy = findfirst(==(yv), axis_values(axY))::Union{Int,Nothing}
        # Fallback if exact match is not present (explicit grid may not be identity)
        ix === nothing && (ix = 1)
        iy === nothing && (iy = 1)
        return string(Symbol(axX.name), lpad(ix, 3, '0'), "_",
                      Symbol(axY.name), lpad(iy, 3, '0'))
    elseif opts.tag_mode === :with_value
        xs = sanitize_name(_format_val(xv; sig=opts.value_sig))
        ys = sanitize_name(_format_val(yv; sig=opts.value_sig))
        return string(Symbol(axX.name), "=", xs, "_", Symbol(axY.name), "=", ys)
    else
        error("GridWorkspaceOptions.tag_mode must be :index_only or :with_value")
    end
end

# Derive per-point workspace fields from base task (job-per-point)
function _derive_point_workspace(base_task::SingleTempoTask, tag::AbstractString, opts::GridWorkspaceOptions)
    job_name = string(opts.point_job_prefix, tag)
    par_out  = task_derive_par_output(base_task, tag)
    return (; work_mode = :jobdir, job_name, par_output = par_out, layout = opts.point_layout)
end

# Extract a NamedTuple of requested keys from GeneralTempoResult
function _extract_namedtuple(res::GeneralTempoResult, keys::Tuple{Vararg{Symbol}})
    vals = Vector{Float64}(undef, length(keys))
    for (i, k) in enumerate(keys)
        v = NaN
        # 1) metrics/extras
        if hasproperty(res, :metrics) && res.metrics !== nothing && haskey(res.metrics, k)
            v = try
                Float64(res.metrics[k])
            catch
                NaN
            end
        elseif hasproperty(res, :extras) && res.extras !== nothing && haskey(res.extras, k)
            v = try
                Float64(res.extras[k])
            catch
                NaN
            end
        # 2) fitted parameter estimates
        elseif hasproperty(res, :param_estimates) && res.param_estimates !== nothing && haskey(res.param_estimates, k)
            v = try
                Float64(res.param_estimates[k].value)
            catch
                NaN
            end
        end
        vals[i] = v
    end
    return NamedTuple{keys}(Tuple(vals))
end

# Save full point result if requested
function _maybe_save_result(res::GeneralTempoResult, tag::AbstractString, opts::GridWorkspaceOptions)
    opts.save_result_jld2 || save_result_jld2(res; filename=string(tag, ".jld2"))
end

# -- Execution ----------------------------------------------------------------

function run_task(task::Adaptive2DGridTask)
    # 0) Build initial grid holder
    grid_init = AdaptiveRefinement2DGrid(task.x, task.y, task.ref_settings)

    # 1) Construct per-point target function for the grid engine
    function target_function(xv, yv; task::Adaptive2DGridTask=task)
        # Format per-point overrides
        overrides = [
            TP(String(task.x.name), xv),
            TP(String(task.y.name), yv),
        ]

        # Build tag, job_name and par_output (job-per-point isolation)
        tag = _point_tag(task.x, task.y, xv, yv, task.opts)
        ws  = _derive_point_workspace(task.base_task, tag, task.opts)

        # Clone base SingleTempoTask with point-specific workspace & overrides
        point_task = task_copy_with(task.base_task;
            override_params_upsert = overrides,
            work_mode  = ws.work_mode,
            job_name   = ws.job_name,
            par_output = ws.par_output,
            layout     = ws.layout,
            io_mirror  = task.opts.io_mirror,
            snapshot_par = task.opts.snapshot_par,
            link_tim     = task.opts.link_tim,
            keep_tmp_on_success = task.opts.keep_tmp_on_success,
            overwrite = task.opts.point_overwrite,
        )

        # Execute and map to NamedTuple required by the grid engine
        println("Run started:  $(task.x.name) = $xv, $(task.y.name) = $yv")
        res = run_task(point_task)

        println("Run finished: $(task.x.name) = $xv, $(task.y.name) = $yv; chisqr = $(res.metrics[:chi2_fit]), pre_post = $(res.metrics[:pre_post_final])")

        # Save full result if requested
        _maybe_save_result(res, tag, task.opts)

        return _extract_namedtuple(res, task.ref_settings.params_to_save)
    end

    # 2) params_function! hook (optional aggregation per refinement step)
    function params_function!(grid::AdaptiveRefinement2DGrid)
        # no-op for now; the engine may call this after each refinement wave
        return nothing
    end

    # 3) Run the adaptive grid calculation (your engine owns the loop)
    grid_refined = calculate_2DGrid(grid_init, target_function, params_function!)

    return grid_refined
end
