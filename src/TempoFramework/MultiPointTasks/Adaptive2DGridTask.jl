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
    grid_root::String                          # absolute OR relative to base task work_dir
    point_job_prefix::String = "grid_points"   # <grid_root>/<point_job_prefix>/<TAG>/ becomes job_name
    results_dirname::String   = "results"      # <grid_root>/<results_dirname>/<TAG>.jld2
    input_dirname::String     = "input"        # optional staging of base inputs once per grid

    tag_mode::Symbol = :with_value             # :with_value | :hash
    tag_value_sig::Int = 6                     # significant digits for values in tags
    tag_sep::String = "__"                     # separator between X and Y parts (for :with_value)
    sanitize::Bool = true                      # filesystem-safe tag strings
    hash_len::Int = 10                         # significant hex chars in hash (for :hash)

    save_results_jld2::Bool = true             # save full GeneralTempoResult per point
    stage_inputs::Bool = true                  # stage base inputs into <grid_root>/<input_dirname>
    stage_inputs_mode::Symbol = :root          # :root | :subdir  (where to place staged inputs)
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

# compact numeric to string for tags (unchanged)
_format_val(x::Real; sig::Int=6) = format_short(x; sig=sig)

# Build a stable, filesystem-safe tag for a point
function _point_tag(axX::GridAxis, axY::GridAxis, xv::Real, yv::Real, opts::GridWorkspaceOptions)
    if opts.tag_mode === :with_value
        xs = _format_val(xv; sig=opts.tag_value_sig)
        ys = _format_val(yv; sig=opts.tag_value_sig)
        tag = string(axX.name, "=", xs, opts.tag_sep, axY.name, "=", ys)
        return opts.sanitize ? sanitize_name(tag) : tag
    elseif opts.tag_mode === :hash
        # Stable textual encoding for floats (17 significant digits)
        xs = @sprintf("%.17g", Float64(xv))
        ys = @sprintf("%.17g", Float64(yv))
        h  = bytes2hex(sha1(String(codeunits(string(axX.name, "|", xs, "|", axY.name, "|", ys)))))[1:opts.hash_len]
        tag = string(axX.name, "_", axY.name, "_h", h)
        return opts.sanitize ? sanitize_name(tag) : tag
    else
        error("GridWorkspaceOptions.tag_mode must be :with_value or :hash")
    end
end

# Helper functions for paths
_grid_results_dir(root::AbstractString, o::GridWorkspaceOptions) = joinpath(root, o.results_dirname)
_grid_point_job_name(o::GridWorkspaceOptions, tag::AbstractString) = joinpath(o.point_job_prefix, tag)
_grid_point_result_path(root::AbstractString, o::GridWorkspaceOptions, tag::AbstractString) = joinpath(_grid_results_dir(root, o), string(tag, ".jld2"))

# Derive per-point workspace fields from base task (job-per-point)
function _derive_point_workspace(base_task::SingleTempoTask, tag::AbstractString, root_dir::AbstractString, opts::GridWorkspaceOptions)
    job_name = _grid_point_job_name(opts, tag)
    par_out  = task_derive_par_output(base_task, tag)
    return (; work_dir = root_dir, job_name, par_output = par_out)
end

# Extract a NamedTuple of requested keys from GeneralTempoResult
function _extract_namedtuple(res::GeneralTempoResult, keys::Tuple{Vararg{Symbol}})
    vals = Vector{Float64}(undef, length(keys))
    @inbounds for (i, k) in enumerate(keys)
        v = NaN
        if !isnothing(res.metrics) && haskey(res.metrics, k)
            v = try
                Float64(res.metrics[k])
            catch
                NaN
            end
        elseif !isnothing(res.param_estimates) && haskey(res.param_estimates, k)
            v = try
                Float64(res.param_estimates[k][1])  # take value from (value, sigma)
            catch
                NaN
            end
        end
        vals[i] = v
    end
    return NamedTuple{keys}(Tuple(vals))
end

# Save full point result if requested
function _maybe_save_result(res::GeneralTempoResult, tag::AbstractString, root_dir::AbstractString, opts::GridWorkspaceOptions)
    opts.save_results_jld2 || return
    mkpath(_grid_results_dir(root_dir, opts))
    save_result_jld2(res; filename = _grid_point_result_path(root_dir, opts, tag))
end

function _save_grid_result(grid::AdaptiveRefinement2DGrid)
    jldsave(filename; grid=grid)
    return nothing
end

# -- Execution ----------------------------------------------------------------

function run_task(task::Adaptive2DGridTask; grid_init = nothing, just_refine = false)::AdaptiveRefinement2DGrid
    # Prepare workspace
    base_root = task_workdir(task.base_task)
    grid_root_abs = isabspath(task.opts.grid_root) ? task.opts.grid_root : joinpath(base_root, task.opts.grid_root)
    mkpath(grid_root_abs)
    if task.opts.stage_inputs
        dest_dir = task.opts.stage_inputs_mode === :root ? grid_root_abs : joinpath(grid_root_abs, task.opts.input_dirname)
        mkpath(dest_dir)
        task_stage_inputs!(task.base_task, dest_dir)
    end

    # 0) Build initial grid holder
    grid_init = isnothing(grid_init) ? AdaptiveRefinement2DGrid(task.x, task.y, task.ref_settings) : grid_init
    keys_tuple = Tuple(task.ref_settings.params_to_save)

    # 1) Construct per-point target function for the grid engine
    function target_function(xv, yv; task::Adaptive2DGridTask=task)
        overrides = [TP(String(task.x.name), xv), TP(String(task.y.name), yv)]
        tag = _point_tag(task.x, task.y, xv, yv, task.opts)
        ws  = _derive_point_workspace(task.base_task, tag, grid_root_abs, task.opts)

        point_task = task_copy_with(task.base_task;
            par_output = ws.par_output,
            override_params_upsert = overrides,
            work_mode = :jobdir,
            job_name  = ws.job_name,
            overwrite = :clean,
            work_dir  = ws.work_dir,
        )

        param_name = task.ref_settings.units[1].name

        println("Run started:  $(task.x.name) = $xv, $(task.y.name) = $yv")
        res = run_task(point_task)
        println("Run finished: $(task.x.name) = $xv, $(task.y.name) = $yv;\n
            $param_name = $(get(res.metrics, param_name, NaN)), pre_post_final = $(get(res.metrics, :pre_post_final, NaN))")

        _maybe_save_result(res, tag, grid_root_abs, task.opts)
        return _extract_namedtuple(res, keys_tuple)
    end

    # 2) params_function! hook (optional aggregation per refinement step)
    function params_function!(grid::AdaptiveRefinement2DGrid)
        # no-op for now; the engine may call this after each refinement wave
        return nothing
    end

    # 3) Run the adaptive grid calculation (your engine owns the loop)
    if just_refine
        grid_refined = refine_2DGrid(grid_init, target_function, params_function!)
    else
        grid_refined = calculate_2DGrid(grid_init, target_function, params_function!)
    end

    return grid_refined
end
