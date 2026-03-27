# TempoFramework/TempoCore/TempoTasks.jl

# Abstract task hierarchy
abstract type AbstractTempoTask end
abstract type SingleTempoTask  <: AbstractTempoTask end
abstract type MultiPointTask   <: AbstractTempoTask end

# -----------------------------------------------------------------------------
# Mandatory API for single tasks
# -----------------------------------------------------------------------------

"""
    run_task(task::SingleTempoTask) -> GeneralTempoResult

Execute the task and return a unified result.
Every concrete `SingleTempoTask` must provide a method.
"""
function run_task(task::SingleTempoTask)::GeneralTempoResult
    error("run_task not implemented for $(typeof(task))")
end

"""
    task_workdir(task::SingleTempoTask) -> AbstractString

Return the working directory where the task runs / writes results.
Used by higher-level wrappers to organize per-node subdirectories.
"""
function task_workdir(task::SingleTempoTask)::AbstractString
    error("task_workdir not implemented for $(typeof(task))")
end

# -----------------------------------------------------------------------------
# Parameter override helpers for building per-node tasks
# -----------------------------------------------------------------------------

"""
    task_with_overrides(task::T, overrides::AbstractVector{TempoParameter};
                        work_dir::AbstractString) -> T where {T<:SingleTempoTask}

Return a copy of `task` (same concrete type) with additional/overriding
parameters `overrides` and a reassigned working directory `work_dir`.

Implementations must merge `overrides` with any existing parameter overrides
present in the task's settings, giving precedence to the new ones.
"""
function task_with_overrides(task::T,
                             overrides::AbstractVector{TempoParameter};
                             work_dir::AbstractString)::T where {T<:SingleTempoTask}
    error("task_with_overrides not implemented for $(T)")
end

"""
    task_with_param(task::T, name::Symbol, value::Float64, flag::Int;
                    work_dir::AbstractString) -> T where {T<:SingleTempoTask}

Convenience wrapper over `task_with_overrides` for a single parameter.
By default constructs `TP(String(name), value; flag=flag)` and delegates.
Usually there is no need to override this method.
"""
function task_with_param(task::T, name::Symbol, value::Float64, flag::Int;
                         work_dir::AbstractString)::T where {T<:SingleTempoTask}
    return task_with_overrides(task, [TP(String(name), value; flag=flag)];
                               work_dir=work_dir)
end

# -----------------------------------------------------------------------------
# Optional: stage inputs for a target workdir
# -----------------------------------------------------------------------------

"""
    task_stage_inputs!(task::SingleTempoTask, dest_dir::AbstractString) -> Nothing

Stage/copy any task-specific input artifacts into `dest_dir` so that the task
can run there using name-only paths. Default implementation does nothing.

Wrappers like prior-marginalization should call this before cloning a task
with a new `work_dir`.
"""
function task_stage_inputs!(::SingleTempoTask, ::AbstractString)
    return nothing
end

"""
    task_copy_with(task::T; kwargs...) -> T where {T<:SingleTempoTask}

Extension hook used by higher-level wrappers to clone a task while preserving
its concrete type. Supported wrapper tasks rely on this hook when they need to
change workspace/layout details or upsert parameter overrides for derived runs.
"""
function task_copy_with(task::T; kwargs...)::T where {T<:SingleTempoTask}
    error("task_copy_with not implemented for $(typeof(task))")
end

# -----------------------------------------------------------------------------
# Extension hook: derive per-node par output filename
# -----------------------------------------------------------------------------

"""
    task_derive_par_output(task::SingleTempoTask, node_tag::AbstractString) -> String

Extension hook used by wrapper tasks to derive a stable per-node/par-output
filename stem from a base task.
"""
function task_derive_par_output(t::SingleTempoTask, node_tag::AbstractString)
    error("task_derive_par_output not implemented for $(typeof(t)) — provide a method for your task type")
end

"""
    save_result_jld2(result; filename::AbstractString) -> Nothing

Optional persistence hook used by wrappers that save per-node task results. Task
families that support result persistence should provide a concrete method.
"""
function save_result_jld2(result; filename::AbstractString)
    error("save_result_jld2 not implemented for $(typeof(result))")
end

# -----------------------------------------------------------------------------
# Multi-point tasks (if/when needed)
# -----------------------------------------------------------------------------

"""
    run_task(task::MultiPointTask)

Entry point for multi-point workflows. Define as needed for your project.
"""
function run_task(task::MultiPointTask)
    error("run_task not implemented for $(typeof(task))")
end

