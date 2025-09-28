# src/TempoFramework/SingleTasks/IterativeTempoTask.jl

struct IterativeSettings
    max_iterations::Int           # maximum number of iterations
    tol_chi2::Float64             # convergence tolerance on chi2 (relative)
    tol_params::Float64           # convergence tolerance on parameters (absolute)
    min_iterations::Int = 1       # minimum number of iterations
    verbose::Bool = true          # print iteration info to stdout
end

"""
    IterativeTempoTask(base_task::SingleTempoTask, settings::IterativeSettings)
        -> IterativeTempoTask   
Wrap a `SingleTempoTask` to run it iteratively, updating the input parameters
after each iteration based on the output best-fit values, until convergence.
The `settings` control the maximum number of iterations and convergence criteria.
"""
struct IterativeTempoTask{T1<:SingleTempoTask} <: SingleTempoTask
    base_task::T1
    settings::IterativeSettings
end



function run_task(task::IterativeTempoTask)::GeneralTempoResult

end

function task_workdir(task::IterativeTempoTask)
    return task.base_task |> task_workdir
end

function task_stage_inputs!(task::IterativeTempoTask, dest_dir::AbstractString)
    return task.base_task |> task_stage_inputs!(dest_dir)
end