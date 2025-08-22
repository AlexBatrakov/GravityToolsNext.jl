abstract type AbstractTempoTask end

abstract type SingleTempoTask <: AbstractTempoTask end

function run_task(task::SingleTempoTask)::GeneralTempoResult
    error("run_task not implemented for $(typeof(task))")
end

abstract type MultiPointTask  <: AbstractTempoTask end

function run_task(task::MultiPointTask)
    error("run_task not implemented for $(typeof(task))")
end