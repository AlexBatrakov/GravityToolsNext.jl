using Test
using GravityToolsNext

@testset "GravityToolsNext" begin
    @testset "Supported surface" begin
        include("package_surface.jl")
    end

    @testset "Wrapper surface" begin
        include("wrapper_surface.jl")
    end

    @testset "TempoCore fixtures" begin
        include("tempo_core_fixtures.jl")
    end

    @testset "TIM documents" begin
        include("tempo_tim_file.jl")
    end

    @testset "Observation preparation" begin
        include("timing_sources.jl")
    end

    @testset "Engine deadlines" begin
        include("tempo_deadlines.jl")
    end

    @testset "Result and task assembly" begin
        include("result_task_assembly.jl")
    end

    @testset "Runner TIM integrity" begin
        include("runner_tim_integrity.jl")
    end

    @testset "Gated integration" begin
        include("gated_integration.jl")
    end

    @testset "AdaptiveGridFramework" begin
        include("adaptive_gridframework.jl")
    end
end
