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

    @testset "Result and task assembly" begin
        include("result_task_assembly.jl")
    end

    @testset "Gated integration" begin
        include("gated_integration.jl")
    end

    @testset "AdaptiveGridFramework" begin
        include("adaptive_gridframework.jl")
    end
end
