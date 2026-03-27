using Test
using GravityToolsNext

@testset "GravityToolsNext" begin
    @testset "Supported surface" begin
        include("package_surface.jl")
    end

    @testset "AdaptiveGridFramework" begin
        include("adaptive_gridframework.jl")
    end
end
