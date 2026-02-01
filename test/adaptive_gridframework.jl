using Test
using GravityToolsNext

@testset "GridAxis" begin
    ax = GridAxis(:x; min=-1.0, max=1.0, N=5, rule=LinRule())
    @test axisvalues(ax) == collect(LinRange(-1.0, 1.0, 5))
    @test linspace(ax) == axisvalues(ax)

    axlog = GridAxis(:a; min=-1e-4, max=-1e-2, N=3, rule=LogRule())
    v = axisvalues(axlog)
    @test all(v .< 0)
    @test isapprox(abs(v[1]), 1e-4; rtol=0, atol=0)
    @test isapprox(abs(v[2]), 1e-3; rtol=1e-12, atol=0)
    @test isapprox(abs(v[3]), 1e-2; rtol=0, atol=1e-15)

    axexp = GridAxis(:θ, [0.0, 1.0, 2.0])
    axexp2 = refine(axexp)
    @test axisvalues(axexp2) == [0.0, 0.5, 1.0, 1.5, 2.0]

    ax2 = refine(ax, 2)
    @test ax2.N == 4 * ax.N - 3
end

@testset "AdaptiveRefinement2DGrid precalculate + selectors" begin
    x = GridAxis(:x; min=-1.0, max=1.0, N=3, rule=LinRule())
    y = GridAxis(:y; min=-2.0, max=2.0, N=3, rule=LinRule())

    ref_sets = RefinementSettings(
        FullUnit(:z);
        desired_refinement_level=0,
        parallel=false,
        params_to_save=(:z,)
    )

    grid = AdaptiveRefinement2DGrid(x, y, ref_sets)

    target_function(xv, yv) = (z = xv^2 + yv^2,)
    params_function!(g) = nothing

    precalculate_2DGrid!(grid, target_function, params_function!)

    @test haskey(grid.vars, :z)
    @test size(grid.vars[:z]) == (3, 3)
    @test grid.status == fill(1, 3, 3)
    @test isapprox(grid.min[:z], minimum(grid.vars[:z]); rtol=0, atol=0)
    @test isapprox(grid.max[:z], maximum(grid.vars[:z]); rtol=0, atol=0)

    @test cell_selector(1, 1, grid) == true
    @test cell_selector(0, 1, grid) == false
    @test cell_selector(1, 0, grid) == false
    @test cell_selector(3, 1, grid) == false  # i_cell must be < x.N

    grid2 = refine(grid)
    @test grid2.x.N == 5
    @test grid2.y.N == 5
    @test size(grid2.ref_level) == (5, 5)
    @test grid2.status[1, 1] == 1
    @test grid2.status[2, 2] == -1
end
