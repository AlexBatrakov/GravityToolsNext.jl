# src/AdaptiveGridFramework/AdaptiveRefinement2DGrid.jl

#--------------------------------------------------------------------------------------------------------------
# AdaptiveRefinement2DGrid.jl
# 
# This file provides the AdaptiveRefinement2DGrid type and associated routines for adaptive refinement
# of 2D grids. It includes types, constructors, cell selection logic for refinement, and both single-core
# and parallel routines for grid calculation and refinement, as well as helper functions.
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# Types & constructors
#--------------------------------------------------------------------------------------------------------------

"""
    AdaptiveRefinement2DGrid

Structure representing a 2D adaptive refinement grid with variables, parameters, axes, refinement settings, and refinement status.
"""
struct AdaptiveRefinement2DGrid{T1 <: AbstractGridRule, T2 <: AbstractGridRule, T3}
    vars::Dict{Symbol,Matrix{Float64}}
    params::Dict{Symbol,Float64}
    min::Dict{Symbol,Float64}
    max::Dict{Symbol,Float64}
    x::GridAxis{T1}
    y::GridAxis{T2}
    ref_sets::RefinementSettings{T3}
    ref_level::Matrix{Int64}
    status::Matrix{Int64}
end

"""
    AdaptiveRefinement2DGrid(x::GridAxis, y::GridAxis, ref_sets)

Construct a new AdaptiveRefinement2DGrid from grid axes and refinement settings.
"""
function AdaptiveRefinement2DGrid(x::GridAxis{T1}, y::GridAxis{T2}, ref_sets::T3) where {T1 <: AbstractGridRule, T2 <: AbstractGridRule, T3}
    vars = Dict{Symbol,Matrix{Float64}}()
    params = Dict{Symbol,Float64}()
    min = Dict{Symbol,Float64}()
    max = Dict{Symbol,Float64}()
    ref_level = [0 for i in 1:x.N, j in 1:y.N]
    status = [-1 for i in 1:x.N, j in 1:y.N]
    return AdaptiveRefinement2DGrid(vars, params, min, max, x, y, ref_sets, ref_level, status)
end

"""
    Base.show(io::IO, ::MIME"text/plain", grid::AdaptiveRefinement2DGrid)

Pretty-print the AdaptiveRefinement2DGrid with variable, parameter, axis, and refinement information.
"""
function Base.show(io::IO, ::MIME"text/plain", grid::AdaptiveRefinement2DGrid)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    iop    = IOContext(io, :indent => indent + 4)
    
    println(io, pad,  "AdaptiveRefinement2DGrid:")
	println(io, spad, "Variables: ", keys(grid.vars))
    println(io, spad, "Parameters: ", grid.params)
    println(io, spad, "Minimal values: ", grid.min)
    println(io, spad, "Maximal values: ", grid.max)
    println(io, spad, "X axis: ", grid.x)
    println(io, spad, "Y axis: ", grid.y)
    show(iop, MIME"text/plain"(), grid.ref_sets)
end

#--------------------------------------------------------------------------------------------------------------
# Cell selectors
#--------------------------------------------------------------------------------------------------------------

"""
    cell_selector(i_cell, j_cell, grid)

Determine if the cell at (i_cell, j_cell) should be refined according to the refinement settings.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid)
    if !((0 < i_cell < grid.x.N) && (0 < j_cell < grid.y.N))
        return false
    end
    combined_case = false
    for ref_unit in grid.ref_sets.units
        combined_case = combined_case || cell_selector(i_cell, j_cell, grid, ref_unit)
    end
    return combined_case::Bool
end

"""
    cell_selector(i_cell, j_cell, grid, ref_unit::LocalMinimaUnit; at_corner=false)

Selects cells containing a unique local minimum within a specified range.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::LocalMinimaUnit; at_corner=false)
    cell     = @view grid.vars[ref_unit.name][i_cell : min(i_cell+1, grid.x.N) , j_cell : min(j_cell+1, grid.y.N)]
    big_cell = @view grid.vars[ref_unit.name][max(i_cell-1, 1) : min(i_cell+2, grid.x.N), max(j_cell-1, 1) : min(j_cell+2, grid.y.N)]
    value_big_cell_min = minimum(big_cell)
    value_cell_min = minimum(cell)
    value_min = grid.min[ref_unit.name]

    # Check that the minimum in big_cell is within the cell
    local_minima_case = value_big_cell_min in cell

    # Check for uniqueness of the minimum within big_cell (not repeated outside cell)
    unique_minimum_case = sum(big_cell .== value_big_cell_min) - sum(cell .== value_big_cell_min) == 0

    # Optionally require the minimum to be at the cell corner
    at_corner_case = grid.vars[ref_unit.name][i_cell,j_cell] == value_big_cell_min

    if ref_unit.from_min
        min_case = ref_unit.min <= value_big_cell_min - value_min
        max_case = ref_unit.max >= value_big_cell_min - value_min
    else
        min_case = ref_unit.min <= value_big_cell_min
        max_case = ref_unit.max >= value_big_cell_min
    end

    return local_minima_case && unique_minimum_case && min_case && max_case && (!at_corner || at_corner_case)
end

"""
    cell_selector(i_cell, j_cell, grid, ref_unit::FullUnit)

Selects cells where the min/max value is within specified bounds.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::FullUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    return min_case && max_case
end

"""
    cell_selector(i_cell, j_cell, grid, ref_unit::DiffUnit)

Selects cells where the difference across the cell exceeds a threshold.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::DiffUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    value_min = grid.min[ref_unit.name]
    if ref_unit.from_min
        min_case = ref_unit.min <= value_cell_min - value_min
        max_case = ref_unit.max >= value_cell_max - value_min
    else
        min_case = ref_unit.min <= value_cell_min
        max_case = ref_unit.max >= value_cell_max
    end
    diff_case = value_cell_max - value_cell_min > ref_unit.diff
    return diff_case && min_case && max_case
end

"""
    cell_selector(i_cell, j_cell, grid, ref_unit::RelDiffUnit)

Selects cells where the relative difference across the cell exceeds a threshold.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::RelDiffUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    value_cell_mean = mean(cell)
    value_min = grid.min[ref_unit.name]
    if ref_unit.from_min
        min_case = ref_unit.min <= value_cell_min - value_min
        max_case = ref_unit.max >= value_cell_max - value_min
        rel_diff_case = (value_cell_max - value_cell_min)  > ref_unit.rel_diff * (value_cell_mean - value_min)
    else
        min_case = ref_unit.min <= value_cell_min
        max_case = ref_unit.max >= value_cell_max
        rel_diff_case = (value_cell_max - value_cell_min) > ref_unit.rel_diff * value_cell_mean
    end
    return rel_diff_case && min_case && max_case
end
    
"""
    cell_selector(i_cell, j_cell, grid, ref_unit::ContourUnit)

Selects cells where the cell crosses a specified contour value.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::ContourUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    value_min = grid.min[ref_unit.name]
    if ref_unit.from_min
        contour_case = any(value_cell_min .< value_min .+ ref_unit.contours .< value_cell_max)
        min_case = ref_unit.min <= value_cell_min - value_min
        max_case = ref_unit.max >= value_cell_max - value_min
    else
        contour_case = any(value_cell_min .< ref_unit.contours .< value_cell_max)
        min_case = ref_unit.min <= value_cell_min
        max_case = ref_unit.max >= value_cell_max
    end
    return contour_case && min_case && max_case
end

"""
    cell_selector(i_cell, j_cell, grid, ref_unit::DiffContourUnit)

Selects cells where the cell crosses a contour and the difference exceeds a threshold.
"""
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::DiffContourUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    value_min = grid.min[ref_unit.name]
    if ref_unit.from_min 
        diff_contour_case = any((value_cell_min .< value_min .+ ref_unit.contours .< value_cell_max) .&& (value_cell_max .- value_cell_min .> ref_unit.diffs))
        min_case = ref_unit.min <= value_cell_min - value_min
        max_case = ref_unit.max >= value_cell_max - value_min
    else
        diff_contour_case = any((value_cell_min .< ref_unit.contours .< value_cell_max) .&& (value_cell_max .- value_cell_min .> ref_unit.diffs))
        min_case = ref_unit.min <= value_cell_min
        max_case = ref_unit.max >= value_cell_max
    end
    return diff_contour_case && min_case && max_case
end

#--------------------------------------------------------------------------------------------------------------
# Refinement helpers
#--------------------------------------------------------------------------------------------------------------

"""
    refine(arr::Matrix)

Refine a matrix by doubling its resolution, inserting -1 in new positions.
"""
function refine(arr::Matrix{T}) where {T}
    arr_refined = fill(-one(T), 2 .* size(arr) .- 1)
    for i in 1:size(arr)[1], j in 1:size(arr)[2]
        arr_refined[2*i-1,2*j-1] = arr[i,j]
    end
    return arr_refined::Matrix{T}
end

"""
    refine(dict::Dict{Symbol,Matrix})

Refine each matrix in a dictionary of matrices.
"""
function refine(dict::Dict{Symbol,Matrix{T}}) where {T}
    dict_refined = Dict{Symbol,Matrix{T}}()
    for (key, value) in dict
        dict_refined[key] = refine(dict[key])
    end
    return dict_refined::Dict{Symbol,Matrix{T}}
end

"""
    refine(grid::AdaptiveRefinement2DGrid)

Return a refined version of the grid with doubled resolution.
"""
function refine(grid::AdaptiveRefinement2DGrid)
    vars_refined = refine(grid.vars)
    params_refined = copy(grid.params)
    min_refined = copy(grid.min)
    max_refined = copy(grid.max)
    x_refined = refine(grid.x)
    y_refined = refine(grid.y)
    ref_sets_refined = grid.ref_sets
    ref_level_refined = refine(grid.ref_level)
    status_refined = refine(grid.status)
    grid_refined = AdaptiveRefinement2DGrid(vars_refined, params_refined, min_refined, max_refined, x_refined, y_refined, ref_sets_refined, ref_level_refined, status_refined)
    return grid_refined
end

#--------------------------------------------------------------------------------------------------------------
# High-level API
#--------------------------------------------------------------------------------------------------------------

"""
    calculate_2DGrid(grid, target_function, params_function!)

Calculate and adaptively refine the grid up to the desired refinement level.
"""
function calculate_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    precalculate_2DGrid!(grid, target_function, params_function!)
    for i in 1:grid.ref_sets.desired_refinement_level
        grid = refine_2DGrid(grid, target_function, params_function!)
    end
    return grid
end

"""
    precalculate_2DGrid!(grid, target_function, params_function!)

Populate the grid with initial values before refinement, using parallel or single-core as needed.
"""
function precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    if grid.ref_sets.parallel == true
        return parallel_precalculate_2DGrid!(grid, target_function, params_function!)
    else 
        return single_core_precalculate_2DGrid!(grid, target_function, params_function!)
    end
end

"""
    refine_2DGrid(grid, target_function, params_function!)

Refine the grid by one level, using parallel or single-core as requested.
"""
function refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    if grid.ref_sets.parallel == true
        return parallel_refine_2DGrid(grid, target_function, params_function!)
    else 
        return single_core_refine_2DGrid(grid, target_function, params_function!)
    end
end

#--------------------------------------------------------------------------------------------------------------
# Single-core
#--------------------------------------------------------------------------------------------------------------

"""
    single_core_precalculate_2DGrid!(grid, target_function, params_function!)

Fill the grid with values using a single core.
"""
function single_core_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    # target_keys = Base.return_types(target_function, Tuple{Float64,Float64})[1].parameters[1]
    for key in grid.ref_sets.params_to_save
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end
    for i in 1:grid.x.N, j in 1:grid.y.N
        target_output = target_function(grid.x.values[i], grid.y.values[j])
        for (key, value) in pairs(target_output)
            grid.vars[key][i, j] = value
        end
    end
    for key in grid.ref_sets.params_to_save
        grid.min[key] = minimum(x->isnan(x) ? +Inf : x, grid.vars[key])
        grid.max[key] = maximum(x->isnan(x) ? -Inf : x, grid.vars[key])
    end
    grid.status .= 1
    params_function!(grid)
    return grid
end

"""
    single_core_refine_2DGrid(grid, target_function, params_function!)

Refine the grid by one level using a single core.
"""
function single_core_refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    grid_refined = refine(grid)

    println("\nRefinement from ($(grid.x.N), $(grid.y.N)) to ($(grid_refined.x.N), $(grid_refined.y.N))")
    interp_counter = 0
    calc_counter = 0
    iterations_counter = 0

    cell_selector_status = fill(false, grid.x.N-1, grid.y.N-1)

    while true
        cells_to_refine = Vector{Tuple{Int64, Int64}}(undef, 0)
        for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
            if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                push!(cells_to_refine, (i_cell, j_cell))
                cell_selector_status[i_cell, j_cell] = true
            end
        end

        n_cells = length(cells_to_refine)

        if n_cells == 0
            break
        end
        iterations_counter += 1

        for cell in cells_to_refine
            i_cell, j_cell = cell
            calc_counter += calculate_cell!(i_cell, j_cell, grid, grid_refined, target_function)
            for cell in [(i_cell-1, j_cell), (i_cell, j_cell-1), (i_cell+1, j_cell), (i_cell, j_cell+1)]
                i_cell, j_cell = cell
                if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                    push!(cells_to_refine, cell)
                    cell_selector_status[i_cell, j_cell] = true
                end
            end
        end
    end

    for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
        if (cell_selector_status[i_cell,j_cell] == false)
            interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
        end
    end

    println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)
    return grid_refined
end

"""
    calculate_cell!(i_cell, j_cell, grid, grid_refined, target_function)

Calculate values for all new points introduced by refining a cell.
"""
function calculate_cell!(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, target_function)
    i_cell_ref = 2*i_cell-1
    j_cell_ref = 2*j_cell-1
    new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1,j_cell:j_cell+1]) + 1
    n_calc = 0

    for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
        if grid_refined.status[i_ref,j_ref] < 1
            is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1) 
            target_output = target_function(grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])
            n_calc += 1
            for (key, value) in pairs(target_output)
                grid_refined.vars[key][i_ref, j_ref] = value
                if !isnan(value)
                    grid_refined.min[key] = value < grid_refined.min[key] ? value : grid_refined.min[key]
                    grid_refined.max[key] = value > grid_refined.max[key] ? value : grid_refined.max[key]
                end
                if is_on_grid
                    grid.vars[key][div(i_ref,2)+1,div(j_ref,2)+1] = value
                    grid.ref_level[div(i_ref,2)+1,div(j_ref,2)+1] = new_ref_level-1
                    if !isnan(value)
                        grid.min[key] = value < grid.min[key] ? value : grid.min[key]
                        grid.max[key] = value > grid.max[key] ? value : grid.max[key]
                    end
                end
            end
            grid_refined.status[i_ref,j_ref] = 1
        end
        if grid_refined.ref_level[i_ref, j_ref] < new_ref_level
            grid_refined.ref_level[i_ref, j_ref] = new_ref_level
        end
    end

    return n_calc
end

"""
    interpolate_cell!(i, j, grid, grid_refined)

Interpolate values for new points in the refined grid that were not directly calculated.
"""
function interpolate_cell!(i::Int64, j::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid)
    i_ref = 2*i-1
    j_ref = 2*j-1
    n_inter = 0
    if grid_refined.status[i_ref, j_ref+1] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref, j_ref+1] = 0.5*(grid.vars[key][i,j]+grid.vars[key][i,j+1])
        end
        grid_refined.ref_level[i_ref, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i,j+1])
        grid_refined.status[i_ref, j_ref+1] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+1, j_ref] == -1
        for key in keys(grid_refined.vars) 
            grid_refined.vars[key][i_ref+1, j_ref] = 0.5*(grid.vars[key][i,j]+grid.vars[key][i+1,j])
        end
        grid_refined.ref_level[i_ref+1, j_ref] = min(grid.ref_level[i,j], grid.ref_level[i+1,j])
        grid_refined.status[i_ref+1, j_ref] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+1, j_ref+2] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref+1, j_ref+2] = 0.5*(grid.vars[key][i,j+1]+grid.vars[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+2] = min(grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+2] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+2, j_ref+1] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref+2, j_ref+1] = 0.5*(grid.vars[key][i+1,j]+grid.vars[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+2, j_ref+1] = min(grid.ref_level[i+1,j], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+2, j_ref+1] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+1, j_ref+1] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref+1, j_ref+1] = 0.25*(grid.vars[key][i,j]+grid.vars[key][i+1,j]+grid.vars[key][i,j+1]+grid.vars[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i+1,j], grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+1] = 0
        n_inter += 1
    end
    return n_inter
end


#--------------------------------------------------------------------------------------------------------------
# Parallel
#--------------------------------------------------------------------------------------------------------------

"""
    parallel_precalculate_2DGrid!(grid, target_function, params_function!)

Fill the grid with values using parallel processing.
"""
function parallel_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    println("\nPrecalculate ($(grid.x.N), $(grid.y.N))")
    np = nprocs()  # determine the number of available processes
    # println(Base.return_types(target_function, Tuple{Float64,Float64}))
    # target_keys = Base.return_types(target_function, Tuple{Float64,Float64})[1].parameters[1]
    for key in grid.ref_sets.params_to_save
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end

    # Counter for indices
    counter = Ref(1)
    lock_obj = ReentrantLock()

    # Function to get the next index for parallel tasks
    function nextidx()
        lock(lock_obj) do
            if counter[] > grid.x.N * grid.y.N
                return nothing  # all indices processed
            end
            idx = counter[]
            counter[] += 1
            return idx
        end
    end

    n_steps = grid.x.N * grid.y.N
    progress = Progress(n_steps, showspeed=true, dt=0.0)
    channel = RemoteChannel(()->Channel{Bool}(), 1)

    @sync begin
        @async while take!(channel)
            lock(lock_obj) do
                next!(progress)
                printstyled("\nTasks completed: $(progress.counter), Total tasks: $(progress.n)\n", color=:green)
            end
        end

        @sync for p in 1:np
            if p != myid() || np == 1
                @async while true
                    idx = nextidx()
                    if idx === nothing
                        break
                    end
                    i = div(idx - 1, grid.y.N) + 1
                    j = mod(idx - 1, grid.y.N) + 1
                    target_output = remotecall_fetch(target_function, p, grid.x.values[i], grid.y.values[j])
                    put!(channel, true)
                    for (key, value) in pairs(target_output)
                        grid.vars[key][i, j] = value
                    end
                end
            end
        end
        put!(channel, false)
    end

    for key in grid.ref_sets.params_to_save
        grid.min[key] = minimum(x -> isnan(x) ? +Inf : x, grid.vars[key])
        grid.max[key] = maximum(x -> isnan(x) ? -Inf : x, grid.vars[key])
    end
    grid.status .= 1
    println("calculations = $n_steps")
    params_function!(grid)
    return grid
end

"""
    parallel_refine_2DGrid(grid, target_function, params_function!)

Refine the grid by one level using parallel processing.
"""
function parallel_refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    grid_refined = refine(grid)
    println("\nRefinement from ($(grid.x.N), $(grid.y.N)) to ($(grid_refined.x.N), $(grid_refined.y.N))")
    
    interp_counter = 0
    calc_counter = 0
    iterations_counter = 0
    np = nprocs()

    refined_points_status = fill(false, grid_refined.x.N, grid_refined.y.N)

    lock_obj = ReentrantLock()
    condition = Condition()

    while true
        points_to_calculate = Vector{Tuple{Int64, Int64}}()
        number_of_points_to_calculate = 0

        lock(lock_obj) do
            number_of_points_to_calculate += update_points_to_calculate!(points_to_calculate, grid, grid_refined, refined_points_status)
        end

        println("Number of points to calculate: $(number_of_points_to_calculate)")
        
        if number_of_points_to_calculate == 0
            GC.gc()  # Clean up before exit if no work remains
            break
        end

        iterations_counter += 1

        tasks_issued = Ref(0)
        tasks_completed = Ref(0)
        all_tasks_completed = Ref(false)

        function next_index()
            index = nothing
            lock(lock_obj) do
                if tasks_issued[] < number_of_points_to_calculate
                    tasks_issued[] += 1
                    index = tasks_issued[]
                end
            end
            return index
        end

        progress = Progress(number_of_points_to_calculate, showspeed=true)
        
        channel = RemoteChannel(()->Channel{Bool}(), 1)

        @sync begin
            @async while take!(channel)
                lock(lock_obj) do
                    ProgressMeter.update!(progress, tasks_completed[])
                    printstyled("\nTasks issued: $(tasks_issued[]), Tasks completed: $(tasks_completed[]), Total tasks: $(number_of_points_to_calculate)\n", color=:green)
                end
            end

            @sync for p in 1:np
                if p != myid() || np == 1
                    @async begin
                        done = false
                        while !done
                            index = next_index()
                            if index === nothing
                                lock(lock_obj) do
                                    if tasks_completed[] >= tasks_issued[] && tasks_issued[] >= number_of_points_to_calculate
                                        all_tasks_completed[] = true
                                        notify(condition)
                                        done = true
                                    end
                                end

                                if !done
                                    wait(condition)
                                    lock(lock_obj) do
                                        if all_tasks_completed[]
                                            done = true
                                        end
                                    end
                                end
                            else
                                i_ref, j_ref = points_to_calculate[index]
                                # println("Before task. Used memory in GB: ", (Sys.total_memory() - Sys.free_memory()) / 1024^3)
                                calculate_point!(p, i_ref, j_ref, grid, grid_refined, target_function, lock_obj)
                                # println("After task. Used memory in GB: ", (Sys.total_memory() - Sys.free_memory()) / 1024^3)
                                lock(lock_obj) do
                                    calc_counter += 1
                                    tasks_completed[] += 1
                                    new_points_to_calculate = update_points_to_calculate!(points_to_calculate, i_ref, j_ref, grid, grid_refined, refined_points_status)
                                    number_of_points_to_calculate += new_points_to_calculate
                                    progress.n += new_points_to_calculate
                                end
                                put!(channel, true)
                                notify(condition)
                                GC.gc()  # Clean up after point calculation
                            end
                        end
                    end
                end
            end
            put!(channel, false)
        end

        GC.gc()  # Clean up after all parallel tasks
    end

    for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
        interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
    end
    GC.gc()  # Clean up after interpolation

    println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)

    GC.gc()  # Clean up before returning result
    return grid_refined
end

#--------------------------------------------------------------------------------------------------------------
# Helpers
#--------------------------------------------------------------------------------------------------------------

"""
    refined_index(i_cell)

Return the index of the refined grid corresponding to the coarse cell index.
"""
refined_index(i_cell::Int64) = (i_cell - 1) * 2 + 1

"""
    update_points_to_calculate!(points_to_calculate, grid, grid_refined, refined_points_status)

Update the list of points to calculate for refinement (parallel version).
"""
function update_points_to_calculate!(points_to_calculate::Vector{Tuple{Int64, Int64}}, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, refined_points_status)
    new_points_to_calculate = 0
    for i_cell in 1:grid.x.N - 1, j_cell in 1:grid.y.N - 1
        if cell_selector(i_cell, j_cell, grid)
            for i_ref in refined_index(i_cell):refined_index(i_cell)+2, j_ref in refined_index(j_cell):refined_index(j_cell)+2
                if refined_points_status[i_ref, j_ref] == false && grid_refined.status[i_ref, j_ref] != 1
                    push!(points_to_calculate, (i_ref, j_ref))
                    refined_points_status[i_ref, j_ref] = true
                    new_points_to_calculate += 1
                end
            end
        end
    end
    return new_points_to_calculate
end

"""
    update_points_to_calculate!(points_to_calculate, i_ref_init, j_ref_init, grid, grid_refined, refined_points_status)

Update the list of points to calculate for a given refined point (parallel version).
"""
function update_points_to_calculate!(points_to_calculate::Vector{Tuple{Int64, Int64}}, i_ref_init::Int, j_ref_init::Int, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, refined_points_status)
    new_points_to_calculate = 0
    is_on_grid = (mod(i_ref_init,2)*mod(j_ref_init,2) == 1)
    
    if is_on_grid
        i_cell_init = div(i_ref_init, 2) + 1
        j_cell_init = div(j_ref_init, 2) + 1
        for (i_cell, j_cell) in ((i_cell_init, j_cell_init), (i_cell_init - 1, j_cell_init), (i_cell_init, j_cell_init - 1), (i_cell_init - 1, j_cell_init - 1))
            if (1 <= i_cell < grid.x.N-1) && (1 <= j_cell < grid.y.N-1) && cell_selector(i_cell, j_cell, grid)
                for i_ref in refined_index(i_cell):refined_index(i_cell)+2, j_ref in refined_index(j_cell):refined_index(j_cell)+2
                    if refined_points_status[i_ref, j_ref] == false && grid_refined.status[i_ref, j_ref] != 1
                        push!(points_to_calculate, (i_ref, j_ref))
                        refined_points_status[i_ref, j_ref] = true
                        new_points_to_calculate += 1
                    end
                end
            end
        end
    end

    return new_points_to_calculate
end

"""
    calculate_point!(p, i_ref, j_ref, grid, grid_refined, target_function, lock_obj)

Remotely calculate the value at a point in the refined grid and update both grids with locking.
"""
function calculate_point!(p, i_ref, j_ref, grid, grid_refined, target_function, lock_obj)
    is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1)
    target_output = remotecall_fetch(target_function, p, grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])

    # Update the global grid and status with a temporary lock
    lock(lock_obj) do
        # println("Debug: locked in 4 on $p")
        grid_refined.status[i_ref, j_ref] = 1
        for (key, value) in pairs(target_output)
            grid_refined.vars[key][i_ref, j_ref] = value
            if !isnan(value)
                grid_refined.min[key] = value < grid_refined.min[key] ? value : grid_refined.min[key]
                grid_refined.max[key] = value > grid_refined.max[key] ? value : grid_refined.max[key]
            end
            if is_on_grid
                grid.vars[key][div(i_ref, 2) + 1, div(j_ref, 2) + 1] = value
                if !isnan(value)
                    grid.min[key] = value < grid.min[key] ? value : grid.min[key]
                    grid.max[key] = value > grid.max[key] ? value : grid.max[key]
                end
            end
        end
    end
    # println("Debug: unlocked in 4 on $p")
end

