#--------------------------------------------------------------------------------------------------------------
# Refinement units

abstract type AbstractRefinementUnit end

struct LocalMinimaUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    from_min::Bool
end

LocalMinimaUnit(name; min = -Inf, max = Inf, from_min = true) = LocalMinimaUnit(name, min, max, from_min)

function Base.show(io::IO, ru::LocalMinimaUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Local mimima unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value: ", ru.min)
    println(io, spad, "Maximum value: ", ru.max)
    print(io,   spad, "From min value: ", ru.from_min)
	return nothing
end

struct FullUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
end

FullUnit(name; min = -Inf, max = Inf) = FullUnit(name, min, max)

function Base.show(io::IO, ru::FullUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Full refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value: ", ru.min)
    print(io,   spad, "Maximum value: ", ru.max)
	return nothing
end

struct DiffUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    diff::Float64
    from_min::Bool
end

DiffUnit(name; min = -Inf, max = Inf, diff = diff, from_min = true) = DiffUnit(name, min, max, diff, from_min)

function Base.show(io::IO, ru::DiffUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Difference refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value: ", ru.min)
    println(io, spad, "Maximum value: ", ru.max)
    println(io, spad, "Maximal difference: ", ru.diff)
    print(io,   spad, "From min value: ", ru.from_min)
	return nothing
end

struct RelDiffUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    rel_diff::Float64
    from_min::Bool
end

RelDiffUnit(name; min = -Inf, max = Inf, rel_diff = rel_diff, from_min = true) = RelDiffUnit(name, min, max, rel_diff, from_min)

function Base.show(io::IO, ru::RelDiffUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Relative difference refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value: ", ru.min)
    println(io, spad, "Maximum value: ", ru.max)
    println(io, spad, "Maximal relative difference: ", ru.rel_diff)
    print(io,   spad, "From min value: ", ru.from_min)
	return nothing
end

struct ContourUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    contours::Vector{Float64}
    from_min::Bool
end

ContourUnit(name; min = -Inf, max = Inf, contours = contours, from_min = true) = ContourUnit(name, min, max, contours, from_min)

function Base.show(io::IO, ru::ContourUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Contour refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value: ", ru.min)
    println(io, spad, "Maximum value: ", ru.max)
    println(io, spad, "Contour levels: ", ru.contours)
    print(io,   spad, "From min value: ", ru.from_min)
	return nothing
end

struct DiffContourUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    diffs::Vector{Float64}
    contours::Vector{Float64}
    from_min::Bool
end

DiffContourUnit(name; min = -Inf, max = Inf, diffs = diffs, contours = contours, from_min = true) = DiffContourUnit(name, min, max, diffs, contours, from_min)

function Base.show(io::IO, ru::DiffContourUnit)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Difference and contour refinement unit:")
    println(io, spad, "Name of variable: ", ru.name)
    println(io, spad, "Minimum value: ", ru.min)
    println(io, spad, "Maximum value: ", ru.max)
    println(io, spad, "Maximal differences: ", ru.diffs)
    println(io, spad, "Contour levels: ", ru.contours)
    print(io,   spad, "From min value: ", ru.from_min)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
# Refinement settings

struct RefinementSettings{T}
    params_to_save::Tuple{Vararg{Symbol}}
    desired_refinement_level::Int64
    parallel::Bool
    units::T
end

function Base.show(io::IO, ref_sets::RefinementSettings{T}) where {T}
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    
    println(io, pad,  "Grid Refinement settings:")
    println(io, spad, "Parameters to save: ", ref_sets.params_to_save)
	println(io, spad, "Desired refinement level: ", ref_sets.desired_refinement_level)
    println(io, spad, "Parallel computation: ", ref_sets.parallel)
    for i in 1:length(ref_sets.units)-1
        println(IOContext(io, :indent => indent+4), ref_sets.units[i])
    end
    print(IOContext(io, :indent => indent+4), ref_sets.units[end])
	return nothing
end

function RefinementSettings(units...; desired_refinement_level::Int64, parallel::Bool, params_to_save::Tuple{Vararg{Symbol}})
    return RefinementSettings(params_to_save, desired_refinement_level, parallel, units)
end

#--------------------------------------------------------------------------------------------------------------
# Refinement 2D grid

abstract type General2DGrid end

struct AdaptiveRefinement2DGrid{T1 <: AbstractRangeRule, T2 <: AbstractRangeRule, T3} <: General2DGrid
    vars::Dict{Symbol,Matrix{Float64}}
    params::Dict{Symbol,Float64}
    min::Dict{Symbol,Float64}
    max::Dict{Symbol,Float64}
    x::RangeVariable{T1}
    y::RangeVariable{T2}
    ref_sets::RefinementSettings{T3}
    ref_level::Matrix{Int64}
    status::Matrix{Int64}
end

function AdaptiveRefinement2DGrid(x::RangeVariable{T1}, y::RangeVariable{T2}, ref_sets::T3) where {T1 <: AbstractRangeRule, T2 <: AbstractRangeRule, T3}
    vars = Dict{Symbol,Matrix{Float64}}()
    params = Dict{Symbol,Float64}()
    min = Dict{Symbol,Float64}()
    max = Dict{Symbol,Float64}()
    ref_level = [0 for i in 1:x.N, j in 1:y.N]
    status = [-1 for i in 1:x.N, j in 1:y.N]
    return AdaptiveRefinement2DGrid(vars, params, min, max, x, y, ref_sets, ref_level, status)
end

function Base.show(io::IO, grid::AdaptiveRefinement2DGrid)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    
    println(io, pad,  "AdaptiveRefinement2DGrid:")
	println(io, spad, "Variables: ", keys(grid.vars))
    println(io, spad, "Parameters: ", grid.params)
    println(io, spad, "Minimal values: ", grid.min)
    println(io, spad, "Maximal values: ", grid.max)
    println(io, spad, "X axis: ", grid.x)
    println(io, spad, "Y axis: ", grid.y)
    print(IOContext(io, :indent => indent+2), grid.ref_sets)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
# Cell selectors for refinement units

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

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::LocalMinimaUnit; at_corner=false)
    cell     = @view grid.vars[ref_unit.name][i_cell : min(i_cell+1, grid.x.N) , j_cell : min(j_cell+1, grid.y.N)]
    big_cell = @view grid.vars[ref_unit.name][max(i_cell-1, 1) : min(i_cell+2, grid.x.N), max(j_cell-1, 1) : min(j_cell+2, grid.y.N)]
    value_big_cell_min = minimum(big_cell)
    value_cell_min = minimum(cell)
    value_min = grid.min[ref_unit.name]

    # Проверка, что минимум в big_cell находится в пределах cell
    local_minima_case = value_big_cell_min in cell

    # Проверка уникальности минимума внутри big_cell
    unique_minimum_case = sum(big_cell .== value_big_cell_min) - sum(cell .== value_big_cell_min) == 0

    at_corner_case = grid.vars[:chisqr][i_cell,j_cell] == value_big_cell_min

    if ref_unit.from_min
        min_case = ref_unit.min <= value_big_cell_min - value_min
        max_case = ref_unit.max >= value_big_cell_min - value_min
    else
        min_case = ref_unit.min <= value_big_cell_min
        max_case = ref_unit.max >= value_big_cell_min
    end

    # if local_minima_case && unique_minimum_case
    #     display(cell)
    #     display(big_cell)
    #     println(value_big_cell_min) 
    # end

    return local_minima_case && unique_minimum_case && min_case && max_case && (!at_corner || at_corner_case)
end

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::FullUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    return min_case && max_case
end

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

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::DiffContourUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    value_min = grid.min[ref_unit.name]
    if ref_unit.from_min 
        diff_contour_case = any((value_cell_min .< value_min .+ ref_unit.contours .< value_cell_max) .&& (value_cell_max .- value_cell_min .> ref_unit.diffs))
        # diff_case    = any(value_cell_max .- value_cell_min > ref_unit.diff)
        min_case = ref_unit.min <= value_cell_min - value_min
        max_case = ref_unit.max >= value_cell_max - value_min
    else
        contour_case = any((value_cell_min .< ref_unit.contours .< value_cell_max) .&& (value_cell_max .- value_cell_min .> ref_unit.diffs))
        min_case = ref_unit.min <= value_cell_min
        max_case = ref_unit.max >= value_cell_max
    end
    return diff_contour_case && min_case && max_case
end

#--------------------------------------------------------------------------------------------------------------
# Supplementary refinement routines

function refine(var::RangeVariable)
    var_refined = typeof(var)(var.name, var.min, var.max, var.N*2-1, var.range_rule)
    return var_refined
end

#=
function refine_1Darray(x::Vector{Float64})
    x_refined = Vector{Float64}(undef, length(x)*2-1)
    for i in 1:length(x)-1
        x_refined[2*i-1] = x[i]
        x_refined[2*i] = 0.5*(x[i+1] + x[i])
    end
    x_refined[end] = x[end]
    return x_refined
end
=#

function refine(arr::Matrix{T}) where {T}
    arr_refined = fill(-one(T), 2 .* size(arr) .- 1)
    for i in 1:size(arr)[1], j in 1:size(arr)[2]
        arr_refined[2*i-1,2*j-1] = arr[i,j]
    end
    return arr_refined::Matrix{T}
end

function refine(dict::Dict{Symbol,Matrix{T}}) where {T}
    dict_refined = Dict{Symbol,Matrix{T}}()
    for (key, value) in dict
        dict_refined[key] = refine(dict[key])
    end
    return dict_refined::Dict{Symbol,Matrix{T}}
end

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
# General routines

function calculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    precalculate_2DGrid!(grid, target_function, params_function!)
    for i in 1:grid.ref_sets.desired_refinement_level
        grid = refine_2DGrid(grid, target_function, params_function!)
    end
    return grid
end

function precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    if grid.ref_sets.parallel == true
        return parallel_precalculate_2DGrid!(grid, target_function, params_function!)
    else 
        return single_core_precalculate_2DGrid!(grid, target_function, params_function!)
    end
end

function refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    if grid.ref_sets.parallel == true
        return parallel_refine_2DGrid(grid, target_function, params_function!)
    else 
        return single_core_refine_2DGrid(grid, target_function, params_function!)
    end
end

#--------------------------------------------------------------------------------------------------------------
# Singe-core routines

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
        n_inter
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
# Parallel routines

function parallel_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    println("\nPrecalculate ($(grid.x.N), $(grid.y.N))")
    np = nprocs()  # определяем количество доступных процессов
    # println(Base.return_types(target_function, Tuple{Float64,Float64}))
    # target_keys = Base.return_types(target_function, Tuple{Float64,Float64})[1].parameters[1]
    for key in grid.ref_sets.params_to_save
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end

    # Счетчик для индексов
    counter = Ref(1)
    lock_obj = ReentrantLock()

    # Функция для получения следующего индекса
    function nextidx()
        lock(lock_obj) do
            if counter[] > grid.x.N * grid.y.N
                return nothing  # все индексы обработаны
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

#=
function pmap_parallel_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)

    target_keys, target_values = target_function(grid.x.values[1], grid.y.values[1], only_keys = true)
    for (i_key, key) in enumerate(target_keys)
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end

    function separate_task(idx::Tuple{Int64,Int64})
        i, j = idx
        calc_keys, calc_values = target_function(grid.x.values[i], grid.y.values[j])
        for (i_key, key) in enumerate(calc_keys)
            grid.vars[key][i, j] = calc_values[i_key]
        end
    end

    ans = pmap(args -> target_function(args...), [(x,y) for x in grid.x, y in grid.y])

    for i in 1:grid.x.N, j in 1:grid.y.N
        calc_keys, calc_values = ans[i,j]
        for (i_key, key) in enumerate(calc_keys)
            grid.vars[key][i, j] = calc_values[i_key]
        end
    end

    grid.status .= 1
    params_function!(grid)
    return grid
end
=#

function parallel_refine_2DGrid_old(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    grid_refined = refine(grid)
    println("\nRefinement from ($(grid.x.N), $(grid.y.N)) to ($(grid_refined.x.N), $(grid_refined.y.N))")
    
    interp_counter = 0
    calc_counter = 0
    iterations_counter = 0
    np = nprocs()

    update_calc_counter(n_calc) = (calc_counter += n_calc)

    cell_selector_status = fill(false, grid.x.N-1, grid.y.N-1)
    lock_obj = ReentrantLock()  # для синхронизации
    condition = Condition()  # Условие для управления задачами

    while true
        iterations_counter += 1
        cells_to_refine = Vector{Tuple{Int64, Int64}}(undef, 0)
        for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
            if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                push!(cells_to_refine, (i_cell, j_cell))
                cell_selector_status[i_cell, j_cell] = true
            end
        end

        n_cells = length(cells_to_refine)
        println("cells to refine: $n_cells")
        
        if n_cells == 0
            break
        end

        tasks_issued = Ref(0)  # Счетчик выданных задач
        tasks_completed = Ref(0)  # Счетчик завершенных задач
        all_done = Ref(false)  # Указатель на завершение всех задач

        function nextidx()
            idx = nothing
            lock(lock_obj) do
                if tasks_issued[] < n_cells
                    tasks_issued[] += 1
                    idx = tasks_issued[]
                end
            end
            return idx
        end

        progress = Progress(n_cells, showspeed=true)  # Инициализация прогрессметра
        
        channel = RemoteChannel(()->Channel{Bool}(), 1)

        @sync begin
            @async while take!(channel)
                lock(lock_obj) do
                    ProgressMeter.update!(progress, tasks_completed[])  # Обновление прогресса завершенных задач
                    printstyled("\nTasks issued: $(tasks_issued[]), Tasks completed: $(tasks_completed[]), Total tasks: $n_cells\n", color=:green)
                end
            end

            @sync for p in 1:np
                if p != myid() || np == 1
                    @async begin
                        done = false  # Флаг для выхода из цикла
                        while !done
                            idx = nextidx()
                            println("Debug: index $idx on $p")
                            if idx === nothing
                                # Проверяем, завершены ли все задачи
                                lock(lock_obj) do
                                    if tasks_completed[] >= tasks_issued[] && tasks_issued[] >= n_cells
                                        all_done[] = true
                                        notify(condition)  # Уведомляем всех, что задачи завершены
                                        done = true  # Устанавливаем флаг для выхода из цикла
                                    end
                                end
                                
                                if !done
                                    println("Debug: waiting on $p")
                                    wait(condition)  # Ожидаем уведомления о новых задачах
                                    lock(lock_obj) do
                                        if all_done[]  # Если все задачи завершены, выходим
                                            done = true
                                        end
                                    end
                                end
                            else
                                i_cell, j_cell = cells_to_refine[idx]
                                println("Debug: started cell ($i_cell, $j_cell) on $p")
                                n_calc = calculate_cell!(p, i_cell, j_cell, grid, grid_refined, target_function)
                                println("Debug: finished cell ($i_cell, $j_cell) on $p")
                                update_calc_counter(n_calc)

                                # Обработка соседних ячеек
                                lock(lock_obj) do
                                    for cell in [(i_cell-1, j_cell), (i_cell, j_cell-1), (i_cell+1, j_cell), (i_cell, j_cell+1)]
                                        i_cell_new, j_cell_new = cell
                                        if i_cell_new >= 1 && i_cell_new < grid.x.N && j_cell_new >= 1 && j_cell_new < grid.y.N
                                            if cell_selector(i_cell_new, j_cell_new, grid) && !cell_selector_status[i_cell_new, j_cell_new]
                                                println("Debug: added new cell ($i_cell_new, $j_cell_new) on $p")
                                                push!(cells_to_refine, cell)
                                                cell_selector_status[i_cell_new, j_cell_new] = true
                                                n_cells += 1
                                                progress.n += 1
                                                notify(condition)  # Уведомляем о новой задаче
                                            end
                                        end
                                    end
                                    tasks_completed[] += 1  # Увеличиваем счетчик завершенных задач
                                end

                                put!(channel, true)
                            end
                        end
                    end
                end
            end
            put!(channel, false)
        end
    end

    # Интерполяция оставшихся ячеек
    for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
        if !cell_selector_status[i_cell, j_cell]
            interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
        end
    end

    println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)
    return grid_refined
end

function calculate_cell!_old(p, i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, target_function)
    i_cell_ref = 2*i_cell - 1
    j_cell_ref = 2*j_cell - 1
    new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1, j_cell:j_cell+1]) + 1
    n_calc = 0
    lock_obj = ReentrantLock()  # для синхронизации

    for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
        # Выполнение основной работы без блокировки
        if grid_refined.status[i_ref, j_ref] < 1
            is_on_grid = (mod(i_ref, 2) * mod(j_ref, 2) == 1)
            target_output = remotecall_fetch(target_function, p, grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])
            n_calc += 1

            # Обновление общей сетки и статуса с кратковременной блокировкой
            lock(lock_obj) do
                grid_refined.status[i_ref, j_ref] = 1
                for (key, value) in pairs(target_output)
                    grid_refined.vars[key][i_ref, j_ref] = value
                    if !isnan(value)
                        grid_refined.min[key] = value < grid_refined.min[key] ? value : grid_refined.min[key]
                        grid_refined.max[key] = value > grid_refined.max[key] ? value : grid_refined.max[key]
                    end

                    if is_on_grid
                        grid.vars[key][div(i_ref, 2) + 1, div(j_ref, 2) + 1] = value
                        grid.ref_level[div(i_ref, 2) + 1, div(j_ref, 2) + 1] = new_ref_level - 1
                        if !isnan(value)
                            grid.min[key] = value < grid.min[key] ? value : grid.min[key]
                            grid.max[key] = value > grid.max[key] ? value : grid.max[key]
                        end
                    end
                end
            end
        end

        # Обновление уровня уточнения
        lock(lock_obj) do
            if grid_refined.ref_level[i_ref, j_ref] < new_ref_level
                grid_refined.ref_level[i_ref, j_ref] = new_ref_level
            end
        end
    end

    return n_calc
end

# function parallel_refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
#     grid_refined = refine(grid)
#     println("\nRefinement from ($(grid.x.N), $(grid.y.N)) to ($(grid_refined.x.N), $(grid_refined.y.N))")
#     interp_counter = 0
#     calc_counter = 0
#     iterations_counter = 0
#     np = nprocs()

#     update_calc_counter(n_calc) = (calc_counter+=n_calc)

#     cell_selector_status = fill(false, grid.x.N-1, grid.y.N-1)

#     while true
#         cells_to_refine = Vector{Tuple{Int64, Int64}}(undef, 0)
#         for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
#             if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
#                 push!(cells_to_refine, (i_cell, j_cell))
#                 cell_selector_status[i_cell, j_cell] = true
#             end
#         end

#         n_cells = length(cells_to_refine)

#         println("cells to refine: $n_cells")

#         if n_cells == 0
#             break
#         end
#         iterations_counter += 1

#         counter_cells = 1
#         nextidx() = (idx=counter_cells; counter_cells+=1; idx)
        
#         progress = Progress(n_cells, showspeed=true)
#         channel = RemoteChannel(()->Channel{Bool}(), 1)

#         add_new_cell() = (n_cells += 1; progress.n += 1)

#         @sync begin
#             @async while take!(channel)
#                 next!(progress)
#                 println("")
#                 # println("$(now()): Ready $(p.counter) from $(p.n)")
#             end
#             @sync begin
#                 for p=1:np
#                     if p != myid() || np == 1
#                         @async begin 
#                             while true
#                                 idx = nextidx()
#                                 if idx > n_cells
#                                     break
#                                 end
#                                 i_cell, j_cell = cells_to_refine[idx]
#                                 n_calc = calculate_cell!(p, i_cell, j_cell, grid, grid_refined, target_function)
#                                 update_calc_counter(n_calc)
#                                 # for cell in [(i_cell-1, j_cell), (i_cell, j_cell-1), (i_cell+1, j_cell), (i_cell, j_cell+1)] # maybe error on edges
#                                 #     i_cell, j_cell = cell
#                                 #     if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
#                                 #         push!(cells_to_refine, cell)
#                                 #         cell_selector_status[i_cell, j_cell] = true
#                                 #         add_new_cell()
#                                 #     end
#                                 # end
#                                 put!(channel, true)
#                             end
#                         end
#                     end
#                 end
#             end
#             put!(channel, false)
#         end
#     end

#     for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
#         if (cell_selector_status[i_cell,j_cell] == false)
#             interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
#         end
#     end

#     println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
#     params_function!(grid)
#     return grid_refined
# end

# function calculate_cell!(p, i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, target_function)
#     i_cell_ref = 2*i_cell-1
#     j_cell_ref = 2*j_cell-1
#     new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1,j_cell:j_cell+1]) + 1
#     n_calc = 0
#     lock = ReentrantLock()  # для синхронизации

#     for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
#         lock(lock) do
#             if grid_refined.status[i_ref,j_ref] < 1
#                 is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1)
#                 grid_refined.status[i_ref,j_ref] = 1
#                 target_output = remotecall_fetch(target_function, p, grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])
#                 n_calc += 1
#                 for (key, value) in pairs(target_output)
#                     grid_refined.vars[key][i_ref, j_ref] = value
#                     if !isnan(value)
#                         grid_refined.min[key] = value < grid_refined.min[key] ? value : grid_refined.min[key]
#                         grid_refined.max[key] = value > grid_refined.max[key] ? value : grid_refined.max[key]
#                     end
#                     if is_on_grid
#                         grid.vars[key][div(i_ref,2)+1,div(j_ref,2)+1] = value
#                         grid.ref_level[div(i_ref,2)+1,div(j_ref,2)+1] = new_ref_level-1
#                         if !isnan(value)
#                             grid.min[key] = value < grid.min[key] ? value : grid.min[key]
#                             grid.max[key] = value > grid.max[key] ? value : grid.max[key]
#                         end
#                     end
#                 end
#             end
#             if grid_refined.ref_level[i_ref, j_ref] < new_ref_level
#                 grid_refined.ref_level[i_ref, j_ref] = new_ref_level
#             end
#         end
#     end

#     return n_calc
# end

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
            GC.gc()  # ✅ Очистка перед выходом, если работы нет
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
                                GC.gc()  # ✅ Очистка после вычисления точки
                            end
                        end
                    end
                end
            end
            put!(channel, false)
        end

        GC.gc()  # ✅ Очистка после завершения всех параллельных задач
    end

    for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
        interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
    end
    GC.gc()  # ✅ Очистка после интерполяции

    println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)

    GC.gc()  # ✅ Очистка перед возвратом результата
    return grid_refined
end

refined_index(i_cell::Int64) = (i_cell - 1) * 2 + 1

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

function calculate_point!(p, i_ref, j_ref, grid, grid_refined, target_function, lock_obj)
    is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1)
    target_output = remotecall_fetch(target_function, p, grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])

    # Обновление общей сетки и статуса с кратковременной блокировкой
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


