using Test
using GravityToolsNext
using JLD2: load
using Distributions: Normal

import GravityToolsNext: run_task, task_workdir, task_copy_with, task_derive_par_output, task_stage_inputs!

Base.@kwdef struct SyntheticWrapperTask <: SingleTempoTask
    work_dir::String
    par_output::String = "synthetic_out.par"
    job_name::Union{Nothing,String} = nothing
    work_mode::Symbol = :inplace
    overwrite::Symbol = :reuse
    temp_dir::Union{Nothing,String} = nothing
    par_input::Union{Nothing,String} = nothing
    snapshot_par::Bool = true
    io_mirror::Any = :none
    override_params::Vector{TempoParameter} = TempoParameter[]
    stage_source::Union{Nothing,String} = nothing
end

struct MissingHookTask <: SingleTempoTask
    work_dir::String
end

task_workdir(task::SyntheticWrapperTask) = task.work_dir
task_workdir(task::MissingHookTask) = task.work_dir

function _override_dict(params::Vector{TempoParameter})
    out = Dict{Symbol,Float64}()
    for p in params
        if p.value !== nothing
            out[p.name_symbol] = Float64(value_as_float(p))
        end
    end
    return out
end

function _synthetic_iteration(; pre_post::Float64=1.0, fit_parameters::Vector{FitParameter}=FitParameter[])
    basic = BasicTempoOutput(
        10.0,
        5.0,
        5.0,
        25.0,
        5,
        5.0,
        pre_post,
        length(fit_parameters),
        16,
        0.0,
        0.0,
        0.0,
    )
    return InternalIterationResult(
        InternalIterationOutput(basic, fit_parameters, TempoOutputError()),
        nothing,
        nothing,
        nothing,
        Dict{Symbol,Any}(),
    )
end

function _synthetic_result(;
    metrics::Dict{Symbol,Float64},
    job_root::AbstractString,
    run_cwd::AbstractString,
    par_out_path::AbstractString,
    success::Bool=true,
    status::Symbol=:ok,
    fit_parameters::Vector{FitParameter}=FitParameter[],
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)
    iter = _synthetic_iteration(
        pre_post = get(metrics, :pre_post_final, 1.0),
        fit_parameters = fit_parameters,
    )
    meta = Dict{Symbol,Any}(
        :job_root => String(job_root),
        :run_cwd => String(run_cwd),
        :par_out_path => String(par_out_path),
        :success => success,
        :status => status,
    )
    merge!(meta, metadata)
    return build_general_tempo_result(
        [iter];
        metadata = meta,
        metrics_hook = (_, _) -> metrics,
    )
end

function _par_output_for(task::SyntheticWrapperTask, tag::AbstractString)
    stem, _ = splitext(basename(task.par_output))
    if endswith(stem, "_out")
        stem = first(stem, lastindex(stem) - 4)
    end
    return string(stem, "_", tag, "_out.par")
end

function task_copy_with(task::SyntheticWrapperTask; kwargs...)
    base_overrides = get(kwargs, :override_params, task.override_params)
    overrides = copy(base_overrides)
    upserts = get(kwargs, :override_params_upsert, TempoParameter[])
    if !isempty(upserts)
        overrides = with_upserted_params(overrides, upserts)
    end

    return SyntheticWrapperTask(
        work_dir = String(get(kwargs, :work_dir, task.work_dir)),
        par_output = String(get(kwargs, :par_output, task.par_output)),
        job_name = get(kwargs, :job_name, task.job_name),
        work_mode = get(kwargs, :work_mode, task.work_mode),
        overwrite = get(kwargs, :overwrite, task.overwrite),
        temp_dir = get(kwargs, :temp_dir, task.temp_dir),
        par_input = get(kwargs, :par_input, task.par_input),
        snapshot_par = get(kwargs, :snapshot_par, task.snapshot_par),
        io_mirror = get(kwargs, :io_mirror, task.io_mirror),
        override_params = overrides,
        stage_source = get(kwargs, :stage_source, task.stage_source),
    )
end

task_derive_par_output(task::SyntheticWrapperTask, tag::AbstractString) = _par_output_for(task, tag)

function task_stage_inputs!(task::SyntheticWrapperTask, dest_dir::AbstractString)
    mkpath(dest_dir)
    task.stage_source === nothing && return nothing
    cp(task.stage_source, joinpath(dest_dir, basename(task.stage_source)); force=true)
    return nothing
end

function run_task(task::SyntheticWrapperTask)::GeneralTempoResult
    overrides = _override_dict(task.override_params)
    x = get(overrides, :X, 0.0)
    y = get(overrides, :Y, 0.0)
    theta = get(overrides, :THETA, 0.0)

    run_cwd = if task.temp_dir !== nothing
        joinpath(task.work_dir, task.temp_dir)
    elseif task.work_mode === :jobdir && task.job_name !== nothing
        joinpath(task.work_dir, task.job_name)
    else
        task.work_dir
    end
    mkpath(run_cwd)

    par_out_path = joinpath(run_cwd, task.par_output)
    write(par_out_path, "synthetic par output\n")

    fit_parameters = FitParameter[
        FitParameter("x_post", x, x + 0.25, 0.1, 0.25, true),
        FitParameter("theta_post", theta, theta + 0.5, 0.1, 0.5, true),
    ]

    metrics = Dict{Symbol,Float64}(
        :chi2_fit => theta^2 + x^2 + y^2 + 1.0,
        :grid_metric => 10.0 * x + y,
        :theta_metric => theta,
        :pre_post_final => 1.0,
    )

    metadata = Dict{Symbol,Any}(
        :task_job_name => task.job_name,
        :task_temp_dir => task.temp_dir,
        :task_par_input => task.par_input,
        :task_snapshot_par => task.snapshot_par,
        :task_overwrite => task.overwrite,
        :task_work_mode => task.work_mode,
        :task_io_mirror => task.io_mirror,
        :task_override_names => [p.name_symbol for p in task.override_params],
    )

    return _synthetic_result(
        metrics = metrics,
        job_root = task.work_dir,
        run_cwd = run_cwd,
        par_out_path = par_out_path,
        fit_parameters = fit_parameters,
        metadata = metadata,
    )
end

function run_task(task::MissingHookTask)::GeneralTempoResult
    workdir = task.work_dir
    mkpath(workdir)
    return _synthetic_result(
        metrics = Dict{Symbol,Float64}(
            :chi2_fit => 1.0,
            :pre_post_final => 1.0,
        ),
        job_root = workdir,
        run_cwd = workdir,
        par_out_path = joinpath(workdir, "missing_hook_out.par"),
    )
end

function _param_value(params::Vector{TempoParameter}, name::Symbol)
    idx = findfirst(p -> p.name_symbol == name, params)
    idx === nothing && return nothing
    return Float64(value_as_float(params[idx]))
end

@testset "Wrapper surface hooks" begin
    fixture_root = normpath(joinpath(@__DIR__, "data", "ddstg_smoke"))

    @testset "BasicTempoTask task_copy_with and task_derive_par_output" begin
        settings = TempoRunSettings(
            work_dir = fixture_root,
            par_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.par",
            tim_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.tim",
            par_output = "j1141-6545_ddstg_gr_pdfb1_white_noise_out.par",
            tempo_version = Tempo2(fixture_root),
            work_mode = :jobdir,
            job_name = "wrapper_surface",
            overwrite = :reuse,
            layout = :split,
            override_params = [TP("F0", 1.0; flag=1)],
        )

        task = BasicTempoTask(settings)
        copied = task_copy_with(
            task;
            work_dir = joinpath(fixture_root, "wrapper_surface_copy"),
            par_output = "custom.par",
            job_name = "wrapper_surface_copy",
            override_params_upsert = [TP("F0", 2.0; flag=0), TP("F1", 3.0; flag=1)],
        )

        @test copied isa BasicTempoTask
        @test copied.settings.paths.work_dir == joinpath(fixture_root, "wrapper_surface_copy")
        @test copied.settings.paths.par_output == "custom.par"
        @test copied.settings.workspace.job_name == "wrapper_surface_copy"
        @test _param_value(task.settings.modifiers.override_params, :F0) == 1.0
        @test _param_value(copied.settings.modifiers.override_params, :F0) == 2.0
        @test _param_value(copied.settings.modifiers.override_params, :F1) == 3.0
        @test task_derive_par_output(task, "node001") == "j1141-6545_ddstg_gr_pdfb1_white_noise_node001_out.par"
        @test task_derive_par_output(copied, "gridA") == "custom_gridA_out.par"
    end

    @testset "Missing hook failures stay crisp" begin
        missing = MissingHookTask(mktempdir())

        err = try
            task_copy_with(missing; job_name = "copy")
            nothing
        catch exc
            exc
        end
        @test err isa ErrorException
        @test occursin("task_copy_with not implemented", sprint(showerror, err))

        err = try
            task_derive_par_output(missing, "node001")
            nothing
        catch exc
            exc
        end
        @test err isa ErrorException
        @test occursin("task_derive_par_output not implemented", sprint(showerror, err))

        prior = AnalyticPrior(Normal(0.0, 1.0))
        settings = PriorMarginalizationSettings(
            parameter = :THETA,
            pin_mode = :fixed,
            prior = prior,
            nodes = ExplicitThetaNodes([0.0]),
        )
        wrapper = PriorMarginalizedTempoTask(missing, settings)

        err = try
            run_task(wrapper)
            nothing
        catch exc
            exc
        end
        @test err isa ErrorException
        @test occursin("task_derive_par_output not implemented", sprint(showerror, err))

        err = try
            save_result_jld2(42; filename = "ignored.jld2")
            nothing
        catch exc
            exc
        end
        @test err isa ErrorException
        @test occursin("save_result_jld2 not implemented", sprint(showerror, err))
    end

    @testset "save_result_jld2 persists under job_root/results" begin
        job_root = mktempdir()
        result = _synthetic_result(
            metrics = Dict{Symbol,Float64}(
                :chi2_fit => 12.5,
                :custom_metric => 3.0,
                :pre_post_final => 1.0,
            ),
            job_root = job_root,
            run_cwd = job_root,
            par_out_path = joinpath(job_root, "synthetic_out.par"),
            fit_parameters = [FitParameter("F0", 1.0, 2.0, 0.1, 1.0, true)],
        )

        save_result_jld2(result; filename = "nested/final_result.jld2")
        saved_path = joinpath(job_root, "results", "nested", "final_result.jld2")

        @test isfile(saved_path)

        loaded = load(saved_path, "res")
        @test loaded isa GeneralTempoResult
        @test loaded.metrics[:custom_metric] == 3.0
        @test loaded.metadata[:job_root] == job_root
        @test loaded.param_estimates[:F0].value == 2.0
    end

    @testset "PriorMarginalizedTempoTask serial wrapper assumptions" begin
        base_dir = mktempdir()
        base_task = SyntheticWrapperTask(
            work_dir = base_dir,
            par_output = "prior_base_out.par",
        )

        settings = PriorMarginalizationSettings(
            parameter = :THETA,
            pin_mode = :fixed,
            prior = AnalyticPrior(Normal(0.0, 0.2)),
            nodes = ExplicitThetaNodes([-1.0, 0.0, 1.0]),
            exec_options = PriorExecutionOptions(
                mode = :chained,
                scheduler = :serial,
                chain_snapshot_par = true,
                node_dir_prefix = "nodes/node_",
                dir_name_mode = :index_only,
            ),
        )

        wrapper = PriorMarginalizedTempoTask(base_task, settings)
        result = run_task(wrapper)

        @test result isa GeneralTempoResult
        @test result.success
        @test result.subresult_type == :prior_nodes
        @test length(result.subresults) == 3
        @test result.metadata[:prior_parameter] == :THETA
        @test result.metrics[:ref_theta] == 0.0
        @test result.metrics[:post_mode_grid] == 0.0
        @test length(result.metadata[:node_workdirs]) == 3

        node1, node2, node3 = result.subresults

        @test node1.metadata[:task_temp_dir] == "nodes/node_001"
        @test node2.metadata[:task_temp_dir] == "nodes/node_002"
        @test node3.metadata[:task_temp_dir] == "nodes/node_003"
        @test node1.metadata[:task_overwrite] == :clean
        @test node2.metadata[:task_overwrite] == :reuse
        @test node3.metadata[:task_overwrite] == :reuse
        @test node1.metadata[:task_par_input] === nothing
        @test node2.metadata[:task_par_input] == relpath(node1.metadata[:par_out_path], base_dir)
        @test node3.metadata[:task_par_input] == relpath(node2.metadata[:par_out_path], base_dir)
        @test node2.metadata[:task_snapshot_par] == true
        @test node3.metadata[:task_snapshot_par] == true
        @test node1.metadata[:task_io_mirror] == (:depth_minus, 2)
        @test node2.metadata[:task_io_mirror] == (:depth_minus, 2)
        @test node3.metadata[:task_io_mirror] == (:depth_minus, 2)
    end

    @testset "Adaptive2DGridTask contract and result shape" begin
        base_dir = mktempdir()
        staged_source = joinpath(base_dir, "synthetic_input.tim")
        write(staged_source, "synthetic tim input\n")

        base_task = SyntheticWrapperTask(
            work_dir = base_dir,
            par_output = "grid_base_out.par",
            stage_source = staged_source,
        )

        x = GridAxis(:X, [1.0, 2.0])
        y = GridAxis(:Y, [10.0, 20.0])
        ref_settings = RefinementSettings(
            FullUnit(:grid_metric);
            desired_refinement_level = 0,
            parallel = false,
            params_to_save = (:grid_metric, :x_post, :missing_metric),
        )
        opts = GridWorkspaceOptions(
            grid_root = "grid_workspace",
            point_job_prefix = "points",
            results_dirname = "saved_results",
            input_dirname = "inputs",
            save_results_jld2 = true,
            stage_inputs = true,
            stage_inputs_mode = :subdir,
        )

        task = Adaptive2DGridTask(
            base_task = base_task,
            x = x,
            y = y,
            ref_settings = ref_settings,
            opts = opts,
        )

        grid = redirect_stdout(devnull) do
            run_task(task)
        end

        @test grid isa AdaptiveRefinement2DGrid
        @test haskey(grid.vars, :grid_metric)
        @test haskey(grid.vars, :x_post)
        @test haskey(grid.vars, :missing_metric)
        @test size(grid.vars[:grid_metric]) == (2, 2)
        @test grid.vars[:grid_metric][1, 1] == 20.0
        @test grid.vars[:grid_metric][2, 2] == 40.0
        @test grid.vars[:x_post][1, 1] == 1.25
        @test grid.vars[:x_post][2, 1] == 2.25
        @test all(isnan, grid.vars[:missing_metric])

        staged_input = joinpath(base_dir, "grid_workspace", "inputs", basename(staged_source))
        @test isfile(staged_input)

        results_dir = joinpath(base_dir, "grid_workspace", "saved_results")
        saved_files = filter(path -> endswith(path, ".jld2"), readdir(results_dir; join=true))
        @test length(saved_files) == 4

        sample = load(first(saved_files), "res")
        @test sample isa GeneralTempoResult
        @test haskey(sample.metrics, :grid_metric)
    end
end
