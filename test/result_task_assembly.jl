using Test
using GravityToolsNext

const GTN = GravityToolsNext

function _result_fixture_paths()
    fixture_root = normpath(joinpath(@__DIR__, "data", "ddstg_smoke"))
    return (
        root = fixture_root,
        par = joinpath(fixture_root, "j1141-6545_ddstg_gr_pdfb1_white_noise.par"),
        tim = joinpath(fixture_root, "j1141-6545_ddstg_gr_pdfb1_white_noise.tim"),
        smoke_out = joinpath(fixture_root, "DDSTG_SMOKE", "output", "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke.out"),
        smoke_par_out = joinpath(fixture_root, "DDSTG_SMOKE", "output", "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par"),
        smoke_residual = joinpath(fixture_root, "DDSTG_SMOKE", "tmp", "residuals_1.dat"),
    )
end

function _copy_result_inputs!(dest_dir::AbstractString)
    f = _result_fixture_paths()
    mkpath(dest_dir)
    par_name = basename(f.par)
    tim_name = basename(f.tim)
    cp(f.par, joinpath(dest_dir, par_name); force=true)
    cp(f.tim, joinpath(dest_dir, tim_name); force=true)
    return (; par_name, tim_name)
end

function _basic_task_settings(work_dir::AbstractString)
    files = _copy_result_inputs!(work_dir)
    return TempoRunSettings(
        work_dir = String(work_dir),
        par_input = files.par_name,
        tim_input = files.tim_name,
        par_output = "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par",
        tempo_version = Tempo2(_result_fixture_paths().root),
        work_mode = :jobdir,
        job_name = "basic-task-job",
        overwrite = :clean,
        layout = :split,
        white_noise_enabled = false,
        save_residuals = true,
    )
end

function _write_fake_basic_task_script(path::AbstractString, fixture_paths)
    script = replace(raw"""#!/bin/sh
set -eu
fixture_out="__FIXTURE_OUT__"
fixture_par_out="__FIXTURE_PAR_OUT__"
fixture_residual="__FIXTURE_RESIDUAL__"
outpar=""
while [ "$#" -gt 0 ]; do
  case "$1" in
    -outpar)
      shift
      outpar="$1"
      ;;
  esac
  shift
done

if [ -n "$outpar" ]; then
  cp "$fixture_par_out" "$outpar"
fi
cp "$fixture_residual" "residuals.dat"
cat "$fixture_out"
""",
        "__FIXTURE_OUT__" => String(fixture_paths.smoke_out),
        "__FIXTURE_PAR_OUT__" => String(fixture_paths.smoke_par_out),
        "__FIXTURE_RESIDUAL__" => String(fixture_paths.smoke_residual),
    )
    write(path, script)
    chmod(path, 0o755)
    return path
end

function _synthetic_iteration_block(;
    header::AbstractString = "Complete fit",
    rms_pre::Float64 = 10.0,
    rms_post::Float64 = 5.0,
    rms_tn::Union{Nothing,Float64} = 4.5,
    chisq::Float64 = 100.0,
    nfree::Int = 10,
    chisq_red::Float64 = 10.0,
    pre_post::Float64 = 2.0,
    nfit::Int = 2,
    npts::Int = 12,
    offset::Float64 = 1.0,
    offset_error::Float64 = 0.1,
    offset_esn::Float64 = 0.3,
    f0_post::Float64 = 1.1,
    f1_post::Float64 = 2.1,
)
    tn_line = rms_tn === nothing ? "" : "RMS post-fit residual TN = $(rms_tn) (us)\n"
    return """
$header
 RMS pre-fit residual = $(rms_pre) (us), RMS post-fit residual = $(rms_post) (us)
 $(tn_line)Fit Chisq = $(chisq) Chisqr/nfree = $(chisq)/$(nfree) = $(chisq_red) pre/post = $(pre_post)
Number of fit parameters: $(nfit)
Number of points in fit = $(npts)
Offset: $(offset) $(offset_error) offset_e*sqrt(n) = $(offset_esn) n = 9
PARAMETER                 Pre-fit                   Post-fit                  Uncertainty   Difference   Fit
-------------------------------------------------------------------------------------------------------------
F0                        1.0                       $(f0_post)                0.01          $(f0_post - 1.0) Y
F1                        2.0                       $(f1_post)                0.02          $(f1_post - 2.0) N
-------------------------------------------------------------------------------------------------------------
"""
end

function _write_residual_fixture(path::AbstractString, rows)
    open(path, "w") do io
        for row in rows
            println(io, join(row, " "))
        end
    end
    return path
end

function _synthetic_tim_entries()
    return [
        TimTOAEntry(1, 3, 100.0, 1400.0, 1.0, "be1"),
        TimTOAEntry(2, 4, 101.0, 1400.0, 1.5, "be2"),
        TimTOAEntry(3, 5, 102.0, 1400.0, 2.0, "be2"),
    ]
end

function _parsed_iteration(; kwargs...)
    out = parse_internal_iteration_tempo_output(_synthetic_iteration_block(; kwargs...), Tempo2)
    @assert !iserror(out.error)
    return out
end

@testset "Result and task assembly tests" begin
    @testset "TempoOutput parser covers richer iteration/result shapes" begin
        multi = string(
            "panic: synthetic preamble before first anchor\n",
            _synthetic_iteration_block(header = "Complete fit", chisq = 90.0, chisq_red = 9.0, f0_post = 1.2),
            "\n",
            _synthetic_iteration_block(
                header = "ss/fs",
                rms_pre = 11.0,
                rms_post = 5.5,
                rms_tn = nothing,
                chisq = 55.0,
                chisq_red = 5.5,
                pre_post = 2.0,
                f0_post = 1.4,
            ),
        )

        parsed = parse_tempo_output(multi, Tempo2)

        @test length(parsed) == 2
        @test parsed[1].error.error_type == :runtime_error
        @test parsed[1].basic !== nothing
        @test parsed[1].fit_parameters !== nothing
        @test !iserror(parsed[2].error)
        @test parsed[2].basic !== nothing
        @test isnan(parsed[2].basic.rms_tn_post_fit_residual_us)
        @test parsed[2][:F0] isa FitParameter
        @test parsed[2][:F0].post_fit == 1.4
        @test parsed[2][:chisqr_red] == 5.5

        missing_table = """
Complete fit
RMS pre-fit residual = 10.0 (us), RMS post-fit residual = 5.0 (us)
Fit Chisq = 10.0 Chisqr/nfree = 10.0/5 = 2.0 pre/post = 2.0
Number of fit parameters: 1
Number of points in fit = 5
Offset: 1.0 0.1 offset_e*sqrt(n) = 0.3 n = 9
"""
        missing_table_out = parse_internal_iteration_tempo_output(missing_table, Tempo2)
        @test missing_table_out.basic !== nothing
        @test missing_table_out.fit_parameters === nothing
        @test missing_table_out.error.error_type == :missing_fit_table

        crash_out = parse_internal_iteration_tempo_output("Segmentation fault\n", Tempo2)
        @test crash_out.error.error_type == :crash
        @test crash_out.basic === nothing
    end

    @testset "Residual/TIM assembly handles boundaries and safe combine failures" begin
        tim_entries = _synthetic_tim_entries()
        residuals = [
            TempoResidualEntry(100.0, 0.5, 1.0, 0.8, 1.0),
            TempoResidualEntry(101.0, 0.6, 1.5, 1.2, 0.0),
            TempoResidualEntry(102.0, 0.7, 2.0, 1.6, NaN),
        ]

        combined = combine_tim_and_residuals(
            tim_entries,
            residuals;
            time_start = 101.0,
            time_finish = 102.0,
        )

        @test getfield.(combined, :in_fit) == [false, true, true]
        @test getfield.(combined, :weight) == [1.0, 0.0, 0.0]
        @test all(isapprox.(getfield.(combined, :red_noise), [0.2, 0.3, 0.4]; atol = 1e-12))

        work_dir = mktempdir()
        good_residual_path = _write_residual_fixture(
            joinpath(work_dir, "residuals_good.dat"),
            [
                (100.0, 1e-6, 2e-6, 1.5e-6, 1e-6),
                (101.0, 1.1e-6, 2.2e-6, 1.6e-6, 2e-6),
                (102.0, 1.2e-6, 2.4e-6, 1.7e-6, 3e-6),
            ],
        )
        out = _parsed_iteration(pre_post = 1.0, f0_post = 1.8)
        iter_result = build_internal_iteration_result(
            out,
            good_residual_path,
            tim_entries;
            save_residuals = true,
            time_start = 101.0,
            time_finish = 102.0,
        )

        @test iter_result.residuals !== nothing
        @test length(iter_result.residuals) == 3
        @test iter_result.stats !== nothing
        @test !in_fit_equals_in_tim(iter_result.stats)
        @test iter_result.stats.in_fit.all.raw.n == 2
        @test iter_result.stats.in_tim.all.raw.n == 3
        @test iter_result.metadata[:stage_residual_read] == :ok
        @test iter_result.metadata[:stage_combine] == :ok
        @test iter_result.metadata[:stage_stats] == :ok
        @test iter_result.metadata[:iter_ok] == true
        @test iter_result.metadata[:combined_mjd_min] == 100.0
        @test iter_result.metadata[:combined_mjd_max] == 102.0

        bad_residual_path = _write_residual_fixture(
            joinpath(work_dir, "residuals_short.dat"),
            [(100.0, 1e-6, 2e-6, 1.5e-6, 1e-6)],
        )
        iter_fail = build_internal_iteration_result(
            out,
            bad_residual_path,
            tim_entries[1:2];
            save_residuals = true,
        )

        @test iter_fail.residuals !== nothing
        @test length(iter_fail.residuals) == 1
        @test iter_fail.stats === nothing
        @test iter_fail.metadata[:stage_residual_read] == :ok
        @test iter_fail.metadata[:stage_combine] == :failed
        @test iter_fail.metadata[:stage_stats] == :skipped
        @test iter_fail.metadata[:iter_ok] == false
        @test iter_fail.metadata[:iter_fail_stage] == :combine
    end

    @testset "GeneralTempoResult and BasicTempoTask assemble supported results" begin
        tim_entries = _synthetic_tim_entries()
        work_dir = mktempdir()
        residual_path = _write_residual_fixture(
            joinpath(work_dir, "residuals_builder.dat"),
            [
                (100.0, 1e-6, 2e-6, 1.5e-6, 1e-6),
                (101.0, 1.1e-6, 2.2e-6, 1.6e-6, 2e-6),
                (102.0, 1.2e-6, 2.4e-6, 1.7e-6, 3e-6),
            ],
        )

        iter1 = build_internal_iteration_result(
            _parsed_iteration(pre_post = 1.0, f0_post = 1.5),
            residual_path,
            tim_entries;
            time_start = 100.0,
            time_finish = 102.0,
        )
        iter2 = build_internal_iteration_result(
            _parsed_iteration(pre_post = 1.0, f0_post = 1.9, chisq = 120.0, chisq_red = 12.0),
            residual_path,
            tim_entries;
            time_start = 100.0,
            time_finish = 102.0,
        )

        built = build_general_tempo_result(
            [iter1, iter2];
            metadata = Dict{Symbol,Any}(:status => :ok, :success => true),
            metrics_hook = (final_iter, conv) -> Dict(
                :extra_metric => length(final_iter.output.fit_parameters),
                :nonfinite_metric => Inf,
                :pre_post_abs => abs(conv.pre_post_final - 1.0),
            ),
        )

        @test built.success
        @test built.status == :ok
        @test built.final_index == 2
        @test built.last_successful_index == 2
        @test built.final === iter2
        @test built.last_successful === iter2
        @test built.residual_stats === iter2.stats
        @test built.white_noise_fit === nothing
        @test built.param_estimates[:F0].value == 1.9
        @test built.metrics[:extra_metric] == 2.0
        @test isnan(built.metrics[:nonfinite_metric])
        @test built.metrics[:pre_post_abs] == 0.0
        @test result_metric(built, :missing_metric) |> isnan
        @test built.convergence.wrms.final_abs_delta == 0.0
        @test built.convergence.chisqr.final_abs_delta == 0.0

        fixture_paths = _result_fixture_paths()
        task_work_dir = mktempdir()
        settings = _basic_task_settings(task_work_dir)
        fake_script_dir = mktempdir()
        fake_script = _write_fake_basic_task_script(
            joinpath(fake_script_dir, "fake_tempo2_basic.sh"),
            fixture_paths,
        )

        task_result = withenv("TEMPO2_CMD" => fake_script) do
            run_task(BasicTempoTask(settings))
        end

        @test task_result isa GeneralTempoResult
        @test task_result.success
        @test task_result.status == :ok
        @test length(task_result.iterations) == 1
        @test task_result.final_index == 1
        @test task_result.last_successful_index == 1
        @test task_result.final === task_result.iterations[1]
        @test task_result.last_successful === task_result.iterations[1]
        @test task_result.final.residuals !== nothing
        @test task_result.final.stats !== nothing
        @test task_result.residual_stats === task_result.final.stats
        @test task_result.white_noise_fit === nothing
        @test task_result.final.metadata[:stage_stats] == :ok
        @test task_result.final.metadata[:stage_white_noise] == :skipped
        @test task_result.final.metadata[:residuals_stored] == true
        @test task_result.metadata[:run_status] == :ok
        @test task_result.metadata[:run_success] == true
        @test task_result.metadata[:status] == :ok
        @test task_result.metadata[:success] == true
        @test task_result.metadata[:n_iter] == 1
        @test task_result.metadata[:wn_enabled] == false
        @test task_result.metadata[:save_residuals] == true
        @test task_result.par_file_final !== nothing
        @test isapprox(task_result.metrics[:chi2_fit_basic], 1836.509373; rtol = 1e-12)
        @test isfinite(task_result.metrics[:wrms_fit])
        @test task_result.param_estimates[:F0].value ≈ 2.53871867521183
    end
end
