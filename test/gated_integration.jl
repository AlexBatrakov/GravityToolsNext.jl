using Test
using GravityToolsNext

const GTN = GravityToolsNext

function _gated_fixture_paths()
    fixture_root = normpath(joinpath(@__DIR__, "data", "ddstg_smoke"))
    return (
        root = fixture_root,
        par = joinpath(fixture_root, "j1141-6545_ddstg_gr_pdfb1_white_noise.par"),
        tim = joinpath(fixture_root, "j1141-6545_ddstg_gr_pdfb1_white_noise.tim"),
        out = "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par",
    )
end

function _copy_gated_inputs!(dest_dir::AbstractString)
    f = _gated_fixture_paths()
    mkpath(dest_dir)
    par_name = basename(f.par)
    tim_name = basename(f.tim)
    cp(f.par, joinpath(dest_dir, par_name); force=true)
    cp(f.tim, joinpath(dest_dir, tim_name); force=true)
    return (; par_name, tim_name)
end

function _truthy_env(value)
    value === nothing && return false
    normalized = lowercase(strip(String(value)))
    return normalized in ("1", "true", "yes", "on")
end

function _tempo2_integration_gate()
    enabled = _truthy_env(get(ENV, "GTN_ENABLE_TEMPO2_INTEGRATION", nothing))
    enabled || return (false, "set GTN_ENABLE_TEMPO2_INTEGRATION=1 to opt into real Tempo2 integration tests")

    fixture_paths = _gated_fixture_paths()
    for path in (fixture_paths.par, fixture_paths.tim)
        isfile(path) || return (false, "required smoke fixture is missing: $(path)")
    end

    try
        tempo = Tempo()
        tempo2 = Tempo2()
        validate(tempo)
        validate(tempo2)

        tempo_ddstg = joinpath(tempo_data_dir(tempo), "data_ddstg")
        tempo2_ddstg = joinpath(tempo_data_dir(tempo2), "data_ddstg")

        isdir(tempo_ddstg) || return (false, "TEMPO data_ddstg directory is missing: $(tempo_ddstg)")
        ispath(tempo2_ddstg) || return (false, "TEMPO2 data_ddstg link is missing: $(tempo2_ddstg)")
        realpath(tempo2_ddstg) == realpath(tempo_ddstg) ||
            return (false, "TEMPO2 data_ddstg link does not point at TEMPO data_ddstg")

        return (true, "ready")
    catch err
        return (false, sprint(showerror, err))
    end
end

function _integration_settings(
    work_dir::AbstractString;
    job_name::AbstractString = "tempo2-integration-smoke",
    white_noise_enabled::Bool = false,
    white_noise_scope::Symbol = :final,
    save_residuals::Bool = true,
)
    files = _copy_gated_inputs!(work_dir)
    return TempoRunSettings(
        work_dir = String(work_dir),
        par_input = files.par_name,
        tim_input = files.tim_name,
        par_output = _gated_fixture_paths().out,
        tempo_version = Tempo2(),
        flags = "",
        nits = 1,
        gain = 1.0,
        write_output = true,
        write_residuals = true,
        save_internal_iterations = false,
        save_residuals = save_residuals,
        white_noise_enabled = white_noise_enabled,
        white_noise_scope = white_noise_scope,
        work_mode = :jobdir,
        job_name = String(job_name),
        overwrite = :clean,
        layout = :split,
        keep_tmp_on_success = false,
        keep_tmp_on_error = true,
        write_manifest = false,
    )
end

@testset "Gated Tempo2 integration" begin
    enabled, reason = _tempo2_integration_gate()
    if !enabled
        @info "Skipping gated Tempo2 integration tests" reason=reason
        @test_skip enabled
    else
        @testset "run_tempo_parsed real Tempo2 smoke path" begin
            work_dir = mktempdir()
            settings = _integration_settings(
                work_dir;
                job_name = "tempo2-run-smoke",
                save_residuals = false,
            )
            run_out = run_tempo_parsed(settings)

            @test run_out isa GTN.TempoRunOutput
            @test run_out.success
            @test run_out.status == :ok
            @test run_out.exit_code == 0
            @test run_out.n_iter == 1
            @test run_out.artifacts.out_path !== nothing
            @test isfile(run_out.artifacts.out_path)
            @test isfile(run_out.artifacts.par_out_path)
            @test any(path -> path !== nothing && isfile(path), run_out.artifacts.residual_paths)
            @test isdir(run_out.artifacts.run_cwd)
            @test run_out.parsed[end].basic !== nothing
            @test run_out.parsed[end].fit_parameters !== nothing

            GTN.cleanup_run!(run_out, settings)
            @test !isdir(run_out.artifacts.run_cwd)
        end

        @testset "BasicTempoTask real Tempo2 smoke path" begin
            work_dir = mktempdir()
            result = run_task(BasicTempoTask(_integration_settings(
                work_dir;
                job_name = "tempo2-basic-task-smoke",
            )))

            @test result isa GeneralTempoResult
            @test result.success
            @test result.status == :ok
            @test result.metadata[:success] == true
            @test result.metadata[:status] == :ok
            @test result.metadata[:run_success] == true
            @test result.metadata[:run_status] == :ok
            @test result.metadata[:n_iter] == 1
            @test result.metadata[:tempo_ver] == Tempo2
            @test result.par_file_final !== nothing
            @test result.final.output.basic !== nothing
            @test result.final.output.fit_parameters !== nothing
            @test result.final.stats !== nothing
            @test result.final.residuals !== nothing
            @test length(result.final.residuals) > 0
            @test haskey(result.param_estimates, :F0)
            @test isfinite(result.param_estimates[:F0].value)
            @test isfinite(result.metrics[:chi2_fit_basic])
            @test isfinite(result.metrics[:wrms_fit])
            @test result.metadata[:out_path] !== nothing
            @test isfile(result.metadata[:out_path])
            @test isfile(result.metadata[:par_out_path])
            @test isdir(result.metadata[:job_root])
            @test isdir(result.metadata[:input_dir])
            @test isdir(result.metadata[:output_dir])
            @test !isdir(result.metadata[:run_cwd])
        end

        @testset "BasicTempoTask white-noise Tempo2 smoke path" begin
            work_dir = mktempdir()
            result = run_task(BasicTempoTask(_integration_settings(
                work_dir;
                job_name = "tempo2-white-noise-smoke",
                white_noise_enabled = true,
                white_noise_scope = :final,
                save_residuals = false,
            )))

            @test result isa GeneralTempoResult
            @test result.success
            @test result.status == :ok
            @test result.metadata[:wn_enabled] == true
            @test result.metadata[:wn_scope] == :final
            @test result.final.stats !== nothing
            @test result.white_noise_fit !== nothing
            @test result.final.metadata[:stage_white_noise] == :ok
            @test !isempty(result.white_noise_fit.by_backend)
            @test isempty(result.white_noise_fit.failed_backends)
            @test result.white_noise_fit.global_stats.n > 0
            @test isfinite(result.metrics[:ad_white_fit])
            @test !isdir(result.metadata[:run_cwd])
        end
    end
end
