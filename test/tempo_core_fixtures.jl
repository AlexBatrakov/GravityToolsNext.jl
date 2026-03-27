using Test
using GravityToolsNext

const GTN = GravityToolsNext

function _tempo_core_fixture_paths()
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

function _copy_smoke_inputs!(dest_dir::AbstractString)
    f = _tempo_core_fixture_paths()
    mkpath(dest_dir)
    par_name = basename(f.par)
    tim_name = basename(f.tim)
    cp(f.par, joinpath(dest_dir, par_name); force=true)
    cp(f.tim, joinpath(dest_dir, tim_name); force=true)
    return (; par_name, tim_name)
end

function _base_settings(work_dir::AbstractString)
    files = _copy_smoke_inputs!(work_dir)
    return TempoRunSettings(
        work_dir = String(work_dir),
        par_input = files.par_name,
        tim_input = files.tim_name,
        par_output = "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par",
        tempo_version = Tempo2(_tempo_core_fixture_paths().root),
        work_mode = :jobdir,
        job_name = "tempo-core-job",
        overwrite = :clean,
        layout = :split,
    )
end

function _param_float(params::Vector{TempoParameter}, name::Symbol)
    p = get_param(params, name)
    p === nothing && return nothing
    return value_as_float(p)
end

function _write_fake_tempo2_script(path::AbstractString, mode::AbstractString, fixture_paths)
    script = replace(raw"""#!/bin/sh
set -eu
mode="__MODE__"
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

copy_outpar() {
  if [ -n "$outpar" ]; then
    cp "$fixture_par_out" "$outpar"
  fi
}

copy_residual() {
  cp "$fixture_residual" "residuals.dat"
}

case "$mode" in
  ok)
    copy_outpar
    copy_residual
    cat "$fixture_out"
    ;;
  parse_failed)
    copy_outpar
    copy_residual
    printf '%s\\n' 'Error: synthetic parser failure'
    ;;
  files_missing)
    copy_outpar
    cat "$fixture_out"
    ;;
  engine_failed)
    printf '%s\\n' 'synthetic engine failure'
    exit 1
    ;;
  *)
    printf '%s\n' "Unknown GTN_FAKE_MODE=$mode" >&2
    exit 2
    ;;
esac
""",
        "__MODE__" => String(mode),
        "__FIXTURE_OUT__" => String(fixture_paths.smoke_out),
        "__FIXTURE_PAR_OUT__" => String(fixture_paths.smoke_par_out),
        "__FIXTURE_RESIDUAL__" => String(fixture_paths.smoke_residual),
    )
    write(path, script)
    chmod(path, 0o755)
    return path
end

@testset "TempoCore fixture tests" begin
    fixture_paths = _tempo_core_fixture_paths()

    @testset "TempoRunSettings validation and copy/update" begin
        settings = TempoRunSettings(
            work_dir = fixture_paths.root,
            par_input = basename(fixture_paths.par),
            tim_input = basename(fixture_paths.tim),
            par_output = "smoke_out.par",
            tempo_version = Tempo2(fixture_paths.root),
            override_params = [TP("F0", 1.0; flag=1), TP("F1", 2.0; flag=0)],
            white_noise_enabled = false,
            white_noise_scope = :final,
            work_mode = :jobdir,
            job_name = "settings-base",
            overwrite = :reuse,
            layout = :split,
            io_mirror = :none,
            verbosity = 1,
        )

        @test default_par_output("nested/example.PAR") == "example_out.par"
        @test validate(settings)
        @test validate_inputs_exist(settings)
        @test settings.logging.verbosity == :warn

        copied = copy_with(
            settings;
            work_dir = joinpath(fixture_paths.root, "settings-copy"),
            par_output = "copy_out.par",
            white_noise_enabled = true,
            white_noise_scope = :all,
            job_name = "settings-copy",
            overwrite = :unique,
            io_mirror = (:depth_minus, 1),
            verbosity = 3,
            override_params_delete = :F1,
            override_params_upsert = [TP("F0", 3.0; flag=-1), TP("F2", 4.0; flag=1)],
        )

        @test copied.paths.work_dir == joinpath(fixture_paths.root, "settings-copy")
        @test copied.paths.par_output == "copy_out.par"
        @test copied.analysis.enabled == true
        @test copied.analysis.scope == :all
        @test copied.workspace.job_name == "settings-copy"
        @test copied.workspace.overwrite == :unique
        @test copied.workspace.io_mirror == (:depth_minus, 1)
        @test copied.logging.verbosity == :debug
        @test _param_float(settings.modifiers.override_params, :F0) == 1.0
        @test _param_float(settings.modifiers.override_params, :F1) == 2.0
        @test _param_float(copied.modifiers.override_params, :F0) == 3.0
        @test get_param(copied.modifiers.override_params, :F1) === nothing
        @test _param_float(copied.modifiers.override_params, :F2) == 4.0

        cleared = copy_with(
            settings;
            override_params_clear = true,
            override_params_upsert = [TP("DM", 116.5; flag=0)],
        )
        @test length(cleared.modifiers.override_params) == 1
        @test _param_float(cleared.modifiers.override_params, :DM) == 116.5

        @test_throws ErrorException RunPaths(
            "relative/workdir",
            basename(fixture_paths.par),
            "out.par",
            basename(fixture_paths.tim),
        )
        @test_throws ErrorException RunPaths(
            fixture_paths.root,
            fixture_paths.par,
            "out.par",
            basename(fixture_paths.tim),
        )
        @test_throws ErrorException RunPaths(
            fixture_paths.root,
            basename(fixture_paths.par),
            "nested/out.par",
            basename(fixture_paths.tim),
        )
        @test_throws ErrorException InputModifiers(TempoParameter[], 10.0, 9.0, false)
        @test_throws ErrorException WorkspaceOptions(io_mirror = (:bad_mode, 1))
        @test_throws ErrorException EngineOptions(Tempo2(fixture_paths.root), " --flag\n--bad", 1, 1.0)
    end

    @testset "TempoParFile update coupling stays cheap and direct" begin
        pf_dir = mktempdir()
        pf = TempoParFile("coupling_test.par", pf_dir)
        upsert_params!(pf, [
            TP("F0", 10.0; flag=1),
            TP("F1", 1.0; flag=1),
            TP("D", 2.0; flag=0),
            TP("DDOT", 1.0; flag=0),
        ])

        update_par_file!(
            pf;
            override_params = [TP("DDOT", 5.0; flag=0)],
            couple_f1_to_ddot = true,
        )

        @test value_as_float(get_param(pf, :DDOT)) == 5.0
        @test value_as_float(get_param(pf, :F1)) == -19.0

        pf2 = copy_par_file(pf; new_name = "coupling_test_2.par", deep_copy = true)
        update_par_file!(
            pf2;
            override_params = [TP("DDOT", 9.0; flag=0), TP("F1", 123.0; flag=1)],
            couple_f1_to_ddot = true,
        )
        @test value_as_float(get_param(pf2, :F1)) == 123.0
    end

    @testset "TempoCore materialization and cleanup" begin
        work_dir = mktempdir()
        settings = copy_with(
            _base_settings(work_dir);
            temp_dir = joinpath("tmp", "deep", "run"),
            io_mirror = (:depth_minus, 1),
        )

        stale_dir = joinpath(work_dir, "tempo-core-job")
        mkpath(stale_dir)
        write(joinpath(stale_dir, "stale.txt"), "stale")

        job = GTN.materialize_job(settings)

        @test job isa GTN.MaterializedJob
        @test job.job_root == joinpath(work_dir, "tempo-core-job")
        @test !isfile(joinpath(job.job_root, "stale.txt"))
        @test job.input_dir == joinpath(job.job_root, "input")
        @test job.output_dir == joinpath(job.job_root, "output")
        @test job.tmp_dir == joinpath(job.job_root, "tmp", "deep", "run")
        @test job.run_cwd == job.tmp_dir
        @test job.par_in_path == joinpath(job.input_dir, "tmp", "deep", "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_in.par")
        @test job.par_out_path == joinpath(job.output_dir, "tmp", "deep", "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_out.par")
        @test isfile(job.par_in_path)
        @test isfile(job.tim_in_path)
        @test isfile(joinpath(job.run_cwd, basename(job.tim_in_path)))
        @test isfile(joinpath(job.run_cwd, "j1141-6545_ddstg_gr_pdfb1_white_noise_smoke_in.par"))
        @test job.notes[:snapshot_par] == true
        @test job.notes[:original_par_archived] == true

        artifacts = RunArtifacts(
            job.job_root,
            job.run_cwd,
            job.input_dir,
            job.output_dir,
            job.tim_in_path,
            job.par_out_path,
            nothing,
            Union{Nothing,String}[],
        )
        GTN.cleanup_run!(artifacts, settings; success=true)
        @test !isdir(job.run_cwd)
        @test isdir(job.input_dir)
        @test isdir(job.output_dir)

        flat_dir = mktempdir()
        flat_settings = copy_with(
            _base_settings(flat_dir);
            work_mode = :inplace,
            layout = :flat,
            cleanup_before_run = false,
        )
        flat_job = GTN.materialize_job(flat_settings)
        write(joinpath(flat_job.job_root, "base_derivatives.dat"), "old")
        write(joinpath(flat_job.job_root, "residuals_1.dat"), "old")
        write(joinpath(flat_job.job_root, basename(flat_job.par_out_path)), "old")
        write(joinpath(flat_job.job_root, "keep_me.txt"), "keep")

        cleanup_old_tempo_files(flat_job)

        @test !isfile(joinpath(flat_job.job_root, "base_derivatives.dat"))
        @test !isfile(joinpath(flat_job.job_root, "residuals_1.dat"))
        @test !isfile(joinpath(flat_job.job_root, basename(flat_job.par_out_path)))
        @test isfile(joinpath(flat_job.job_root, "keep_me.txt"))
    end

    @testset "TempoCore parser and status fixtures" begin
        smoke_text = read(fixture_paths.smoke_out, String)
        parsed = parse_tempo_output(smoke_text, Tempo2)

        @test length(parsed) == 1
        @test !iserror(parsed[end].error)
        @test parsed[end].basic !== nothing
        @test isapprox(parsed[end].basic.chisqr, 1836.509373; rtol=1e-12)
        @test parsed[end].basic.number_of_fit_parameters == 17
        @test parsed[end].fit_parameters !== nothing
        @test length(parsed[end].fit_parameters) > 10
        f0 = only(filter(p -> p.name_symbol == :F0, parsed[end].fit_parameters))
        @test isapprox(f0.post_fit, 2.53871867521183; rtol=1e-12)

        err_section = "Error: lost connection to synthetic backend\n"
        err_out = parse_internal_iteration_tempo_output(err_section, Tempo2)
        @test err_out.error.error_type == :lost_connection
        @test err_out.basic === nothing

        nan_table = """
PARAMETER                 Pre-fit                   Post-fit                  Uncertainty   Difference   Fit
-------------------------------------------------------------------------------------------------------------
F0                        1                         1                         NaN           0            Y
-------------------------------------------------------------------------------------------------------------
"""
        params, perr = parse_fit_parameters(nan_table, Tempo2)
        @test params === nothing
        @test perr.error_type == :nan_in_fit_param

        function run_fake_status(mode::AbstractString)
            work_dir = mktempdir()
            settings = _base_settings(work_dir)
            fake_script_dir = mktempdir()
            fake_script = _write_fake_tempo2_script(
                joinpath(fake_script_dir, "fake_tempo2.sh"),
                mode,
                fixture_paths,
            )
            return withenv(
                "TEMPO2_CMD" => fake_script,
            ) do
                run_tempo_parsed(settings)
            end
        end

        ok_out = run_fake_status("ok")
        @test ok_out.success
        @test ok_out.status == :ok
        @test ok_out.exit_code == 0
        @test ok_out.n_iter == 1
        @test isfile(ok_out.artifacts.par_out_path)
        @test ok_out.artifacts.out_path !== nothing
        @test isfile(ok_out.artifacts.out_path)
        @test ok_out.artifacts.residual_paths[1] !== nothing
        @test isfile(ok_out.artifacts.residual_paths[1])

        parse_failed = run_fake_status("parse_failed")
        @test !parse_failed.success
        @test parse_failed.status == :parse_failed
        @test parse_failed.exit_code == 0
        @test parse_failed.parsed[end].error.error_type == :explicit_error

        files_missing = run_fake_status("files_missing")
        @test !files_missing.success
        @test files_missing.status == :files_missing
        @test files_missing.exit_code == 0
        @test files_missing.artifacts.residual_paths[1] === nothing

        engine_failed = run_fake_status("engine_failed")
        @test !engine_failed.success
        @test engine_failed.status == :engine_failed
        @test engine_failed.exit_code == 1
    end
end
