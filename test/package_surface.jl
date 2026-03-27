using Test
using GravityToolsNext

@testset "Supported package surface" begin
    exported = Set(names(GravityToolsNext))

    for sym in (
        :get_tempo_command,
        :parse_internal_iteration_tempo_output,
        :is_converged,
        :task_copy_with,
        :task_derive_par_output,
        :save_result_jld2,
        :Adaptive2DGridTask,
        :GridWorkspaceOptions,
    )
        @test sym in exported
    end

    @test :has_converged ∉ exported
    @test :print_backend_table ∉ exported

    fixture_root = normpath(joinpath(@__DIR__, "data", "ddstg_smoke"))
    settings = TempoRunSettings(
        work_dir = fixture_root,
        par_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.par",
        tim_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.tim",
        par_output = "j1141-6545_ddstg_gr_pdfb1_white_noise_out.par",
        tempo_version = Tempo2(fixture_root),
        work_mode = :jobdir,
        job_name = "surface_test",
        overwrite = :reuse,
        layout = :split,
    )

    task = BasicTempoTask(settings)
    copied = task_copy_with(task; job_name = "surface_test_copy")

    @test copied isa BasicTempoTask
    @test copied.settings.workspace.job_name == "surface_test_copy"
    @test task_derive_par_output(task, "node001") == "j1141-6545_ddstg_gr_pdfb1_white_noise_node001_out.par"
end
