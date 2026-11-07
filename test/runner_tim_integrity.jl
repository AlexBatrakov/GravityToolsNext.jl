using Test, GravityToolsNext, SHA

runner_tim_quote(s)="'"*replace(s,"'"=>"'\\''")*"'"
const RUNNER_TIM_ORIGINAL="FORMAT 1\nobs 1400 53590.0 2 pks -be A\n"
const RUNNER_TIM_CHANGED=replace(RUNNER_TIM_ORIGINAL,"53590.0"=>"53590.1")

function runner_tim_fixture(root;layout=:split,temp_dir=nothing,io_mirror=:none)
    write(joinpath(root,"model.par"),"PSRJ PROBE\nF0 2.5\nPEPOCH 53590\n")
    write(joinpath(root,"input observations.tim"),RUNNER_TIM_ORIGINAL)
    marker=joinpath(root,"engine-started")
    fake=joinpath(root,"read-tim.sh")
    write(fake,"#!/bin/sh\nset -eu\nprintf started >> "*runner_tim_quote(marker)*"\ncat \"\$3\"\n")
    chmod(fake,0o755)
    settings=TempoRunSettings(;work_dir=root,par_input="model.par",tim_input="input observations.tim",
        par_output="model_out.par",tempo_version=Tempo2(root),work_mode=:jobdir,job_name="probe",
        overwrite=:error,layout,temp_dir,io_mirror,link_tim=false,cleanup_before_run=false,
        write_output=false,write_residuals=false)
    return (;settings,fake,marker)
end

function runner_tim_copy!(job)
    actual=joinpath(job.run_cwd,basename(job.tim_in_path))
    if actual!=job.tim_in_path
        bytes=read(actual);rm(actual);write(actual,bytes)
    end
    return actual
end

@testset "Executed TIM logical path and archive consistency" begin
    G=GravityToolsNext
    for (layout,temp_dir,io_mirror) in [(:flat,nothing,:none),(:split,nothing,:none),
            (:split,"nested/work",:full),(:split,"nested/work",1)],
            storage in (:linked,:copy),entry in (:job,:prepared)
        mktempdir() do root
            f=runner_tim_fixture(root;layout,temp_dir,io_mirror)
            withenv("TEMPO2_CMD"=>f.fake) do
                job=G.materialize_job(f.settings)
                actual=storage==:copy ? runner_tim_copy!(job) : joinpath(job.run_cwd,basename(job.tim_in_path))
                if storage==:linked && layout==:split
                    @test islink(actual)
                end
                prepared=G._prepare_tempo_run(job)
                @test prepared.tim_path==abspath(job.run_cwd,prepared.cmd.exec[4])
                @test prepared.cmd.dir==job.run_cwd
                @test prepared.tim_sha256==bytes2hex(sha256(read(actual)))
                @test read(actual)==read(job.tim_in_path)
                @test prepared.diagnostics[:executed_tim_path]==actual
                @test prepared.diagnostics[:archive_tim_path]==job.tim_in_path
                @test prepared.diagnostics[:archive_tim_sha256]==prepared.tim_sha256
                out,err,code=entry==:job ? run_tempo_raw(job) : run_tempo_raw(prepared)
                @test code==0 && isempty(err) && out==RUNNER_TIM_ORIGINAL
                @test read(f.marker,String)=="started"
                @test job.notes[:run_diagnostics][:executed_tim_sha256]==prepared.tim_sha256
            end
        end
    end
    for change in (:modified_copy,:replace_symlink,:retarget,:missing,:broken_link,:archive_edit,:both_changed,:par_edit)
        mktempdir() do root
            f=runner_tim_fixture(root)
            withenv("TEMPO2_CMD"=>f.fake) do
                job=G.materialize_job(f.settings)
                actual=joinpath(job.run_cwd,basename(job.tim_in_path))
                change in (:modified_copy,:archive_edit,:both_changed) && runner_tim_copy!(job)
                prepared=G._prepare_tempo_run(job)
                oldtarget=realpath(actual)
                if change==:par_edit
                    write(prepared.init_path,read(prepared.init_path,String)*"F1 0\n")
                elseif change==:archive_edit
                    write(job.tim_in_path,RUNNER_TIM_CHANGED)
                elseif change==:both_changed
                    write(actual,RUNNER_TIM_CHANGED);write(job.tim_in_path,RUNNER_TIM_CHANGED)
                elseif change==:modified_copy
                    write(actual,RUNNER_TIM_CHANGED)
                else
                    rm(actual)
                    if change==:replace_symlink
                        write(actual,RUNNER_TIM_CHANGED)
                    elseif change==:retarget
                        other=joinpath(root,"other.tim");write(other,RUNNER_TIM_CHANGED);symlink(other,actual)
                        @test isfile(oldtarget) && read(oldtarget,String)==RUNNER_TIM_ORIGINAL
                    elseif change==:broken_link
                        symlink(joinpath(root,"absent.tim"),actual)
                    end
                end
                archive_before=read(job.tim_in_path)
                err=try;run_tempo_raw(prepared);nothing;catch e;e;end
                @test err isa ArgumentError
                @test occursin(change==:par_edit ? "PAR" : change==:archive_edit ? "archive TIM" : "executed TIM",sprint(showerror,err))
                @test !ispath(f.marker) # No process launched, even if a fake would exit zero.
                @test read(job.tim_in_path)==archive_before
                @test !haskey(job.notes,:run_diagnostics)
                if change in (:modified_copy,:replace_symlink,:retarget,:both_changed)
                    @test read(actual,String)==RUNNER_TIM_CHANGED # No silent repair.
                end
            end
        end
    end
    for entry in (:prepare,:raw),missing in (false,true)
        mktempdir() do root
            f=runner_tim_fixture(root)
            withenv("TEMPO2_CMD"=>f.fake) do
                job=G.materialize_job(f.settings);actual=joinpath(job.run_cwd,basename(job.tim_in_path));rm(actual)
                missing || write(actual,RUNNER_TIM_CHANGED)
                err=try;entry==:prepare ? G._prepare_tempo_run(job) : run_tempo_raw(job);nothing;catch e;e;end
                @test err isa ArgumentError
                @test occursin(missing ? "executed TIM" : "differs from archive",sprint(showerror,err))
                @test !ispath(f.marker)
                @test read(job.tim_in_path,String)==RUNNER_TIM_ORIGINAL
                @test !isfile(joinpath(job.run_cwd,"model_init.par"))
                missing || @test read(actual,String)==RUNNER_TIM_CHANGED
            end
        end
    end
    # chmod cannot deny root's CAP_DAC_OVERRIDE. Report that capability explicitly
    # instead of claiming an unreadable-file guard was exercised under root.
    if Sys.isunix()
        mktempdir() do root
            f=runner_tim_fixture(root)
            withenv("TEMPO2_CMD"=>f.fake) do
                job=G.materialize_job(f.settings);actual=runner_tim_copy!(job);p=G._prepare_tempo_run(job)
                chmod(actual,0o000)
                try
                    readable=try;read(actual);true;catch;false;end
                    if readable
                        @info "Unreadable-file permission probe unavailable: process can read mode-000 files"
                        @test_skip !readable
                    else
                        @test_throws ArgumentError run_tempo_raw(p)
                    end
                    @test !ispath(f.marker)
                finally
                    chmod(actual,0o600)
                end
            end
        end
    end
end

# Uses the existing neutral fake-output/Basic fixtures from result_task_assembly.
@testset "Parsed and Basic executed TIM consumers" begin
    G=GravityToolsNext
    for consumer in (:parsed,:basic),outcome in (:ok,:nonzero,:timeout,:changed_input)
        mktempdir() do root
            settings=copy_with(_basic_task_settings(root);timeout_s=outcome==:timeout ? 1 : nothing,
                keep_tmp_on_error=true,keep_tmp_on_success=true)
            marker=joinpath(root,"engine-started")
            fake=_write_fake_basic_task_script(joinpath(root,"fake.sh"),_result_fixture_paths())
            text=read(fake,String)
            text=replace(text,"set -eu\n"=>"set -eu\nprintf started >> "*runner_tim_quote(marker)*"\ninput=\"\$3\"\n")
            write(fake,text)
            open(fake,"a") do io
                print(io,outcome==:nonzero ? "\nexit 7\n" : outcome==:timeout ? "\nexec sleep 30\n" : outcome==:changed_input ? "\nprintf '\\nC changed by fake\\n' >> \"\$input\"\n" : "\nexit 0\n")
            end
            run_consumer()=consumer==:parsed ? run_tempo_parsed(settings) : run_task(BasicTempoTask(settings))
            withenv("TEMPO2_CMD"=>fake) do
                if outcome==:changed_input
                    @test_throws ArgumentError run_consumer()
                else
                    result=run_consumer()
                    @test result.success==(outcome==:ok)
                    @test result.status==(outcome==:ok ? :ok : :engine_failed)
                    if consumer==:basic
                        @test result.metadata[:run_success]==(outcome==:ok)
                        @test isdir(result.metadata[:run_cwd])
                    else
                        @test result.exit_code==(outcome==:ok ? 0 : 1)
                        @test isdir(result.artifacts.run_cwd)
                        @test read(result.artifacts.tim_path)==read(joinpath(result.artifacts.run_cwd,basename(result.artifacts.tim_path)))
                    end
                end
                @test read(marker,String)=="started"
            end
        end
    end
    for consumer in (:parsed,:basic)
        mktempdir() do root
            settings=copy_with(_basic_task_settings(root);cleanup_before_run=false,overwrite=:reuse)
            job=G.materialize_job(settings);actual=runner_tim_copy!(job)
            write(actual,RUNNER_TIM_CHANGED);archive=read(job.tim_in_path)
            marker=joinpath(root,"engine-started");fake=joinpath(root,"never.sh")
            write(fake,"#!/bin/sh\nprintf started > "*runner_tim_quote(marker)*"\n");chmod(fake,0o755)
            withenv("TEMPO2_CMD"=>fake) do
                @test_throws ArgumentError consumer==:parsed ? run_tempo_parsed(settings) : run_task(BasicTempoTask(settings))
                @test !ispath(marker)
                @test read(actual,String)==RUNNER_TIM_CHANGED && read(job.tim_in_path)==archive
            end
        end
    end
end
