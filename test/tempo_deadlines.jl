using Test, GravityToolsNext

@testset "Finite engine deadlines" begin
    G=GravityToolsNext
    @test_throws ArgumentError WorkspaceOptions(timeout_s=0)
    @test_throws ArgumentError WorkspaceOptions(timeout_s=-1)
    @test_throws ArgumentError WorkspaceOptions(timeout_s=true)
    if Sys.isunix()
        ok=G._capture_tempo_process(`sh -c 'printf success; printf warning >&2'`,2)
        @test ok.code==0 && ok.status==:ok && !ok.timed_out
        @test ok.stdout=="success" && ok.stderr=="warning"
        failure=G._capture_tempo_process(`sh -c 'printf partial; printf problem >&2; exit 7'`,2)
        @test failure.code==1 && failure.native_exit==7 && failure.status==:engine_failed
        @test startswith(failure.stdout,"partial") && failure.stderr=="problem"
        normal=G._capture_tempo_process(`sh -c 'printf old; exit 3'`,nothing)
        @test normal.code==1 && startswith(normal.stdout,"old")
        @test !normal.timed_out
        mktempdir() do root
            # Both shell and child ignore TERM; the owned group must receive KILL.
            childfile=joinpath(root,"child.pid")
            script=joinpath(root,"hang.sh")
            write(script,"""#!/bin/sh
trap '' TERM
printf before-timeout
printf error-before-timeout >&2
sh -c 'trap "" TERM; exec sleep 30' &
printf '%s' "\$!" > "\$1"
wait
""")
            started=time()
            result=G._capture_tempo_process(`sh $script $childfile`,1)
            @test time()-started<5
            @test result.code==1 && result.timed_out && result.status==:timeout
            @test startswith(result.stdout,"before-timeout")
            @test result.stderr=="error-before-timeout"
            child=parse(Int,read(childfile,String))
            # A descendant may remain as a zombie until the OS reaper collects it;
            # it must never remain executing. Julia reaps its direct child.
            if Sys.islinux() && isfile("/proc/$child/stat")
                @test split(read("/proc/$child/stat",String))[3]=="Z"
            else
                @test ccall(:kill,Cint,(Cint,Cint),child,0)!=0
            end
        end
        descendants=G._capture_tempo_process(`sh -c 'sleep 30 & printf child-started'`,2)
        @test descendants.code==1 && descendants.status==:descendants_terminated
        @test startswith(descendants.stdout,"child-started")
        missing=G._capture_tempo_process(`/definitely/no/engine`,1)
        @test missing.code==1 && missing.status==:engine_failed
        @test occursin("ERROR",missing.stdout)
        # A separate process group is untouched by timeout cleanup.
        unrelated=run(`sleep 30`;wait=false)
        try
            G._capture_tempo_process(`sh -c 'sleep 30'`,1)
            @test process_running(unrelated)
        finally
            kill(unrelated);wait(unrelated)
        end
    end
end
