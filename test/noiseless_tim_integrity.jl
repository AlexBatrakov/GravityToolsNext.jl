using Test, GravityToolsNext, SHA

@testset "Template actual/archive/candidate TIM identity" begin
    G=GravityToolsNext
    mktempdir() do root
        source,par=template_test_inputs(root);marker=joinpath(root,"started")
        fake=template_fake_executable(root,"ok";marker)
        t=withenv("TEMPO2_CMD"=>fake) do
            create_noiseless_template(source,par;workspace=joinpath(root,"job"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
        end
        @test t.metadata["format_version"]==3 && t.metadata["initial_margin_policy"]=="diagnostic_only"
        @test read(marker,String)==repeat("started",length(t.metadata["history"]))
        for h in t.metadata["history"]
            q=h["query"];actual=h["executed_input_tim"];archive=h["input_tim"]
            bytes=t.artifacts[q*"/candidate.tim"];hash=bytes2hex(sha256(bytes))
            @test t.artifacts[actual]==t.artifacts[archive]==bytes
            invocation=String(copy(t.artifacts[q*"/invocation.json"]))
            execution=String(copy(t.artifacts[q*"/execution.json"]))
            for key in ["tim_sha256","executed_tim_sha256","archive_tim_sha256"]
                @test occursin("\""*key*"\":"*G._to_json(hash),invocation)
            end
            for (key,path) in [("executed_tim_path",actual),("archive_tim_path",archive)]
                original=joinpath(root,"job",basename(q),relpath(path,q))
                @test occursin("\""*key*"\":"*G._to_json(original),invocation)
                @test occursin("\""*key*"\":"*G._to_json(original),execution)
            end
        end
        final=last(t.metadata["history"])
        @test t.artifacts[final["executed_input_tim"]]==t.artifacts["final.tim"]
        bundle=save_noiseless_template(t,joinpath(root,"bundle"))
        rm(joinpath(root,"job");recursive=true);rm(fake)
        loaded=withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(bundle);end
        @test loaded.metadata==t.metadata && loaded.artifacts==t.artifacts
        @test loaded.lineage==t.lineage
        save_noiseless_template(loaded,joinpath(root,"resaved"))
        incomplete=deepcopy(t)
        delete!(first(incomplete.metadata["history"]),"executed_input_tim")
        @test_throws ArgumentError G._validate_noiseless_template(incomplete)
        for which in (:actual,:archive,:history_path,:archive_alias)
            bad=deepcopy(t);h=first(bad.metadata["history"])
            if which==:history_path
                h["executed_input_tim"]="../invalid.tim"
            elseif which==:archive_alias
                h["executed_input_tim"]=h["input_tim"]
            else
                path=h[which==:actual ? "executed_input_tim" : "input_tim"]
                append!(bad.artifacts[path],codeunits("C changed snapshot\n"))
                bad.metadata["artifact_sha256"][path]=bytes2hex(sha256(bad.artifacts[path]))
            end
            @test_throws ArgumentError save_noiseless_template(bad,joinpath(root,"bad"))
            @test !ispath(joinpath(root,"bad"))
        end
    end
    for mode in ["tim_replace","tim_retarget","archive_edit","candidate_edit"]
        mktempdir() do root
            source,par=template_test_inputs(root);before=G._tim_render(source.document);par_before=read(par)
            marker=joinpath(root,"started");fake=template_fake_executable(root,mode;marker)
            err=try
                withenv("TEMPO2_CMD"=>fake) do
                    create_noiseless_template(source,par;workspace=joinpath(root,"job"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
                end
                nothing
            catch e;e;end
            @test err isa IdealizationError && err.status==:validation_failed
            @test read(marker,String)=="started"
            @test G._tim_render(source.document)==before && read(par)==par_before
            @test !isdir(joinpath(root,"job/q002"))
            @test isfile(joinpath(root,"job/failure.txt"))
            @test occursin(mode=="archive_edit" ? "archive TIM" : mode=="candidate_edit" ? "differs from candidate" : "executed TIM",sprint(showerror,err))
        end
    end
end
