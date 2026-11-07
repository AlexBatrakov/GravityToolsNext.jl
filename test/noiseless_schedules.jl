using Test, GravityToolsNext, SHA, JLD2

function schedule_test_source(;mode=nothing)
    src=regular_schedule("55825.0";count=4,cadence_s=1.03,observatory="pks",frequency_mhz=1400,
        uncertainty_us=[7,2,3,5],backend=["A","B","A","B"],mode)
    # Selection retains original generated IDs and moves the phase anchor to row 2.
    return select_toas(src;method=:time_window,time_start="55825.00001",time_finish="55826.0")
end

@testset "Generated template settings and effective modes" begin
    G=GravityToolsNext
    for kw in [(max_bootstrap_corrections=4,),(max_bootstrap_corrections=true,),
               (mode_switch_tolerance_s=1e-12,),(mode_switch_tolerance_s=big"1e-400",)]
        @test_throws ArgumentError IdealizationSettings(;kw...)
    end
    @test IdealizationSettings(max_bootstrap_corrections=0).max_bootstrap_corrections==0
    mktempdir() do root
        _,par=template_test_inputs(root)
        for mode in [nothing,0,1]
            src=schedule_test_source(;mode)
            before=G._tim_render(src.document)
            @test src.document.origin==:generated && first(src.lineage).source_ordinal==2
            @test !occursin(r"(?m)^MODE\s",G._template_evaluator_text(src))
            @test G._tim_render(src.document)==before
            fake=template_fake_executable(root,"ok")
            t=withenv("TEMPO2_CMD"=>fake) do
                create_noiseless_template(src,par;workspace=joinpath(root,"mode-$(mode)"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
            end
            @test t.lineage==src.lineage && t.metadata["source_origin"]=="generated"
            @test t.metadata["source_ids"]==[G._template_lineage_record(x) for x in src.lineage]
            b=t.metadata["bootstrap"]
            @test b["anchor_source_id"]["source_ordinal"]==2
            @test b["anchor_epoch_mjd"]==string(first(toa_entries(src.document)).epoch_mjd)
            @test b["entrance_limit_s"]==2e-7 # sigma_min is 2 us in this selection.
            @test all(h->h["pulse_numbers"]==b["initial_pulse_numbers"],t.metadata["history"])
            @test [h["mode"] for h in t.metadata["history"]]==[0,0,1]
            @test parse(BigFloat,b["mode_switch_max_difference_s"])==0
            @test occursin("JUMP -be A",String(copy(t.artifacts["reference.par"])))
            @test G._template_reference_mode(String(copy(t.artifacts["bootstrap.par"])),1)==String(copy(t.artifacts["reference.par"]))
            @test G._tim_render(src.document)==before && read(par,String)==TEMPLATE_TEST_PAR
            bundle=save_noiseless_template(t,joinpath(root,"bundle-$(mode)"))
            relocated=bundle*"-moved";mv(bundle,relocated)
            rm(joinpath(root,"mode-$(mode)");recursive=true)
            rm(fake)
            loaded=withenv("TEMPO2_CMD"=>"/missing/engine") do;load_noiseless_template(relocated);end
            @test loaded.metadata==t.metadata && loaded.lineage==t.lineage
            @test G._tim_render(loaded.document)==G._tim_render(t.document)
            @test [e.tokens[4] for e in toa_entries(loaded.document)]==[e.tokens[4] for e in toa_entries(src.document)]
            for mutate in [x->(x.metadata["bootstrap"]["anchor_output_index"]=2),
                           x->(x.metadata["bootstrap"]["initial_margin_cycles"]="0.5"),
                           x->(x.metadata["bootstrap"]["initial_pulse_numbers"][1]="999"),
                           x->(x.metadata["bootstrap"]["mode_switch_max_difference_s"]="1e-8"),
                           x->(x.metadata["history"][1]["mode"]=1),
                           x->(x.metadata["history"][1]["corrections"]=8),
                           x->(x.metadata["source_ids"][1]["source_origin"]="file"),
                           x->(x.metadata["epoch_displacements_s"][1]="0")]
                changed=deepcopy(t);mutate(changed)
                @test_throws Exception save_noiseless_template(changed,joinpath(root,"invalid"))
                @test !ispath(joinpath(root,"invalid"))
            end
            changed=deepcopy(t)
            changed.metadata["format_version"]=1
            @test_throws ArgumentError save_noiseless_template(changed,joinpath(root,"invalid"))
            changed=deepcopy(t)
            changed.artifacts["bootstrap.par"]=Vector{UInt8}(codeunits(replace(String(copy(t.artifacts["bootstrap.par"])),"2.5"=>"2.6")))
            changed.metadata["artifact_sha256"]["bootstrap.par"]=G._template_sha(changed.artifacts["bootstrap.par"])
            @test_throws ArgumentError save_noiseless_template(changed,joinpath(root,"invalid"))
        end
        existing=tim_source(joinpath(root,"original.tim"))
        f=joinpath(root,"mode0.tim");write(f,replace(G._tim_render(existing.document),"FORMAT 1"=>"FORMAT 1\nMODE 0"))
        @test_throws ArgumentError create_noiseless_template(tim_source(f),par;workspace=joinpath(root,"no-launch"),tempo=Tempo2(root))
        @test !ispath(joinpath(root,"no-launch"))
    end
end

@testset "Explicit schedules and generated failure paths" begin
    G=GravityToolsNext
    mktempdir() do root
        _,par=template_test_inputs(root)
        src=explicit_schedule(["55825.0000119","55825.0000235","55825.0000349"];
            observatory="pks",frequency_mhz=1400,uncertainty_us=[2,3,4],backend=["A","B","A"])
        fake=template_fake_executable(root,"ok")
        t=withenv("TEMPO2_CMD"=>fake) do
            create_noiseless_template(src,par;workspace=joinpath(root,"explicit"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
        end
        @test t.metadata["source_settings"]["kind"]=="explicit_schedule"
        @test last(t.metadata["history"])["mode"]==1
        @test last(t.metadata["history"])["max_abs_s"]<=1e-9
        @test t.metadata["bootstrap"]["anchor_source_id"]["source_origin"]=="generated"
        @test_throws ArgumentError save_noiseless_template(t,par)
        @test_throws ArgumentError create_noiseless_template(src,par;workspace=joinpath(root,"explicit"),tempo=Tempo2(root))
        h=first(t.metadata["history"]);q=h["query"]
        entries=toa_entries(G._template_read_bytes(t.artifacts[q*"/candidate.tim"],128))
        raw=String(copy(t.artifacts[q*"/tmp/rows.tsv"]))
        @test_throws ArgumentError G._template_residual_table(replace(raw,"T000001"=>"X000001"),entries)
        @test_throws ArgumentError G._template_residual_table(replace(raw,"T000001"=>"T000002"),entries)
        # Change metadata AND its outer checksum: load must check scientific
        # consistency rather than merely trust an internally matching checksum.
        bundle=save_noiseless_template(t,joinpath(root,"corrupt-metadata"))
        record=JLD2.load(joinpath(bundle,"template.jld2"),"record")
        record["metadata"]["bootstrap"]["anchor_epoch_mjd"]="55826.0"
        JLD2.jldsave(joinpath(bundle,"template.jld2");record)
        checks=joinpath(bundle,"checksums.tsv")
        lines=readlines(checks)
        index=only(findall(l->endswith(l,"\ttemplate.jld2"),lines))
        lines[index]=bytes2hex(sha256(read(joinpath(bundle,"template.jld2"))))*"\ttemplate.jld2"
        write(checks,join(lines,'\n')*"\n")
        @test_throws ArgumentError load_noiseless_template(bundle)
    end
    for (mode,status,settings) in [
            # This old ambiguity fixture never converges. It must now exhaust
            # its bounded corrections rather than fail the preliminary margin.
            ("ambiguity",:bootstrap_nonconvergence,IdealizationSettings(timeout_s=10)),
            ("nonconvergence",:bootstrap_nonconvergence,IdealizationSettings(timeout_s=10,max_bootstrap_corrections=1)),
            ("ok",:bootstrap_nonconvergence,IdealizationSettings(timeout_s=10,max_bootstrap_corrections=0)),
            ("parity",:mode_switch_failed,IdealizationSettings(timeout_s=10)),
            ("pulse",:validation_failed,IdealizationSettings(timeout_s=10)),
            ("wrongweight",:validation_failed,IdealizationSettings(timeout_s=10)),
            ("offset",:validation_failed,IdealizationSettings(timeout_s=10)),
            ("nonzero",:engine_failed,IdealizationSettings(timeout_s=10)),
            ("timeout",:timeout,IdealizationSettings(timeout_s=1)),
            ("nonfinite",:validation_failed,IdealizationSettings(timeout_s=10)),
            ("red",:validation_failed,IdealizationSettings(timeout_s=10)),
            ("missing",:validation_failed,IdealizationSettings(timeout_s=10))]
        mktempdir() do root
            _,par=template_test_inputs(root);src=schedule_test_source()
            before=G._tim_render(src.document)
            fake=template_fake_executable(root,mode)
            err=try
                withenv("TEMPO2_CMD"=>fake) do
                    create_noiseless_template(src,par;workspace=joinpath(root,"job"),tempo=Tempo2(root),settings)
                end
                nothing
            catch e;e;end
            @test err isa IdealizationError
            @test err.status==status
            @test isfile(joinpath(root,"job","failure.txt"))
            @test G._tim_render(src.document)==before
            queries=filter(n->occursin(r"^q\d+$",n),readdir(joinpath(root,"job")))
            @test length(queries)<=settings.max_bootstrap_corrections+1+settings.max_corrections+1
            mode=="ambiguity" && @test length(queries)==settings.max_bootstrap_corrections+1
            mode in ("timeout","nonzero") && @test occursin("partial-output",read(joinpath(root,"job",first(queries),"stdout.txt"),String))
        end
    end
end
