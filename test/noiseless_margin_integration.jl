using Test, GravityToolsNext, SHA

# Opt-in, separately budgeted acceptance of the formerly margin-screened case.
# No real calls or extra skipped testsets occur in the default suite.
if get(ENV,"GTN_ENABLE_MARGIN_INTEGRATION","0")=="1"
    @testset "Diagnostic margin bounded real integration" begin
        G=GravityToolsNext
        budget=ENV["GTN_NOISELESS_BUDGET_DIR"]
        root=mktempdir(budget;prefix="margin-integration-",cleanup=false)
        fixture=joinpath(@__DIR__,"data/ddstg_smoke/j1141-6545_ddstg_gr_pdfb1_white_noise")
        epochs=["53590.0","53590.013","53590.17","53591.23","53592.41","53593.67",
            "53595.11","53596.37","53597.83","53599.19","53601.53","53603.07",
            "53605.31","53607.79","53609.43","53611.123456"]
        source=explicit_schedule(epochs;observatory="pks",frequency_mhz=1400,
            uncertainty_us=[2,5,3,7,4,6,3,8,2,5,4,7,3,6,4,5],
            backend=[isodd(i) ? "PDFB1" : "PDFB2" for i in eachindex(epochs)])
        par=joinpath(root,"generated.par")
        write(par,read(fixture*".par",String)*"\nJUMP -be PDFB2 0.0002 -1\nTNEF -be PDFB2 1.3\n")
        source_before=G._tim_render(source.document);par_before=read(par)
        workspace=joinpath(root,"generated")
        t=_with_template_budget(budget,workspace,10;limit=16) do
            create_noiseless_template(source,par;workspace)
        end
        m=t.metadata;b=m["bootstrap"];h=m["history"]
        @test m["format_version"]==3 && m["initial_margin_policy"]=="diagnostic_only"
        @test 0<parse(BigFloat,b["initial_margin_cycles"])<BigFloat("0.01")
        @test m["source_origin"]=="generated" && t.lineage==source.lineage
        @test all(x->x["pulse_numbers"]==b["initial_pulse_numbers"],h)
        @test parse(BigFloat,b["mode_switch_max_difference_s"])<=1e-14
        @test b["entrance_limit_s"]==2e-7
        @test last(h)["mode"]==1 && last(h)["max_abs_s"]<=1e-9 && last(h)["max_weighted_centered_s"]<=1e-9
        @test [e.tokens[4] for e in toa_entries(t.document)]==[e.tokens[4] for e in toa_entries(source.document)]
        @test G._tim_render(source.document)==source_before && read(par)==par_before
        @test occursin("JUMP -be PDFB2",String(copy(t.artifacts["reference.par"])))
        bundle=save_noiseless_template(t,joinpath(root,"bundle"));moved=bundle*"-moved";mv(bundle,moved)
        mv(workspace,workspace*"-retained")
        loaded=try
            withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(moved);end
        finally
            mv(workspace*"-retained",workspace)
        end
        @test loaded.metadata==m && loaded.lineage==source.lineage
        @test G._tim_render(loaded.document)==G._tim_render(t.document)
        withenv("TEMPO2_CMD"=>"/no/engine") do;save_noiseless_template(loaded,joinpath(root,"resaved"));end
        G._write_json(joinpath(root,"generated-measured.json"),m)
        verify=joinpath(root,"independent");mkdir(verify)
        write(joinpath(verify,"reference.par"),t.artifacts["reference.par"])
        cp(joinpath(moved,"final.tim"),joinpath(verify,"final.tim"))
        q=_with_template_budget(budget,verify,1;limit=16) do
            G._template_query(verify,joinpath(verify,"final.tim"),String(copy(t.artifacts["reference.par"])),
                IdealizationSettings(),Tempo2(),1,b["initial_pulse_numbers"];stage="independent_verification")
        end
        @test q.table.max_abs_s<=1e-9 && q.table.max_weighted_centered_s<=1e-9
        @test q.table.pulses==b["initial_pulse_numbers"]
        @test read(q.prepared.tim_path)==read(q.job.tim_in_path)==read(joinpath(moved,"final.tim"))
        G._write_json(joinpath(root,"independent-measured.json"),Dict("max_abs_s"=>q.table.max_abs_s,
            "max_weighted_centered_s"=>q.table.max_weighted_centered_s,"pulse_numbers"=>q.table.pulses,
            "final_tim_sha256"=>bytes2hex(sha256(read(q.prepared.tim_path)))))

        entries=toa_entries(tim_source(fixture*".tim").document)
        indices=round.(Int,range(1,length(entries),length=24))
        indices[7],indices[8]=indices[8],indices[7];push!(indices,indices[12])
        input=joinpath(root,"existing.tim")
        write(input,"FORMAT 1\nMODE 1\n"*join((chomp(entries[i].raw) for i in indices),"\n")*"\n")
        existing=select_toas(tim_source(input);method=:stratified,count=25,seed=42)
        par=joinpath(root,"existing.par")
        write(par,read(fixture*".par",String)*"\nJUMP -fe MULT_1 0.0002 -1\nTNRedAmp -14\nTNRedGam 4\nTNRedC 3\nTNRedFLow 0\nTNRedCorner 0\nTNsubtractRed 1\n")
        input_before=read(input);par_before=read(par);workspace=joinpath(root,"existing")
        t=_with_template_budget(budget,workspace,5;limit=16) do
            create_noiseless_template(existing,par;workspace,settings=IdealizationSettings(max_corrections=4))
        end
        m=t.metadata;h=m["history"]
        @test m["format_version"]==3 && m["initial_margin_policy"]=="diagnostic_only"
        @test m["bootstrap"]===nothing && all(x->x["mode"]==1,h)
        @test all(x->x["pulse_numbers"]==first(h)["pulse_numbers"],h)
        @test last(h)["max_abs_s"]<=1e-9 && last(h)["max_weighted_centered_s"]<=1e-9
        @test all(x->String(copy(t.artifacts[x["query"]*"/tmp/param.labels"]))=="0 param_ZERO 0\n",h)
        @test !occursin("TNRed",String(copy(t.artifacts["reference.par"])))
        @test occursin("JUMP -fe MULT_1",String(copy(t.artifacts["reference.par"])))
        @test read(input)==input_before && read(par)==par_before && t.lineage==existing.lineage
        saved=save_noiseless_template(t,joinpath(root,"existing-bundle"))
        reloaded=withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(saved);end
        @test reloaded.metadata==m && reloaded.lineage==t.lineage
        G._write_json(joinpath(root,"existing-measured.json"),m)
        println("Margin integration evidence: ",root)
    end
end
