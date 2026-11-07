using Test, GravityToolsNext, SHA

# A persistent reservation bounds reruns, including failures and interrupted
# attempts. An unresolved reservation is never reset automatically.
function _reserve_template_calls(budget_dir,cap;limit=24)
    1<=limit<=24 && 1<=cap<=limit || error("Invalid packet reservation")
    isdir(budget_dir) || error("GTN_NOISELESS_BUDGET_DIR must be an existing evidence directory")
    lockdir=joinpath(budget_dir,"engine-budget.lock")
    mkdir(lockdir)
    try
        state=joinpath(budget_dir,"engine-budget.txt")
        ledger=joinpath(budget_dir,"invocations.jsonl")
        used=isfile(state) ? parse(Int,strip(read(state,String))) : (isfile(ledger) ? length(readlines(ledger)) : 0)
        used+cap<=limit || error("Real TEMPO2 budget would exceed $limit calls; do not reset the ledger")
        write(state,string(used+cap))
        return used
    catch
        rm(lockdir)
        rethrow()
    end
end

function _with_template_budget(f,budget,workspace,cap;limit=24)
    used=_reserve_template_calls(budget,cap;limit)
    try
        return f()
    finally
        calls=_record_template_calls(budget,used,workspace)
        println("Closing packet calls: ",used+calls,"/",limit,"; workspace=",workspace)
    end
end

function _record_template_calls(budget_dir,used,workspace)
    # Discover invocation records, independent of stage/job naming. Reservations
    # conservatively charge a recorded attempt even if its status is incomplete.
    jobs=String[]
    if isdir(workspace)
        for (dir,_,names) in walkdir(workspace)
            "invocation.json" in names && push!(jobs,dir)
        end
    end
    sort!(jobs)
    for (i,dir) in enumerate(jobs)
        p=joinpath(dir,"invocation.json")
        invocation=read(p,String)
        count_match=match(r"\"toa_count\":(\d+)",invocation)
        count_match===nothing && error("Invocation must record its TOA count")
        count=parse(Int,only(count_match.captures))
        1<=count<=128 || error("Packet TOA budget violated")
        open(joinpath(budget_dir,"invocations.jsonl"),"a") do io
            println(io,GravityToolsNext._to_json(Dict("number"=>used+i,"stage"=>"template_integration",
                "job"=>dir,"invocation_json"=>invocation,"count"=>count,
                "execution_json"=>isfile(joinpath(dir,"execution.json")) ? read(joinpath(dir,"execution.json"),String) : nothing,
                "engine_sha256"=>bytes2hex(sha256(read(tempo_cmd_path(Tempo2())))))))
        end
    end
    write(joinpath(budget_dir,"engine-budget.txt"),string(used+length(jobs)))
    rm(joinpath(budget_dir,"engine-budget.lock"))
    return length(jobs)
end

@testset "Binary-free operational budget recorder" begin
    mktempdir() do root
        workspace=joinpath(root,"workspace");mkdir(workspace)
        job=joinpath(workspace,"bootstrap-mode0");mkdir(job)
        GravityToolsNext._write_json(joinpath(job,"invocation.json"),Dict("toa_count"=>16,"mode"=>0,"stage"=>"bootstrap"))
        GravityToolsNext._write_json(joinpath(job,"execution.json"),Dict("native_exit"=>0))
        fake=joinpath(root,"fake-binary");write(fake,"not executed")
        used=_reserve_template_calls(root,2)
        @test used==0 && isdir(joinpath(root,"engine-budget.lock"))
        calls=withenv("TEMPO2_CMD"=>fake) do;_record_template_calls(root,used,workspace);end
        @test calls==1
        @test strip(read(joinpath(root,"engine-budget.txt"),String))=="1"
        @test occursin("bootstrap-mode0",read(joinpath(root,"invocations.jsonl"),String))
        @test !ispath(joinpath(root,"engine-budget.lock"))
        @test_throws ErrorException _reserve_template_calls(root,16;limit=16)
        @test_throws ErrorException _reserve_template_calls(root,1;limit=25)
        @test strip(read(joinpath(root,"engine-budget.txt"),String))=="1"
        @test_throws ErrorException _reserve_template_calls(root,24)
        @test strip(read(joinpath(root,"engine-budget.txt"),String))=="1"
        @test !ispath(joinpath(root,"engine-budget.lock"))
    end
end

@testset "Gated noiseless template integration" begin
    enabled=get(ENV,"GTN_ENABLE_NOISELESS_INTEGRATION","0")=="1" && "existing" in split(get(ENV,"GTN_S3_CASES","existing,regular,ambiguous,explicit"),',')
    if !enabled
        @test_skip enabled
    else
        budget=get(ENV,"GTN_NOISELESS_BUDGET_DIR","")
        isdir(budget) || error("Explicit persistent GTN_NOISELESS_BUDGET_DIR is required")
        root=mktempdir(budget;prefix="template-integration-",cleanup=false)
        fixture=joinpath(@__DIR__,"data/ddstg_smoke/j1141-6545_ddstg_gr_pdfb1_white_noise")
        source_all=tim_source(fixture*".tim")
        entries=toa_entries(source_all.document)
        indices=round.(Int,range(1,length(entries),length=24))
        indices[7],indices[8]=indices[8],indices[7]
        push!(indices,indices[12])
        input=joinpath(root,"input.tim")
        write(input,"FORMAT 1\nMODE 1\n"*join((chomp(entries[i].raw) for i in indices),"\n")*"\n")
        source=select_toas(tim_source(input);method=:stratified,count=25,seed=42)
        par=joinpath(root,"input.par")
        # These supported red controls must be removed from the evaluator;
        # a local deterministic selector is preserved, with its -1 fit disabled.
        write(par,read(fixture*".par",String)*"\nJUMP -fe MULT_1 0.0002 -1\nTNRedAmp -14\nTNRedGam 4\nTNRedC 3\nTNRedFLow 0\nTNRedCorner 0\nTNsubtractRed 1\n")
        before_tim=read(input);before_par=read(par)
        workspace=joinpath(root,"job")
        used=_reserve_template_calls(budget,6)
        calls=0
        try
            t=create_noiseless_template(source,par;workspace,
                resource_paths=[joinpath(tempo_data_dir(Tempo2()),"ephemeris/DE430.1950.2050")])
            @test length(t.metadata["history"])<=6
            @test haskey(t.metadata["code_identity"]["git"],"head")
            final=last(t.metadata["history"])
            @test final["max_abs_s"]<=1e-9
            @test final["max_weighted_centered_s"]<=1e-9
            @test all(h->h["pulse_numbers"]==first(t.metadata["history"])["pulse_numbers"],t.metadata["history"])
            @test !occursin("TNRed",String(copy(t.artifacts["reference.par"])))
            @test occursin("JUMP -fe MULT_1",String(copy(t.artifacts["reference.par"])))
            @test all(h->String(copy(t.artifacts[h["query"]*"/tmp/param.labels"]))=="0 param_ZERO 0\n",t.metadata["history"])
            @test read(input)==before_tim && read(par)==before_par
            @test length(unique(e.id for e in toa_entries(t.document)))==25
            save_noiseless_template(t,joinpath(root,"bundle"))
            relocated=joinpath(root,"relocated");mv(joinpath(root,"bundle"),relocated)
            # Rename the old job; loading cannot fall back to original absolute paths.
            mv(workspace,joinpath(root,"retained-job"))
            try
                loaded=withenv("TEMPO2_CMD"=>"/unavailable") do;load_noiseless_template(relocated);end
                @test GravityToolsNext._tim_render(loaded.document)==GravityToolsNext._tim_render(t.document)
                @test loaded.metadata==t.metadata
                @test loaded.lineage==t.lineage
            finally
                mv(joinpath(root,"retained-job"),workspace)
            end
            GravityToolsNext._write_json(joinpath(root,"measured.json"),Dict("history"=>t.metadata["history"],"workspace"=>workspace))
            println("Real noiseless template verified: ",final["max_abs_s"]*1e9," ns; evidence=",root)
        finally
            calls=_record_template_calls(budget,used,workspace)
            println("Real-engine calls in template test: ",calls,"; total ledger usage=",used+calls,"/24")
        end
    end
end

if get(ENV,"GTN_ENABLE_S3_CLOSURE","0")=="1"
    get(ENV,"GTN_ENABLE_NOISELESS_INTEGRATION","0")=="1" || error("Closing matrix requires the existing-TIM regression gate too")
    @testset "Gated generated schedule closing matrix" begin
        G=GravityToolsNext
        budget=ENV["GTN_NOISELESS_BUDGET_DIR"]
        root=mktempdir(budget;prefix="schedule-integration-",cleanup=false)
        fixture=joinpath(@__DIR__,"data/ddstg_smoke/j1141-6545_ddstg_gr_pdfb1_white_noise.par")
        regular=regular_schedule("53590.0";count=16,cadence_s=600,observatory="pks",
            frequency_mhz=1400,uncertainty_us=2,backend="PDFB1",mode=1)
        epochs=["53590.0","53590.013","53590.17","53591.23","53592.41","53593.67",
            "53595.11","53596.37","53597.83","53599.19","53601.53","53603.07",
            "53605.31","53607.79","53609.43","53611.123456"]
        ambiguous=explicit_schedule(epochs;observatory="pks",frequency_mhz=1400,
            uncertainty_us=[2,5,3,7,4,6,3,8,2,5,4,7,3,6,4,5],
            backend=[isodd(i) ? "PDFB1" : "PDFB2" for i in eachindex(epochs)])
        # A second independent literal schedule, fixed before its first query. The
        # formerly margin-screened schedule above retains its original epochs.
        # Neither schedule is tuned using residuals.
        explicit=explicit_schedule(["53595.0","53595.08","53597.17","53598.03","53599.49","53601.11",
            "53602.73","53604.27","53606.41","53609.09","53612.81","53616.37"];
            observatory="pks",frequency_mhz=1400,uncertainty_us=[2,5,3,7,4,6,3,8,2,5,4,7],
            backend=[isodd(i) ? "PDFB1" : "PDFB2" for i in 1:12])
        cases=split(get(ENV,"GTN_S3_CASES","existing,regular,ambiguous,explicit"),',')
        for (name,source) in [("regular",regular),("ambiguous",ambiguous),("explicit",explicit)]
            name in cases || continue
            par=joinpath(root,name*".par")
            write(par,read(fixture,String)*(name!="regular" ? "\nJUMP -be PDFB2 0.0002 -1\nTNEF -be PDFB2 1.3\n" : ""))
            before=G._tim_render(source.document);par_before=read(par)
            workspace=joinpath(root,name)
            t=_with_template_budget(budget,workspace,10) do
                create_noiseless_template(source,par;workspace)
            end
            h=t.metadata["history"];b=t.metadata["bootstrap"]
            @test t.metadata["initial_margin_policy"]=="diagnostic_only"
            name=="ambiguous" && @test parse(BigFloat,b["initial_margin_cycles"])<BigFloat("0.01")
            @test t.metadata["source_origin"]=="generated" && t.lineage==source.lineage
            @test all(x->x["pulse_numbers"]==b["initial_pulse_numbers"],h)
            @test last(h)["mode"]==1 && last(h)["max_abs_s"]<=1e-9 && last(h)["max_weighted_centered_s"]<=1e-9
            @test parse(BigFloat,b["mode_switch_max_difference_s"])<=1e-14
            @test b["entrance_limit_s"]==2e-7
            @test [e.tokens[4] for e in toa_entries(t.document)]==[e.tokens[4] for e in toa_entries(source.document)]
            @test G._tim_render(source.document)==before && read(par)==par_before
            @test length(filter(x->x["mode"]==0,h))<=4 && length(filter(x->x["mode"]==1,h))<=6
            name=="explicit" && @test occursin("JUMP -be PDFB2",String(copy(t.artifacts["reference.par"])))
            bundle=save_noiseless_template(t,joinpath(root,name*"-bundle"))
            relocated=bundle*"-moved";mv(bundle,relocated)
            mv(workspace,workspace*"-retained")
            try
                loaded=withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(relocated);end
                @test loaded.metadata==t.metadata && loaded.lineage==source.lineage
                @test G._tim_render(loaded.document)==G._tim_render(t.document)
                withenv("TEMPO2_CMD"=>"/no/engine") do
                    save_noiseless_template(loaded,joinpath(root,name*"-resaved"))
                end
                G._write_json(joinpath(root,name*"-measured.json"),loaded.metadata)
            finally
                mv(workspace*"-retained",workspace)
            end
            println("Real ",name," schedule verified: ",last(h)["max_abs_s"]*1e9," ns; bundle=",relocated)
            if name=="explicit"
                # One independent MODE 1 query of the final named bundle TIM.
                verify=joinpath(root,"independent");mkdir(verify)
                write(joinpath(verify,"reference.par"),t.artifacts["reference.par"])
                cp(joinpath(relocated,"final.tim"),joinpath(verify,"final.tim"))
                q=_with_template_budget(budget,verify,1) do
                    G._template_query(verify,joinpath(verify,"final.tim"),String(copy(t.artifacts["reference.par"])),
                        IdealizationSettings(),Tempo2(),1,b["initial_pulse_numbers"];stage="independent_verification")
                end
                @test q.table.max_abs_s<=1e-9 && q.table.max_weighted_centered_s<=1e-9
                @test read(q.prepared.tim_path)==read(q.job.tim_in_path)==read(joinpath(relocated,"final.tim"))
                G._write_json(joinpath(root,"independent-measured.json"),Dict("max_abs_s"=>q.table.max_abs_s,
                    "max_weighted_centered_s"=>q.table.max_weighted_centered_s,"pulse_numbers"=>q.table.pulses,
                    "final_tim_sha256"=>bytes2hex(sha256(read(q.prepared.tim_path)))))
            end
        end
    end
end
