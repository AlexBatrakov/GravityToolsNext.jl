using Test, GravityToolsNext, JLD2, SHA

function margin_test_source()
    GravityToolsNext._tim_precision(128) do
        anchor=BigFloat("55825.0")
        epochs=[anchor+x/86400 for x in [BigFloat(0),BigFloat("1.199"),BigFloat("2.07")]]
        explicit_schedule(epochs;observatory="pks",frequency_mhz=1400,
            uncertainty_us=[2,3,4],backend=["A","B","A"])
    end
end

function margin_test_rewrite_record(bundle,mutate)
    file=joinpath(bundle,"template.jld2")
    record=JLD2.load(file,"record");mutate(record);JLD2.jldsave(file;record)
    checks=joinpath(bundle,"checksums.tsv");lines=readlines(checks)
    i=only(findall(l->endswith(l,"\ttemplate.jld2"),lines))
    lines[i]=bytes2hex(sha256(read(file)))*"\ttemplate.jld2"
    write(checks,join(lines,'\n')*"\n")
end

@testset "Diagnostic-only margin in the current format" begin
    G=GravityToolsNext
    mktempdir() do root
        _,par=template_test_inputs(root);source=margin_test_source()
        before=G._tim_render(source.document);fake=template_fake_executable(root,"ok")
        t=withenv("TEMPO2_CMD"=>fake) do
            create_noiseless_template(source,par;workspace=joinpath(root,"job"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
        end
        m=t.metadata;b=m["bootstrap"]
        @test m["format_version"]==3 && m["initial_margin_policy"]=="diagnostic_only"
        @test Set(keys(m["settings"]))==Set(string.(fieldnames(IdealizationSettings)))
        @test 0<parse(BigFloat,b["initial_margin_cycles"])<BigFloat("0.01")
        @test all(h->h["pulse_numbers"]==b["initial_pulse_numbers"],m["history"])
        @test [h["mode"] for h in m["history"]]==[0,0,1]
        @test b["mode_switch_max_difference_s"]=="0.0"
        q=m["history"][2]["query"]
        @test first(toa_entries(G._template_read_bytes(t.artifacts[q*"/candidate.tim"],128))).epoch_mjd==first(toa_entries(source.document)).epoch_mjd
        @test last(m["history"])["max_abs_s"]<=1e-9 && last(m["history"])["max_weighted_centered_s"]<=1e-9
        @test G._tim_render(source.document)==before && read(par,String)==TEMPLATE_TEST_PAR
        target=save_noiseless_template(t,joinpath(root,"bundle"));moved=target*"-moved";mv(target,moved)
        rm(joinpath(root,"job");recursive=true);rm(fake)
        loaded=withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(moved);end
        @test loaded.metadata==m && loaded.lineage==source.lineage
        @test G._tim_render(loaded.document)==G._tim_render(t.document)
        withenv("TEMPO2_CMD"=>"/no/engine") do;save_noiseless_template(loaded,joinpath(root,"resaved"));end
        @test JLD2.load(joinpath(root,"resaved/template.jld2"),"record")["version"]==3

        for mutate in [x->delete!(x.metadata,"initial_margin_policy"),
                       x->(x.metadata["initial_margin_policy"]="enforced"),
                       x->(x.metadata["initial_margin_policy"]=false),
                       x->(x.metadata["settings"]["unknown_setting"]=1),
                       x->(x.metadata["bootstrap"]["initial_margin_cycles"]="NaN"),
                       x->delete!(x.metadata["bootstrap"],"initial_margin_cycles"),
                       x->delete!(x.metadata["bootstrap"],"initial_pulse_numbers"),
                       x->(x.metadata["format_version"]=false)]
            invalid=deepcopy(t);mutate(invalid)
            @test_throws Exception save_noiseless_template(invalid,joinpath(root,"invalid"))
            @test !ispath(joinpath(root,"invalid"))
        end
        for version in (1,2,4,true,3.0)
            invalid=deepcopy(t);invalid.metadata["format_version"]=version
            @test_throws ArgumentError save_noiseless_template(invalid,joinpath(root,"unsupported"))
        end
        # Recompute the outer checksum so malformed format/shape records reach
        # the validator rather than merely failing their file hash.
        for (i,(mutate,message)) in enumerate([
                (r->delete!(r["metadata"],"initial_margin_policy"),"initial_margin_policy"),
                (r->(r["metadata"]["initial_margin_policy"]="unknown"),"Unsupported initial margin policy"),
                (r->(r["version"]=1),"Unsupported bundle format"),
                (r->(r["version"]=2),"Unsupported bundle format"),
                (r->(r["version"]=4),"Unsupported bundle format"),
                (r->(r["version"]=true),"Unsupported bundle format"),
                (r->(r["version"]=3.0),"Unsupported bundle format"),
                (r->(r["metadata"]["format_version"]=2),"record/metadata versions differ"),
                (r->delete!(r["metadata"]["settings"],"max_bootstrap_corrections"),"max_bootstrap_corrections"),
                (r->delete!(first(r["metadata"]["history"]),"executed_input_tim"),"executed_input_tim"),
                (r->delete!(first(r["metadata"]["history"]),"mode"),"missing required field 'mode'"),
                (r->delete!(r["metadata"],"source_ids"),"source_ids"),
                (r->delete!(r,"version"),"missing required field 'version'"),
                (r->delete!(r,"lineage"),"missing required field 'lineage'")])
            bundle=save_noiseless_template(t,joinpath(root,"corrupt-$i"))
            margin_test_rewrite_record(bundle,mutate)
            err=try
                withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(bundle);end
                nothing
            catch e;e;end
            @test err isa ArgumentError
            @test occursin(message,sprint(showerror,err))
        end
    end
    mktempdir() do root
        _,par=template_test_inputs(root);fake=template_fake_executable(root,"pulse")
        err=try
            withenv("TEMPO2_CMD"=>fake) do
                create_noiseless_template(margin_test_source(),par;workspace=joinpath(root,"job"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
            end
            nothing
        catch e;e;end
        @test err isa IdealizationError && err.status==:validation_failed
        @test occursin("Pulse assignment changed",read(joinpath(root,"job/failure.txt"),String))
        @test !isdir(joinpath(root,"job/q003"))
    end
end
