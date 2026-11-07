using Test, GravityToolsNext, SHA

const TEMPLATE_TEST_PAR="""PSRJ TEST
BINARY DDSTG
MODEL_MODE GR
F0 2.5 -1
F1 -1e-14 1
PEPOCH 55825
START 55824
FINISH 55826
JUMP -be A 0.0002 1
REFPHS TZR
MODE 1
TNEF -be A 1.2
TNEQ -be A -6
"""

function template_test_inputs(root)
    par=joinpath(root,"original.par");write(par,TEMPLATE_TEST_PAR)
    doc=regular_schedule("55825.0";count=3,cadence_s=1,observatory="pks",frequency_mhz=1400,
        uncertainty_us=[2,3,4],backend="A").document
    for (e,d) in zip(toa_entries(doc),[1e-4,-8e-5,-2e-5]);shift_toa_epoch!(doc,e.id,d);end
    tim=joinpath(root,"original.tim");write_tim_document(doc,tim)
    return tim_source(tim),par
end

function template_fake_executable(root,mode;marker=nothing)
    shellquote(s)="'"*replace(s,"'"=>"'\\''")*"'"
    path=joinpath(root,"fake-$(mode).sh")
    julia=Base.julia_cmd().exec[1]
    helper=joinpath(@__DIR__,"helpers","template_engine.jl")
    home_line=haskey(ENV,"HOME") ? "export HOME="*shellquote(ENV["HOME"])*"\n" : ""
    marker_line=marker===nothing ? "" : "printf started >> "*shellquote(marker)*"\n"
    write(path,"#!/bin/sh\n"*marker_line*home_line*"export JULIA_DEPOT_PATH="*shellquote(first(DEPOT_PATH))*"\nexport JULIA_LOAD_PATH=@stdlib\nexec "*shellquote(julia)*" --startup-file=no "*shellquote(helper)*" "*shellquote(mode)*" \"\$@\"\n")
    chmod(path,0o755)
    return path
end

@testset "Noiseless template preparation" begin
    G=GravityToolsNext
    @test_throws ArgumentError IdealizationSettings(max_corrections=6)
    @test_throws ArgumentError IdealizationSettings(timeout_s=0)
    @test_throws ArgumentError IdealizationSettings(max_step_phase=NaN)
    @test_throws ArgumentError IdealizationSettings(tolerance_s=Inf)
    @test_throws ArgumentError IdealizationSettings(precision_bits=80)
    model=G._prepare_template_par(TEMPLATE_TEST_PAR)
    @test occursin("MEAN",model.text)
    @test !occursin("TNEF",model.text) && !occursin("TNEQ",model.text)
    ps=G._tim_precision(256) do;G._template_read_par(model.text);end
    @test all(p->p.flag in (nothing,0),ps)
    @test get_param(ps,Symbol("JUMP -be A")).value==BigFloat("0.0002")
    @test get_param(ps,:F1).value==G._tim_precision(256) do;BigFloat("-1e-14");end
    red=G._prepare_template_par(TEMPLATE_TEST_PAR*"TNRedAmp -14\nTNRedGam 4\nTNRedC 3\nTNRedFLow 0\nTNRedCorner 0\nTNsubtractRed 1\n")
    @test !occursin("TNRed",red.text)
    unknown=G._prepare_template_par(TEMPLATE_TEST_PAR*"DDSTG_FUTURE 3.25 1\n")
    @test occursin("DDSTG_FUTURE",unknown.text)
    for line in ["TNFuture 1","TNDMAmp -14","TN_QpPeriod 1","ECORR -be A 1","WAVE_OM 1","SIFUNC 1","RNAMP 1","CONSTRAIN F0","TRACK -2","START 55824 1","JUMP -gtnid T000001 0.1 0","FDJUMP -gtnid T000001 0.1 0"]
        @test_throws Exception G._prepare_template_par(TEMPLATE_TEST_PAR*line*"\n")
    end
    @test_throws ArgumentError G._prepare_template_par(TEMPLATE_TEST_PAR*"f0 3\n")
    @test_throws ArgumentError G._prepare_template_par(TEMPLATE_TEST_PAR*"F3 1e-20 0 extra trailing\n")
    mktempdir() do root
        source,par=template_test_inputs(root)
        generated=regular_schedule("55825.0";count=3,cadence_s=1,observatory="pks",frequency_mhz=1400,uncertainty_us=2)
        @test length(G._validate_template_source(generated,IdealizationSettings(),model.selectors))==3
        disguised=tim_source(generated.document)
        @test_throws ArgumentError create_noiseless_template(disguised,par;workspace=joinpath(root,"generated"),tempo=Tempo2(root))
        @test !ispath(joinpath(root,"generated"))
        for flag in ["-to","-mystery","-gtnid"]
            text=replace(GravityToolsNext._tim_render(source.document),"-be A"=>"-be A $flag 1")
            file=joinpath(root,replace(flag,"-"=>"")*".tim");write(file,text)
            @test_throws ArgumentError create_noiseless_template(tim_source(file),par;workspace=joinpath(root,"bad"),tempo=Tempo2(root))
        end
    end
end

@testset "Noiseless templates with fake engine" begin
    mktempdir() do root
        source,par=template_test_inputs(root)
        before=GravityToolsNext._tim_render(source.document);par_before=read(par)
        fake=template_fake_executable(root,"ok")
        template=withenv("TEMPO2_CMD"=>fake) do
            create_noiseless_template(source,par;workspace=joinpath(root,"job"),tempo=Tempo2(root),settings=IdealizationSettings(timeout_s=10))
        end
        @test template isa NoiselessTemplate
        @test length(template.metadata["history"])==2
        @test last(template.metadata["history"])["max_abs_s"]<1e-9
        @test GravityToolsNext._tim_render(source.document)==before && read(par)==par_before
        @test all(e->last(e.flags[end]) in ["T000001","T000002","T000003"],toa_entries(template.document))
        @test template.lineage==source.lineage
        @test !template.metadata["resource_provenance_complete"]
        @test haskey(template.metadata,"code_identity")
        @test haskey(template.metadata["code_identity"]["git"],"head")
        target=save_noiseless_template(template,joinpath(root,"bundle"))
        @test_throws ArgumentError save_noiseless_template(template,target)
        @test_throws ArgumentError save_noiseless_template(template,par)
        link=joinpath(root,"alias");symlink(target,link)
        @test_throws ArgumentError save_noiseless_template(template,link)
        @test_throws ArgumentError load_noiseless_template(link)
        relocated=joinpath(root,"relocated");mv(target,relocated)
        # Delete all old job artifacts and disable the executable to prove reuse.
        rm(joinpath(root,"job");recursive=true);rm(fake)
        loaded=withenv("TEMPO2_CMD"=>"/no/engine") do;load_noiseless_template(relocated);end
        @test GravityToolsNext._tim_render(loaded.document)==GravityToolsNext._tim_render(template.document)
        @test loaded.metadata==template.metadata
        @test loaded.lineage==template.lineage
        @test [e.tokens[3] for e in toa_entries(loaded.document)]==template.metadata["final_epochs_mjd"]
        save_noiseless_template(loaded,joinpath(root,"resaved"))
        # A changed verified document or saved artifact must be rejected.
        changed=deepcopy(template);shift_toa_epoch!(changed.document,toa_entries(changed.document)[1].id,1e-9)
        @test_throws ArgumentError save_noiseless_template(changed,joinpath(root,"changed"))
        @test !ispath(joinpath(root,"changed"))
        nested=joinpath(root,"resaved","linked");symlink(par,nested)
        @test_throws ArgumentError load_noiseless_template(joinpath(root,"resaved"))
        rm(nested)
        invalid=joinpath(root,"resaved","checksums.tsv")
        saved=read(invalid,String)
        write(invalid,saved*repeat("0",64)*"\t../escape\n")
        @test_throws ArgumentError load_noiseless_template(joinpath(root,"resaved"))
        write(invalid,saved)
        write(joinpath(relocated,"final.tim"),"corrupt")
        @test_throws ArgumentError load_noiseless_template(relocated)
        extra=joinpath(root,"resaved","extra");write(extra,"extra")
        @test_throws ArgumentError load_noiseless_template(joinpath(root,"resaved"))
        rm(extra)
        rm(joinpath(root,"resaved","checksums.tsv"))
        @test_throws ArgumentError load_noiseless_template(joinpath(root,"resaved"))
    end
    mktempdir() do root
        source,par=template_test_inputs(root)
        fake=template_fake_executable(root,"phase_guard")
        err=try
            withenv("TEMPO2_CMD"=>fake) do
                create_noiseless_template(source,par;workspace=joinpath(root,"loose-tolerance"),tempo=Tempo2(root),
                    settings=IdealizationSettings(tolerance_s=1.0,timeout_s=10))
            end
            nothing
        catch e;e;end
        @test err isa IdealizationError && err.status==:phase_guard
    end
    for mode in ["timeout","nonzero","missing","duplicate","deleted","red","nonfinite","sat","error","fit","constraint","init","pulse","phase_guard","nonconvergence"]
        mktempdir() do root
            source,par=template_test_inputs(root)
            fake=template_fake_executable(root,mode)
            err=try
                withenv("TEMPO2_CMD"=>fake) do
                    create_noiseless_template(source,par;workspace=joinpath(root,"failure"),tempo=Tempo2(root),
                        settings=IdealizationSettings(max_corrections=1,timeout_s=mode=="timeout" ? 1 : 10))
                end
                nothing
            catch e;e;end
            @test err isa IdealizationError
            @test isfile(joinpath(root,"failure","failure.txt"))
            @test err.status==(mode=="timeout" ? :timeout : mode=="nonzero" ? :engine_failed : mode=="phase_guard" ? :phase_guard : mode=="nonconvergence" ? :nonconvergence : :validation_failed)
        end
    end
end
