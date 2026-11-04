using Test
using GravityToolsNext
using Random

_src_big(x,bits=256)=BigFloat(x,RoundNearest;precision=bits)
_src_ids(s)=[e.id for e in toa_entries(s.document)]
_src_ord(s)=[l.source_ordinal for l in s.lineage]
function _src_fixture(dir; flags=["-be A","-be A","-be B","-be B","-be A","-be B"])
    epochs=["55825.3","55825.1","55825.2","55825.2","55825.4","55825.0"]
    text="C preamble\r\nFORMAT 1\r\nMODE 1\r\n# leading\r\n"
    for i in eachindex(epochs)
        text*="obs\t 1400.00 $(epochs[i])\t2.50 pks $(flags[i]) -x one -x two\r\n"
        i==3 && (text*="# middle\r\n\t\r\n")
    end
    text*="# trailing without final newline"
    path=joinpath(dir,"source.tim");write(path,text)
    return tim_source(path),text
end

@testset "Observation sources and selection" begin
    @testset "TIM composition and current epochs" begin
        mktempdir() do dir
            source,text=_src_fixture(dir)
            doc=source.document
            entries=toa_entries(doc)
            shift_toa_epoch!(doc,entries[2].id,1e-9)
            current=toa_entries(doc)[2].epoch_mjd
            selected=select_toas(source;method=:stride,stride=2,offset=2)
            @test _src_ord(selected)==[2,4,6]
            @test _src_ids(selected)==_src_ids(source)[[2,4,6]]
            @test [x.output_index for x in selected.lineage]==[1,2,3]
            @test toa_entries(selected.document)[1].epoch_mjd==current
            @test toa_entries(selected.document)[1].tokens[3]=="55825.1"
            @test toa_entries(selected.document)[1].uncertainty_us==2.5
            @test toa_entries(selected.document)[1].flags==entries[2].flags
            @test selected.settings.kind==:tim
            @test selected.settings.source_sha256==doc.source_sha256
            set_toa_epoch!(selected.document,entries[2].id,_src_big("55826.0"))
            @test toa_entries(doc)[2].epoch_mjd==current
            output=write_tim_document(selected.document,joinpath(dir,"subset.tim"))
            rendered=read(output,String)
            @test occursin("C preamble\r\nFORMAT 1\r\nMODE 1\r\n# leading\r\n",rendered)
            @test occursin("# middle\r\n\t\r\n",rendered)
            @test endswith(rendered,"# trailing without final newline")
            @test length(toa_entries(read_tim_document(output)))==3
            @test read(doc.source_path,String)==text
            @test_throws ArgumentError write_tim_document(selected.document,doc.source_path;overwrite=true)
            alias=joinpath(dir,"source-alias.tim");hardlink(doc.source_path,alias)
            @test_throws ArgumentError write_tim_document(selected.document,alias;overwrite=true)
            @test all(e.id.origin==:file for e in entries)
            # Direct composition takes IDs as a set, retaining source order.
            subdoc=subset_tim_document(doc,[entries[6].id,entries[2].id])
            @test [e.ordinal for e in toa_entries(subdoc)]==[2,6]
            @test_throws ArgumentError subset_tim_document(doc,[entries[1].id,entries[1].id])
            @test_throws ArgumentError subset_tim_document(doc,[TimTOAID("foreign",1)])
            @test isempty(toa_entries(subset_tim_document(doc,TimTOAID[])))
            @test_throws ArgumentError tim_source(subset_tim_document(doc,TimTOAID[]))
            owned=tim_source(doc)
            shift_toa_epoch!(owned.document,entries[1].id,2e-9)
            @test toa_entries(doc)[1].epoch_mjd==entries[1].epoch_mjd
        end
    end

    @testset "Counts, order, windows and lineage chains" begin
        mktempdir() do dir
            source,_=_src_fixture(dir)
            @test _src_ord(select_toas(source))==collect(1:6)
            @test _src_ord(select_toas(source;method=:stride,stride=2))==[1,3,5]
            @test _src_ord(select_toas(source;method=:stride,stride=20,offset=3))==[3]
            @test _src_ord(select_toas(source;method=:time_window,time_start="55825.1",time_finish="55825.3"))==[1,2,3,4]
            @test _src_ord(select_toas(source;method=:time_window,time_finish="55825.1"))==[2,6]
            @test _src_ord(select_toas(source;method=:time_window,time_start="55825.4"))==[5]
            entries=toa_entries(source.document)
            shift_toa_epoch!(source.document,entries[3].id,1e-9)
            # Uses current MJD, not old tokens: only the unedited duplicate is at .2.
            @test _src_ord(select_toas(source;method=:time_window,time_start="55825.2",time_finish="55825.2"))==[4]
            random=select_toas(source;method=:random,fraction=0.5,seed=4)
            @test length(random.lineage)==3
            @test length(select_toas(source;method=:random,fraction=0.75,seed=4).lineage)==4
            @test issorted(_src_ord(random))
            @test length(unique(_src_ids(random)))==3
            nested=select_toas(select_toas(source;method=:stride,stride=2);method=:stride,stride=2)
            @test _src_ord(nested)==[1,5]
            @test _src_ids(nested)==_src_ids(source)[[1,5]]
            @test length(nested.settings.selections)==2
            @test isempty(source.settings.selections)
            @test last(nested.settings.selections).selected_input_indices==[1,3]
        end
    end

    @testset "Stratification and explicit RNG ownership" begin
        mktempdir() do dir
            source,_=_src_fixture(dir)
            a=select_toas(source;method=:random,count=3,rng=Xoshiro(97))
            b=select_toas(source;method=:random,count=3,rng=Xoshiro(97))
            @test _src_ids(a)==_src_ids(b)
            state=copy(last(a.settings.selections).rng.state_before)
            @test _src_ids(select_toas(source;method=:random,count=3,rng=state))==_src_ids(a)
            @test last(a.settings.selections).rng.seed===nothing
            @test _src_ids(select_toas(source;method=:random,count=3,rng=MersenneTwister(97)))==
                _src_ids(select_toas(source;method=:random,count=3,rng=MersenneTwister(97)))
            seeded=select_toas(source;method=:random,count=3,seed=97)
            @test _src_ids(seeded)==_src_ids(a)
            @test last(seeded.settings.selections).rng.seed==97
            rng=Xoshiro(11);before=copy(rng)
            select_toas(source;method=:random,count=2,rng=rng)
            @test rand(copy(rng),UInt64)!=rand(copy(before),UInt64)
            # Invalid requests are validated before any RNG draw.
            rng=Xoshiro(11);before=copy(rng)
            @test_throws ArgumentError select_toas(source;method=:random,count=100,rng=rng)
            @test rand(copy(rng),UInt64)==rand(copy(before),UInt64)
            strat=select_toas(source;method=:stratified,count=3,seed=42)
            @test last(strat.settings.selections).group_counts==Dict("A"=>2,"B"=>1) # lexical tie
            @test issorted(_src_ord(strat))
            @test _src_ids(select_toas(source;method=:stratified,count=3,rng=Xoshiro(42)))==_src_ids(strat)
            explicit=select_toas(source;method=:stratified,group_counts=Dict("A"=>0,"B"=>2),seed=5)
            @test all(i in [3,4,6] for i in _src_ord(explicit))
            @test length(explicit.lineage)==2
            @test _src_ids(select_toas(source;method=:random,count=6,seed=3))==_src_ids(source)
            @test _src_ids(select_toas(source;method=:stratified,count=6,seed=3))==_src_ids(source)
            # Group extraction rejects every missing/duplicate occurrence, even equal.
            for flag in ("-x A","-be A -be A","-be A -be B")
                broken,_=_src_fixture(dir;flags=fill(flag,6))
                @test_throws ArgumentError select_toas(broken;method=:stratified,count=2,seed=1)
            end
            # No use of the ambient RNG, and no coupling to another owned stream.
            value=fetch(@async begin
                Random.seed!(123);expected=rand(UInt64)
                Random.seed!(123)
                select_toas(source;method=:random,count=2,seed=4)
                explicit_schedule(["55825.0"];observatory="pks",frequency_mhz=1,uncertainty_us=1)
                (expected,rand(UInt64))
            end)
            @test value[1]==value[2]
            independent=Xoshiro(20);reference=copy(independent)
            select_toas(source;method=:random,count=3,rng=Xoshiro(19))
            @test rand(independent)==rand(reference)
            imbalance=explicit_schedule(fill("55825.0",10);observatory="pks",frequency_mhz=1400,
                uncertainty_us=1,backend=["C","C","C","C","C","B","B","B","A","A"])
            @test last(select_toas(imbalance;method=:stratified,count=4,seed=1).settings.selections).group_counts==
                Dict("A"=>1,"B"=>1,"C"=>2)
        end
    end

    @testset "Impossible and ambiguous selections" begin
        mktempdir() do dir
            source,_=_src_fixture(dir)
            requests=[(;method=:bogus),(;method=:random,count=2),(;method=:random,count=0,seed=1),
                (;method=:random,count=-1,seed=1),(;method=:random,count=7,seed=1),
                (;method=:random,count=1,fraction=0.5,seed=1),(;method=:random,fraction=0.1,seed=1),
                (;method=:random,fraction=NaN,seed=1),(;method=:random,fraction=1.1,seed=1),
                (;method=:random,count=2,rng=Xoshiro(1),seed=1),(;method=:random,count=2,seed=-1),
                (;method=:all,seed=1),(;method=:all,count=1),(;method=:all,stride=2),
                (;method=:stride,stride=0),(;method=:stride,stride=2,offset=7),
                (;method=:stride,stride=-1),(;method=:time_window),
                (;method=:time_window,time_start="55826.0"),(;method=:time_window,time_start="55825.4",time_finish="55825.0"),
                (;method=:time_window,time_start=55825.0),
                (;method=:stratified,count=2,group_counts=Dict("A"=>1,"B"=>1),seed=1),
                (;method=:stratified,group_counts=Dict("A"=>1),seed=1),
                (;method=:stratified,group_counts=Dict("A"=>1,"B"=>1,"C"=>0),seed=1),
                (;method=:stratified,group_counts=Dict("A"=>4,"B"=>0),seed=1),
                (;method=:stratified,group_counts=Dict("A"=>0,"B"=>0),seed=1)]
            for kw in requests; @test_throws ArgumentError select_toas(source;kw...); end
        end
    end

    @testset "Explicit schedules, systems, generated origin" begin
        epochs=["55825.2","55825.1","55825.1"]
        flags=["-x"=>"one","-x"=>"two"]
        source=explicit_schedule(epochs;observatory=["pks","@","pks"],frequency_mhz=[1400,0,700],
            uncertainty_us=[1,2,3],backend=["A","B","A"],flags=flags,filenames="same",mode=1)
        entries=toa_entries(source.document)
        @test [e.epoch_mjd for e in entries]==_src_big.(epochs,128)
        @test source.document.origin==:generated
        @test source.document.source_path===nothing
        @test source.document.source_realpath===nothing
        @test source.document.source_file_id===nothing
        @test all(e.id.origin==:generated for e in entries)
        @test length(unique(_src_ids(source)))==3
        @test entries[1].filename==entries[2].filename
        @test entries[1].flags==[flags; "-be"=>"A"]
        @test source.settings.observations[2].uncertainty_us=="2"
        flags[1]="-x"=>"mutated"
        @test toa_entries(source.document)[1].flags[1]==("-x"=>"one")
        @test _src_ord(select_toas(source;method=:time_window,time_finish="55825.1"))==[2,3]
        systems=Dict("A"=>(observatory="pks",frequency_mhz=1400,uncertainty_us=2,flags=["-fe"=>"one"]),
            "B"=>(observatory="@",frequency_mhz=0,uncertainty_us=3))
        generated=explicit_schedule(epochs;backend=["B","A","B"],systems=systems,flags=[["-n"=>"1"],Pair{String,String}[],["-n"=>"3"]])
        @test [e.freq_mhz for e in toa_entries(generated.document)]==[0,1400,0]
        @test toa_entries(generated.document)[2].flags==["-fe"=>"one","-be"=>"A"]
        @test generated.settings.systems==systems
        mktempdir() do dir
            path=write_tim_document(source.document,joinpath(dir,"generated.tim"))
            readback=read_tim_document(path)
            @test readback.origin==:file
            @test only(unique(e.id.origin for e in toa_entries(readback)))==:file
            @test _src_ids(tim_source(readback))!=_src_ids(source)
            @test all(a.epoch_mjd==b.epoch_mjd for (a,b) in zip(entries,toa_entries(readback)))
            # Source IDs remain separately available after writing.
            @test source.lineage[2].source_id==entries[2].id
        end
        @test all(first(p)!="-pn" for e in entries for p in e.flags)
    end

    @testset "Explicit precision and regular cadence" begin
        for bits in (128,256)
            setprecision(BigFloat,80) do
                setrounding(BigFloat,RoundDown) do
                    for start in ("10000.000000001","55825.123456789012345678901234","99999.99999999")
                        source=regular_schedule(start;count=10,cadence_s=_src_big("1e-9"),precision_bits=bits,
                            observatory="pks",frequency_mhz=1400,uncertainty_us=1,backend="A")
                        entries=toa_entries(source.document)
                        @test length(entries)==10
                        @test all(precision(e.epoch_mjd)==bits for e in entries)
                        for (i,e) in enumerate(entries)
                            err=setprecision(BigFloat,256) do
                                setrounding(BigFloat,RoundNearest) do
                                    reference=_src_big(start)+BigFloat(i-1)*_src_big("1e-9")/BigFloat(86400)
                                    abs(e.epoch_mjd-reference)*BigFloat(86400)
                                end
                            end
                            @test err<=_src_big("1e-18")
                        end
                        @test source.settings.kind==:regular_schedule
                        @test source.settings.count==10
                        @test isempty(source.settings.selections)
                        @test precision(BigFloat)==80
                        @test rounding(BigFloat)==RoundDown
                    end
                end
            end
        end
        source=explicit_schedule([_src_big("55825.0"),_src_big("55825.0")];observatory="pks",frequency_mhz=1,uncertainty_us=1)
        @test length(unique(_src_ids(source)))==2
        @test _src_ids(explicit_schedule([55825];observatory="@",frequency_mhz=0,uncertainty_us=1))==
            _src_ids(explicit_schedule(["55825.0"];observatory="@",frequency_mhz=0,uncertainty_us=1))
    end

    @testset "Invalid schedules fail without writes" begin
        base=(observatory="pks",frequency_mhz=1400,uncertainty_us=1)
        @test_throws ArgumentError explicit_schedule(String[];base...)
        @test_throws ArgumentError explicit_schedule([55825.0];base...)
        @test_throws ArgumentError explicit_schedule(["55825.0"];observatory="pks")
        @test_throws ArgumentError explicit_schedule(["55825.0"];base...,backend=["A","B"])
        @test_throws ArgumentError explicit_schedule(["55825.0"];base...,backend="A",flags=["-be"=>"A"])
        @test_throws ArgumentError explicit_schedule(["55825.0"];base...,flags=["-x"])
        @test_throws ArgumentError explicit_schedule(["55825.0"];observatory="two sites",frequency_mhz=1,uncertainty_us=1)
        @test_throws TimFormatError explicit_schedule(["55825.0"];observatory="pks",frequency_mhz=-1,uncertainty_us=1)
        @test_throws TimFormatError explicit_schedule(["55825.0"];observatory="pks",frequency_mhz=1,uncertainty_us=0)
        @test_throws ArgumentError explicit_schedule(["55825.0"];backend="missing",systems=Dict())
        @test_throws ArgumentError explicit_schedule(["55825.0"];base...,backend="A",systems=Dict("A"=>base))
        @test_throws ArgumentError regular_schedule("55825.0";count=0,cadence_s=1,base...)
        for cadence in (0,-1,NaN,Inf)
            @test_throws ArgumentError regular_schedule("55825.0";count=2,cadence_s=cadence,base...)
        end
        @test_throws ArgumentError regular_schedule("99999.9";count=2,cadence_s=86400,base...)
        @test_throws ArgumentError explicit_schedule(["55825.0"];base...,precision_bits=80)
        @test_throws ArgumentError create_tim_document([(epoch_mjd="55825.0",)])
        @test_throws ArgumentError explicit_schedule(["55825.0","55825.1"];base...,filenames=["obs","#hidden"])
        @test_throws ArgumentError explicit_schedule(["55825.0"];base...,filenames="C")
        @test_throws ArgumentError create_tim_document([(epoch_mjd="55825.0",observatory="pks",
            frequency_mhz=1400,uncertainty_us=1,flags=Dict("-be"=>"A"))])
    end
end
