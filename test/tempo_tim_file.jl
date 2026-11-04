using Test
using GravityToolsNext
using SHA

_tim_test_big(s, bits=256) = BigFloat(s, RoundNearest; precision=bits)
_tim_test_row(epoch="55825.123456789012345678901234"; name="obs", flags="") =
    "$name 1400.00 $epoch 2.500 pks$flags"

function _tim_test_document(dir, text; name="source.tim", bits=128)
    path = joinpath(dir, name)
    write(path, text)
    return read_tim_document(path; precision_bits=bits)
end

@testset "TIM document I/O" begin
    @testset "Byte preservation, metadata and identity" begin
        mktempdir() do dir
            row = _tim_test_row(; flags=" -be A -x one -x two -pn -900719925474099300000 -pn +00042 -q -1.25e-3")
            for (i, (eol, mode, final)) in enumerate(Iterators.product(("\n", "\r\n"), ("", "MODE 0", "MODE 1"), (false, true)))
                lines = ["C comment", "# UTF-8 café comment", " \t", "FORMAT 1"]
                isempty(mode) || push!(lines, mode)
                append!(lines, [row, row])
                text = join(lines, eol) * (final ? eol : "")
                doc = _tim_test_document(dir, text; name="input$i.tim")
                entries = toa_entries(doc)
                @test length(entries) == 2
                @test entries[1].source_line == length(lines)-1
                @test entries[1].ordinal == 1 && entries[2].ordinal == 2
                @test entries[1].id != entries[2].id
                @test length(Set(e.id for e in entries)) == 2
                @test entries[1].id.source_sha256 == bytes2hex(sha256(text))
                @test toa_entries(read_tim_document(doc.source_path))[1].id == entries[1].id
                @test entries[1].freq_mhz == 1400
                @test entries[1].uncertainty_us == 2.5
                @test entries[1].observatory == "pks"
                @test entries[1].tokens[2] == "1400.00"
                @test entries[1].flags == ["-be"=>"A", "-x"=>"one", "-x"=>"two",
                    "-pn"=>"-900719925474099300000", "-pn"=>"+00042", "-q"=>"-1.25e-3"]
                @test precision(entries[1].epoch_mjd) == 128
                output = write_tim_document(doc, joinpath(dir, "output$i.tim"))
                @test read(output) == Vector{UInt8}(codeunits(text))
                @test read(doc.source_path, String) == text
            end
            # A FORMAT-only header must not lose the first TOA like the legacy reader.
            doc = _tim_test_document(dir, "FORMAT 1\n" * _tim_test_row(; name="café") * "\n")
            @test only(toa_entries(doc)).filename == "café"
            empty_doc = _tim_test_document(dir, "FORMAT 1"; name="empty.tim")
            @test isempty(toa_entries(empty_doc))
            @test read(write_tim_document(empty_doc, joinpath(dir,"empty_copy.tim")), String) == "FORMAT 1"
            zero_freq = _tim_test_document(dir, "FORMAT 1\nobs 0 55825.0 1 @"; name="infinite.tim")
            @test only(toa_entries(zero_freq)).freq_mhz == 0
        end
    end

    @testset "Only the epoch token changes" begin
        mktempdir() do dir
            original_epoch = "55825.123456789012345678901234567890"
            text = "FORMAT 1\r\n\t obs\t 1400.00  $original_epoch\t2.500   pks  -be A\t-pn -0042  \r\n" *
                _tim_test_row() * "\n"
            doc = _tim_test_document(dir, text)
            before = toa_entries(doc)
            target = _tim_test_big("55825.123456789012357252975308641964")
            @test set_toa_epoch!(doc, before[1].id, target) === doc
            after = toa_entries(doc)
            @test after[1].id == before[1].id
            @test after[2].epoch_mjd == before[2].epoch_mjd
            @test after[1].flags == before[1].flags
            @test after[1].tokens == before[1].tokens
            @test before[1].epoch_mjd != after[1].epoch_mjd # snapshots do not change
            out = write_tim_document(doc, joinpath(dir,"edited.tim"))
            serialized = read(out,String)
            reread = toa_entries(read_tim_document(out))
            token = reread[1].tokens[3]
            @test length(split(token,'.')[2]) == 30
            @test serialized == replace(text, original_epoch=>token; count=1)
            @test reread[1].id != after[1].id # fresh document, new content hash
            @test read(doc.source_path,String) == text
            # Mutating the returned metadata vectors cannot change the document.
            after[1].flags[1] = "-be"=>"changed"
            after[1].tokens[3] = "broken"
            @test toa_entries(doc)[1].flags[1] == ("-be"=>"A")
            @test read(write_tim_document(doc,joinpath(dir,"edited_again.tim")),String) == serialized
        end
    end

    @testset "Explicit precision and nanosecond site delays" begin
        mktempdir() do dir
            for ambient in (80, 192), bits in (128, 256),
                epoch in ("10000.0000000000001", "55825.123456789012345678901234", "99999.9999999999999"),
                sign in (-1, 1)
                setprecision(BigFloat, ambient) do
                    setrounding(BigFloat, RoundDown) do
                        doc = _tim_test_document(dir, "FORMAT 1\n"*_tim_test_row(epoch); bits=bits)
                        entry = only(toa_entries(doc))
                        @test entry.epoch_mjd == _tim_test_big(epoch,bits)
                        delay = sign * 1e-9
                        reference = setprecision(BigFloat,256) do
                            setrounding(BigFloat,RoundNearest) do
                                _tim_test_big(epoch) + BigFloat(delay)/BigFloat(86400)
                            end
                        end
                        shift_toa_epoch!(doc,entry.id,delay)
                        edited = only(toa_entries(doc))
                        @test sign*(edited.epoch_mjd-entry.epoch_mjd) > 0
                        output = write_tim_document(doc,joinpath(dir,"shifted.tim"); overwrite=true)
                        recovered = only(toa_entries(read_tim_document(output; precision_bits=bits)))
                        error_s = setprecision(BigFloat,256) do
                            abs(recovered.epoch_mjd-reference)*BigFloat(86400)
                        end
                        @test error_s <= _tim_test_big("1e-18")
                        @test precision(BigFloat) == ambient
                        @test rounding(BigFloat) == RoundDown
                    end
                end
            end
        end
    end

    @testset "1000 real write/read cycles" begin
        mktempdir() do dir
            for bits in (128,256)
                doc = _tim_test_document(dir,"FORMAT 1\n"*_tim_test_row(); bits=bits)
                reference = setprecision(BigFloat,256) do
                    _tim_test_big(only(toa_entries(doc)).tokens[3]) + _tim_test_big("1e-9")/BigFloat(86400)
                end
                shift_toa_epoch!(doc,only(toa_entries(doc)).id,_tim_test_big("1e-9"))
                canonical = ""
                for i in 1:1000
                    # Alternate targets because each new source is protected.
                    path=joinpath(dir, isodd(i) ? "cycle_a.tim" : "cycle_b.tim")
                    write_tim_document(doc,path; overwrite=true)
                    next_doc=read_tim_document(path; precision_bits=bits)
                    entry=only(toa_entries(next_doc))
                    error_s=setprecision(BigFloat,256) do
                        abs(entry.epoch_mjd-reference)*BigFloat(86400)
                    end
                    @test error_s <= _tim_test_big("1e-15")
                    i == 1 ? (canonical=entry.tokens[3]) : (@test entry.tokens[3] == canonical)
                    # Exercise edited serialization on every cycle, not only raw passthrough.
                    set_toa_epoch!(next_doc,entry.id,entry.epoch_mjd)
                    doc=next_doc
                end
                @test canonical == "55825.123456789012357252975308"
            end
        end
    end

    @testset "Boundary formatting and accumulated site shifts" begin
        mktempdir() do dir
            for bits in (128,256)
                doc=_tim_test_document(dir,"FORMAT 1\n"*_tim_test_row("55825.0"; name="café"); bits=bits)
                id=only(toa_entries(doc)).id
                for _ in 1:1000
                    shift_toa_epoch!(doc,id,_tim_test_big("1e-9"))
                end
                reference=setprecision(BigFloat,256) do
                    BigFloat(55825)+_tim_test_big("1e-6")/BigFloat(86400)
                end
                out=joinpath(dir,"accumulated.tim")
                write_tim_document(doc,out;overwrite=true)
                recovered=only(toa_entries(read_tim_document(out;precision_bits=bits)))
                @test setprecision(BigFloat,256) do
                    abs(recovered.epoch_mjd-reference)*BigFloat(86400) < _tim_test_big("1e-18")
                end
                # The shortest allowed input must also support the largest
                # representable in-domain epoch without rounding it to 100000.
                edge=prevfloat(_tim_test_big("100000",bits))
                set_toa_epoch!(doc,id,edge)
                write_tim_document(doc,out;overwrite=true)
                recovered=only(toa_entries(read_tim_document(out;precision_bits=bits)))
                @test recovered.epoch_mjd == edge
                @test length(split(recovered.tokens[3],'.')[2]) >= 24
                @test recovered.filename == "café"
            end
        end
    end

    @testset "Strict errors with physical lines" begin
        mktempdir() do dir
            directives = ("INCLUDE missing.tim", "END", "SKIP", "NOSKIP", "TIME 1", "PHASE 1", "JUMP",
                "EFAC 2", "EQUAD 1", "GLOBAL_EFAC 2", "T2EFAC -be A 2", "T2EQUAD -be A 1",
                "SIGMA 1", "ESET 1", "EMAX 1", "EFLOOR 1", "FMIN 1", "FMAX 1", "INFO a", "PROFILE_DIR /missing", "HEAD")
            bad_rows = ("I "*_tim_test_row(), _tim_test_row(;name="Crab"), _tim_test_row(;name="Island"),
                "c comment", " # indented", " C comment", "obs 1400 55825.0", _tim_test_row("55825"),
                _tim_test_row("5.5825e4"), _tim_test_row("+55825.0"), _tim_test_row("9999.9"),
                _tim_test_row("100000.0"), _tim_test_row("NaN"), _tim_test_row("Inf"),
                "obs NaN 55825.0 1 pks", "obs -1 55825.0 1 pks", "obs 1 55825.0 0 pks",
                "obs 1 55825.0 Inf pks", "obs 1 55825.0 -1 pks", _tim_test_row(;flags=" -be"),
                _tim_test_row(;flags=" -be -x"), _tim_test_row(;flags=" -pn 1.0"),
                _tim_test_row(;flags=" -pn 1e3"), _tim_test_row(;flags=" -x \"quoted\""),
                _tim_test_row(;flags=" # inline"), _tim_test_row(;flags=" -x \"two words\""),
                _tim_test_row()*"\r", "MODE 2", "FORMAT 2", "FORMAT 1")
            for row in (directives...,bad_rows...)
                text="FORMAT 1\n# note\n"*row*"\n"
                # Bare CR is unsupported except as part of CRLF.
                endswith(row,'\r') && (text="FORMAT 1\n# note\n"*row)
                path=joinpath(dir,"invalid.tim");write(path,text)
                err=try read_tim_document(path); nothing catch e; e end
                @test err isa TimFormatError
                if err isa TimFormatError
                    @test err.line == 3
                    @test occursin("$path:3:",sprint(showerror,err))
                end
                @test read(path,String) == text
            end
            for text in ("", "MODE 1\n", _tim_test_row()*"\n", "format 1\n")
                @test_throws TimFormatError _tim_test_document(dir,text)
            end
            @test_throws TimFormatError _tim_test_document(dir,"FORMAT 1\n"*_tim_test_row()*"\nMODE 1\n")
            @test_throws TimFormatError _tim_test_document(dir,"FORMAT 1\nMODE 1\nMODE 0\n")
            @test_throws TimFormatError _tim_test_document(dir,"FORMAT 1\nobs\0 1 55825.0 1 pks")
            @test_throws ArgumentError read_tim_document(joinpath(dir,"invalid.tim"); precision_bits=80)
        end
    end

    @testset "Failed edits restore state and MPFR context" begin
        mktempdir() do dir
            doc=_tim_test_document(dir,"FORMAT 1\n"*_tim_test_row())
            entry=only(toa_entries(doc))
            setprecision(BigFloat,80) do
                setrounding(BigFloat,RoundUp) do
                    @test_throws TimFormatError set_toa_epoch!(doc,entry.id,_tim_test_big("9999"))
                    @test_throws TimFormatError set_toa_epoch!(doc,entry.id,_tim_test_big("100000"))
                    @test_throws TimFormatError set_toa_epoch!(doc,entry.id,_tim_test_big("NaN"))
                    @test_throws ArgumentError set_toa_epoch!(doc,entry.id,BigFloat(55825))
                    @test_throws ArgumentError set_toa_epoch!(doc,entry.id,55825.0)
                    @test_throws ArgumentError shift_toa_epoch!(doc,entry.id,Inf)
                    @test_throws TimFormatError shift_toa_epoch!(doc,entry.id,-1e20)
                    @test_throws ArgumentError shift_toa_epoch!(doc,TimTOAID("other",1),1e-9)
                    @test only(toa_entries(doc)).epoch_mjd == entry.epoch_mjd
                    @test precision(BigFloat) == 80
                    @test rounding(BigFloat) == RoundUp
                    @test_throws TimFormatError _tim_test_document(dir,"FORMAT 1\n"*_tim_test_row("55825e0"); name="bad.tim")
                    @test precision(BigFloat) == 80
                    @test rounding(BigFloat) == RoundUp
                end
            end
        end
    end

    @testset "Source aliases and non-destructive publication" begin
        mktempdir() do dir
            original="FORMAT 1\n"*_tim_test_row()
            doc=_tim_test_document(dir,original)
            source=doc.source_path
            link=joinpath(dir,"link.tim"); symlink(source,link)
            hard=joinpath(dir,"hard.tim"); hardlink(source,hard)
            parent_alias=joinpath(dir,"parent"); symlink(dir,parent_alias)
            for overwrite in (false,true), path in (source,joinpath(dir,".","source.tim"),link,hard,joinpath(parent_alias,"source.tim"))
                @test_throws ArgumentError write_tim_document(doc,path;overwrite=overwrite)
                @test read(source,String) == original
            end
            # Reading through a symlink also protects both its path and its target.
            linked_doc=read_tim_document(link)
            @test_throws ArgumentError write_tim_document(linked_doc,source;overwrite=true)
            @test_throws ArgumentError write_tim_document(linked_doc,link;overwrite=true)
            # The saved inode protects original hardlinks even after source replacement.
            moved=joinpath(dir,"moved.tim"); mv(source,moved)
            write(source,original*"\n")
            @test_throws ArgumentError write_tim_document(doc,moved;overwrite=true)
            @test_throws ArgumentError write_tim_document(doc,hard;overwrite=true)
            target=joinpath(dir,"output.tim"); write(target,"KEEP")
            @test_throws ArgumentError write_tim_document(doc,target)
            @test read(target,String) == "KEEP"
            write_tim_document(doc,target;overwrite=true)
            @test read(target,String) == original
            broken=joinpath(dir,"broken.tim"); symlink(joinpath(dir,"absent.tim"),broken)
            @test_throws ArgumentError write_tim_document(doc,broken;overwrite=true)
            @test !ispath(joinpath(dir,"absent.tim"))
            @test_throws ArgumentError write_tim_document(doc,dir;overwrite=true)
            before=Set(readdir(dir))
            @test_throws Exception write_tim_document(doc,joinpath(dir,repeat("x",300)))
            @test Set(readdir(dir)) == before # includes cleanup after publish failure
            @test_throws Exception write_tim_document(doc,joinpath(dir,"missing","out.tim"))
            @test !isdir(joinpath(dir,"missing"))
            # Validate a corrupt candidate before touching an existing target.
            push!(doc.nodes, GravityToolsNext.TimRawLine(:unknown,"\nINCLUDE missing.tim\n"))
            before=Set(readdir(dir))
            @test_throws TimFormatError write_tim_document(doc,target;overwrite=true)
            @test read(target,String) == original
            @test Set(readdir(dir)) == before
        end
    end

    @testset "Legacy analysis reader remains separate" begin
        mktempdir() do dir
            text="FORMAT 1\nMODE 1\n"*_tim_test_row(;flags=" -be PDFB1")*"\n"
            doc=_tim_test_document(dir,text)
            legacy=only(read_tim_file(doc.source_path))
            @test legacy isa TimTOAEntry
            @test legacy.toa isa Float64
            @test legacy.backend == "PDFB1"
            @test only(toa_entries(doc)).epoch_mjd isa BigFloat
        end
    end
end
