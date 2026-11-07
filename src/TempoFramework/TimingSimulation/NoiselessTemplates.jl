"""
    IdealizationSettings(; tolerance_s=1e-9, max_corrections=5, timeout_s=60,
                          precision_bits=128, max_step_phase=0.01,
                          max_bootstrap_corrections=3,
                          mode_switch_tolerance_s=1e-14)

Bounded local mean-reference idealization. Both arithmetic-centered and
weighted-centered residual maxima must meet `tolerance_s`. Each engine query
has a POSIX deadline; at most five corrections are permitted. The dimensionless
step limit is a conservative local-domain guard, not a universal phase guarantee;
engine-derived pulse numbers must also remain unchanged between queries.
Generated schedules first use a separate MODE 0 initialization with a fixed
first-selected-TOA anchor. Its half-cycle margin is recorded in cycles as a
diagnostic, without an initial rejection threshold. Final acceptance
always uses MODE 1. These limits do not trigger bootstrap for existing TIMs.
"""
struct IdealizationSettings
    tolerance_s::Float64
    max_corrections::Int
    timeout_s::Int
    precision_bits::Int
    max_step_phase::Float64
    max_bootstrap_corrections::Int
    mode_switch_tolerance_s::Float64
    function IdealizationSettings(;tolerance_s::Real=1e-9,max_corrections::Integer=5,
            timeout_s::Integer=60,precision_bits::Integer=128,max_step_phase::Real=0.01,
            max_bootstrap_corrections::Integer=3,
            mode_switch_tolerance_s::Real=1e-14)
        isfinite(tolerance_s) && 0<tolerance_s && 0<Float64(tolerance_s) && isfinite(Float64(tolerance_s)) ||
            throw(ArgumentError("tolerance_s must be positive and finite"))
        !(max_corrections isa Bool) && 0<=max_corrections<=5 || throw(ArgumentError("max_corrections must be in 0:5"))
        !(timeout_s isa Bool) && 1<=timeout_s<=60 || throw(ArgumentError("timeout_s must be in 1:60"))
        _tim_bits(precision_bits)
        isfinite(max_step_phase) && 0<max_step_phase<=0.1 || throw(ArgumentError("max_step_phase must be in (0,0.1]"))
        !(max_bootstrap_corrections isa Bool) && 0<=max_bootstrap_corrections<=3 ||
            throw(ArgumentError("max_bootstrap_corrections must be in 0:3"))
        isfinite(mode_switch_tolerance_s) && 0<mode_switch_tolerance_s<=1e-14 && 0<Float64(mode_switch_tolerance_s)<=1e-14 ||
            throw(ArgumentError("mode_switch_tolerance_s must be in (0,1e-14]"))
        new(Float64(tolerance_s),Int(max_corrections),Int(timeout_s),Int(precision_bits),Float64(max_step_phase),
            Int(max_bootstrap_corrections),Float64(mode_switch_tolerance_s))
    end
end

"""A failed idealization. `workspace` retains inputs, partial output and diagnostics.
No unconverged or failed engine output is returned as a successful template.
"""
struct IdealizationError <: Exception
    status::Symbol
    message::String
    workspace::String
end
Base.showerror(io::IO,e::IdealizationError)=print(io,e.status,": ",e.message," (workspace: ",e.workspace,")")

"""
    NoiselessTemplate

A verified serialized TIM (`document`), original `lineage`, inspectable
`metadata`, and owned artifact bytes. Create with [`create_noiseless_template`](@ref).
The harmless reserved `-gtnid` flags are deliberately retained in the final TIM;
they are collision-checked and cannot be deterministic JUMP selectors.
Saving/loading is explicit and does not execute TEMPO2. Treat fields as read-only;
saving rejects a document or artifact set modified after verification.
"""
struct NoiselessTemplate
    document::TimDocument
    lineage::Vector{TOALineage}
    metadata::Dict{String,Any}
    artifacts::Dict{String,Vector{UInt8}}
end

const _TEMPLATE_ID_FLAG="-gtnid"
const _TEMPLATE_METADATA_FLAGS=Set(split("-bw -gof -tobs -nbin -fe -tmplt -be -nch -epoch -f -pn"))
const _TEMPLATE_ACTIVE_FLAGS=Set(split("-pnadd -padd -addsat -to -phase -time -dm -dmo -dme -efac -equad"))
const _TEMPLATE_NOISE_GROUPS=Set(["TNEF","TNEQ","T2EFAC","T2EQUAD"])
const _TEMPLATE_NOISE_SCALARS=Set(["TNGLOBALEF","TNGLOBALEQ","ADDTNGLOBALEQ",
    "TNREDAMP","TNREDGAM","TNREDC","TNREDFLOW","TNREDCORNER","TNSUBTRACTRED"])
const _TEMPLATE_DIAGNOSTICS=Set(["CHI2","CHI2R","NTOA","TRES"])

_template_sha(bytes)=bytes2hex(sha256(bytes))

# Strictly validate the subset before using the existing PAR reader: it otherwise
# warns/skips malformed lines and coalesces duplicate parameter names.
function _template_read_par(text::String)
    parameters=TempoParameter[]; seen=Set{Symbol}()
    for (line_no,raw) in enumerate(split(text,'\n'))
        line=strip(_strip_inline_comment(raw))
        (isempty(line) || occursin(r"^[Cc](?:\s|$)",line)) && continue
        tokens=split(line)
        p=extract_tempo_parameter_from_line(line)
        words=split(p.name)
        p.name=join(vcat(uppercase(words[1]),words[2:end])," ")
        _TEMPLATE_ID_FLAG in words && throw(ArgumentError("Reserved diagnostic flag used by a model selector"))
        p.value===nothing && throw(ArgumentError("Unsupported PAR syntax at line $line_no"))
        p.name_symbol in seen && throw(ArgumentError("Duplicate PAR parameter $(p.name)"))
        push!(seen,p.name_symbol)
        name_tokens=length(split(p.name))
        represented=1+(p.flag===nothing ? 0 : 1)+(p.uncertainty===nothing ? 0 : 1)
        length(tokens)==name_tokens+represented || throw(ArgumentError("Unparsed PAR fields at line $line_no"))
        length(tokens)<=name_tokens+3 || throw(ArgumentError("Unsupported multi-value PAR line $line_no"))
        p.value isa Real && !isfinite(p.value) && throw(ArgumentError("Nonfinite PAR value $(p.name)"))
        if uppercase(first(tokens))=="JUMP"
            length(tokens) in (4,5,6) && startswith(tokens[2],"-") &&
                tryparse(BigFloat,tokens[3])===nothing && p.value isa Real ||
                throw(ArgumentError("Only JUMP -flag string_selector value [fit [uncertainty]] is supported"))
            tokens[2]==_TEMPLATE_ID_FLAG && throw(ArgumentError("Reserved diagnostic flag used by a JUMP"))
        end
        push!(parameters,p)
    end
    isempty(parameters) && throw(ArgumentError("Empty reference PAR"))
    return parameters
end

function _template_noise_kind(p::TempoParameter)
    h=uppercase(first(split(p.name)))
    h in _TEMPLATE_DIAGNOSTICS && return :diagnostic
    h in _TEMPLATE_NOISE_GROUPS && return :noise
    h in _TEMPLATE_NOISE_SCALARS && return :noise
    unsupported=startswith(h,"TN") || startswith(h,"ADDTN") || startswith(h,"WAVE") ||
        occursin("IFUNC",h) || occursin("ECORR",h) ||
        h in ("RNAMP","RNIDX","EFAC","EQUAD","DMMODEL","CMMODEL","DMOFF","CMOFF",
            "CONSTRAIN","CONSTRAINT","CONSTRAINT_EFACTOR","WHITENOISEMODEL","COVARFUNC")
    unsupported && throw(ArgumentError("Unsupported noise/whitening control $(p.name)"))
    return :deterministic
end

function _prepare_template_par(original::String)
    _tim_precision(256) do
        ps=_template_read_par(original)
        model=Dict(uppercase(p.name)=>p for p in ps)
        get(model,"BINARY",TP("BINARY","")).value=="DDSTG" &&
            get(model,"MODEL_MODE",TP("MODEL_MODE","")).value=="GR" ||
            throw(ArgumentError("Noiseless templates currently support the validated DDSTG GR model only"))
        haskey(model,"TRACK") && throw(ArgumentError("Active TRACK is not supported"))
        for p in ps
            h=uppercase(first(split(p.name)))
            h in ("PHASE","PHASEJUMP","PULSE","NPHASE","PNADD") && throw(ArgumentError("Unsupported phase control $(p.name)"))
            h in ("START","FINISH") && p.flag ∉ (nothing,0) &&
                throw(ArgumentError("Active START/FINISH selection is not supported"))
            h=="MODE" && p.value!=1 && throw(ArgumentError("Only weighted MODE 1 is supported"))
        end
        f0=get(model,"F0",TP("F0",0)).value
        f0 isa Real && f0>0 && isfinite(Float64(f0)) && Float64(f0)>0 || throw(ArgumentError("Positive F0 is required"))
        transformed=String[]; kept=TempoParameter[]; selectors=Set{String}()
        for p in ps
            kind=_template_noise_kind(p)
            if kind!==:deterministic
                p.value isa Real || throw(ArgumentError("Invalid numeric noise/diagnostic value $(p.name)"))
                push!(transformed,"Remove $(kind): $(get_par_file_representation(p))")
                continue
            end
            startswith(uppercase(p.name),"JUMP ") && push!(selectors,split(p.name)[2])
            q=deepcopy(p)
            if q.flag!==nothing && q.flag!=0
                push!(transformed,"Clear fit flag $(q.name): $(q.flag) -> 0")
                q.flag=0
            end
            push!(kept,q)
        end
        pf=TempoParFile("reference.par","",Dict(p.name_symbol=>p for p in kept),[p.name_symbol for p in kept])
        for p in [TP("REFPHS","MEAN"),TP("MODE",1),TP("NITS",1),TP("GAIN",1)]
            push!(transformed,"Set $(p.name) = $(p.value)")
            upsert_param!(pf,p)
        end
        text=join((get_par_file_representation(p) for p in params_as_vector(pf)),"\n")*"\n"
        roundtrip=Dict(p.name=>p for p in _template_read_par(text))
        for p in kept
            p.name in ("REFPHS","MODE","NITS","GAIN") && continue
            p.value==roundtrip[p.name].value || throw(ArgumentError("PAR serialization changes $(p.name); unsupported precision/syntax"))
        end
        return (text=text,transformations=transformed,selectors=selectors,f0_hz=Float64(f0))
    end
end

function _validate_template_source(source::ObservationSource,settings,selectors)
    (source.settings.kind==:tim && source.document.origin==:file ||
        source.settings.kind in (:explicit_schedule,:regular_schedule) && source.document.origin==:generated) ||
        throw(ArgumentError("Expected a file-origin TIM or genuine explicit/regular schedule"))
    source.settings.origin==source.document.origin || throw(ArgumentError("Source origin mismatch"))
    entries=toa_entries(source.document)
    1<=length(entries)<=128 || throw(ArgumentError("Initial template support requires 1:128 TOAs"))
    length(source.lineage)==length(entries) && length(unique(e.id for e in entries))==length(entries) ||
        throw(ArgumentError("Invalid observation lineage/identities"))
    [x.output_index for x in source.lineage]==collect(eachindex(entries)) || throw(ArgumentError("Invalid lineage order"))
    [x.source_id for x in source.lineage]==[e.id for e in entries] || throw(ArgumentError("Source IDs differ from lineage"))
    settings.precision_bits>=source.document.precision_bits || throw(ArgumentError("Template precision must not reduce source precision"))
    for node in source.document.nodes
        source.document.origin==:file && node isa TimRawLine && node.kind==:mode && occursin(r"MODE\s+0",node.raw) &&
            throw(ArgumentError("TIM MODE 0 is unsupported"))
    end
    allowed=union(_TEMPLATE_METADATA_FLAGS,selectors)
    for e in entries
        length(e.flags)<=32 || throw(ArgumentError("Too many flags for the validated stock output path"))
        for token in e.tokens
            ncodeunits(token)<=80 && !occursin('|',token) || throw(ArgumentError("Unsupported long/output-delimited TIM token"))
        end
        for (flag,value) in e.flags
            lowercase(flag) in _TEMPLATE_ACTIVE_FLAGS && throw(ArgumentError("Unsupported active timing flag $flag"))
            flag==_TEMPLATE_ID_FLAG && throw(ArgumentError("Reserved diagnostic flag collision"))
            flag in allowed || throw(ArgumentError("Unvalidated timing/metadata flag $flag"))
        end
    end
    return entries
end
