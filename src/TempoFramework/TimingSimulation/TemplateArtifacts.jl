const _NOISELESS_TEMPLATE_FORMAT = 3

function _template_new_directory(path::AbstractString)
    target=abspath(path)
    (ispath(target) || islink(target)) && throw(ArgumentError("Target already exists: $target"))
    isdir(dirname(target)) || throw(ArgumentError("Target parent must exist"))
    mkdir(target) # Atomic create-only claim; never clean or reuse caller paths.
    return target
end

function _template_plain(value)
    value===nothing && return nothing
    value isa Union{AbstractString,Bool,Integer,Float64} && return value
    value isa BigFloat && return string(value)
    value isa Symbol && return string(value)
    value isa NamedTuple && return Dict(string(k)=>_template_plain(v) for (k,v) in pairs(value))
    value isa AbstractDict && return Dict(string(k)=>_template_plain(v) for (k,v) in value)
    value isa Pair && return [_template_plain(first(value)),_template_plain(last(value))]
    value isa Union{AbstractVector,Tuple} && return [_template_plain(v) for v in value]
    value isa AbstractRNG && return Dict("type"=>string(typeof(value)),"state_repr"=>repr(value),
        "replay_policy"=>"Informational Julia-version-bound state; actual selected IDs are retained")
    return Dict("type"=>string(typeof(value)),"repr"=>repr(value))
end

function _template_file_identity(path)
    try
        return Dict{String,Any}("path"=>abspath(path),"realpath"=>realpath(path),"sha256"=>_template_sha(read(path)),"readable"=>true)
    catch err
        return Dict{String,Any}("path"=>String(path),"readable"=>false,"error"=>sprint(showerror,err))
    end
end

function _template_code_identity()
    root=realpath(joinpath(@__DIR__,"../../.."))
    hashes=Dict{String,Any}()
    for (dir,_,names) in walkdir(joinpath(root,"src"))
        for name in names
            p=joinpath(dir,name);hashes[relpath(p,root)]=_template_sha(read(p))
        end
    end
    for name in ["Project.toml","Manifest.toml"]
        isfile(joinpath(root,name)) && (hashes[name]=_template_sha(read(joinpath(root,name))))
    end
    gitinfo=Dict{String,Any}()
    for (name,args) in [("head",["rev-parse","HEAD"]),("status",["status","--porcelain"]),
                       ("tracked_source_diff",["diff","--binary","--","src","Project.toml"])]
        try
            # An explicit repository avoids ownership auto-discovery without
            # changing global Git settings. Installed source archives may have no Git.
            dotgit=joinpath(root,".git")
            options=isdir(dotgit) ? ["--git-dir="*dotgit,"--work-tree="*root] :
                ["-c","safe.directory="*root,"-C",root]
            gitinfo[name]=read(Cmd(vcat(["git"],options,args)),String)
        catch err
            gitinfo[name*"_error"]=sprint(showerror,err)
        end
    end
    return Dict("source_root"=>root,"file_sha256"=>hashes,"git"=>gitinfo,
        "policy"=>"Hashes include current untracked/dirty source; git diff alone does not capture new files")
end

_template_lineage_record(x::TOALineage)=Dict("output_index"=>x.output_index,"source_ordinal"=>x.source_ordinal,
    "source_sha256"=>x.source_id.source_sha256,"source_line"=>x.source_id.source_line,"source_origin"=>string(x.source_id.origin))

function _template_provenance(source,par_path,settings,tempo,resources,model,history,artifacts;bootstrap=nothing)
    exe=tempo_cmd_path(tempo)
    resolved=isabspath(exe) ? exe : something(Sys.which(exe),exe)
    runtime=tempo_data_dir(tempo)
    identities=[_template_file_identity(resolved)]
    push!(identities,_template_file_identity(joinpath(runtime,"lib","libtempo2.so.2")))
    plugdir=joinpath(runtime,"plugins")
    if isdir(plugdir)
        for n in sort(readdir(plugdir))
            startswith(n,"general2_") && endswith(n,"_plug.t2") && push!(identities,_template_file_identity(joinpath(plugdir,n)))
        end
    end
    final=last(history)
    return Dict{String,Any}("format_version"=>_NOISELESS_TEMPLATE_FORMAT,"initial_margin_policy"=>"diagnostic_only",
        "verified"=>true,"created_utc"=>string(Dates.now(Dates.UTC)),
        "julia_version"=>string(VERSION),"machine"=>Sys.MACHINE,"active_project"=>something(Base.active_project(),""),
        "settings"=>Dict(string(k)=>getfield(settings,k) for k in fieldnames(IdealizationSettings)),
        "source_ids"=>[_template_lineage_record(x) for x in source.lineage],
        "bootstrap"=>bootstrap,
        "source_settings"=>_template_plain(source.settings),"source_content_sha256"=>source.document.source_sha256,
        "source_origin"=>string(source.document.origin),"source_par_path"=>abspath(par_path),
        "source_epochs_mjd"=>[string(e.epoch_mjd) for e in toa_entries(source.document)],
        "final_epochs_mjd"=>[e.tokens[3] for e in toa_entries(_template_read_bytes(artifacts["final.tim"],settings.precision_bits))],
        "epoch_displacements_s"=>_tim_precision(256) do
            [string((b.epoch_mjd-a.epoch_mjd)*86400) for (a,b) in zip(toa_entries(source.document),
                toa_entries(_template_read_bytes(artifacts["final.tim"],settings.precision_bits)))]
        end,
        "diagnostic_ids"=>["T"*lpad(i,6,'0') for i in 1:length(source.lineage)],
        "flag_policy"=>"Retain harmless -gtnid in verified/persisted TIM; reserved flag cannot select any JUMP",
        "phase_policy"=>bootstrap===nothing ? "REFPHS MEAN; no TRACK; source -pn is inactive metadata; stable engine pulse numbers required" :
            "REFPHS MEAN; nearest branch relative to first selected generated TOA; fixed bootstrap anchor and pulse vector; no TRACK/pn rebasing; not absolute TZR locking",
        "residual_definition"=>"general2 post, seconds, postfit recomputation with unweighted arithmetic mean removal",
        "offset_definition"=>"param.vals prints fitted param_ZERO in seconds (%g); table weighted mean is a separate inference",
        "span_policy"=>"Unflagged START/FINISH retained in controlled init; engine may update extrema within each fresh process; all TOAs must remain present",
        "ddstg_policy"=>"DDSTG GR pre/post evaluations include derived/bookkeeping state; only param_ZERO may fit; never chain output models",
        "transformations"=>model.transformations,"history"=>history,"final_query"=>final["query"],
        "engine_identity"=>identities,"resource_identity"=>[_template_file_identity(p) for p in resources],
        "resource_provenance_complete"=>false,"resource_policy"=>"Explicit partial inventory, not an engine access trace; default clock/EOP/DDSTG/ephemeris dependencies are not exhaustively tracked",
        "code_identity"=>_template_code_identity(),"artifact_sha256"=>Dict(p=>_template_sha(v) for (p,v) in artifacts))
end

# Parse already-owned bytes without touching the old job or trusting its paths.
function _template_read_bytes(bytes,bits)
    return _parse_tim_text(String(copy(bytes)),"<noiseless template>",nothing,nothing,bits;origin=:generated)
end

function _template_safe_relative(name::AbstractString)
    !isempty(name) && !isabspath(name) && !occursin('\\',name) &&
        all(p->!isempty(p) && p ∉ (".",".."),split(name,'/')) || throw(ArgumentError("Unsafe bundle path: $name"))
    return String(name)
end

# Required current fields are never reconstructed from an older record shape.
function _template_require_fields(record,fields,context)
    _template_require(record isa AbstractDict,"Expected a dictionary for $context")
    for field in fields
        _template_require(haskey(record,field),"Incomplete $context: missing required field '$field'")
    end
end

function _template_settings_from_record(m)
    _template_require_fields(m,("initial_margin_policy","settings"),"template metadata")
    _template_require(m["initial_margin_policy"]=="diagnostic_only","Unsupported initial margin policy")
    raw=m["settings"];fields=string.(fieldnames(IdealizationSettings))
    _template_require_fields(raw,fields,"template settings")
    _template_require(Set(keys(raw))==Set(fields),"Settings do not match the current template format")
    return IdealizationSettings(; (Symbol(k)=>v for (k,v) in raw)...)
end

# Prove stage relationships from owned bytes without engine or workspace access.
function _validate_template_stages(template)
    m=template.metadata;a=template.artifacts;bits=m["settings"]["precision_bits"]
    settings=_template_settings_from_record(m)
    generated=m["source_origin"]=="generated"
    _template_require(m["source_origin"] in ("file","generated"),"Invalid source origin")
    _template_require(m["source_settings"]["origin"]==m["source_origin"] &&
        m["source_settings"]["kind"] in (generated ? ["regular_schedule","explicit_schedule"] : ["tim"]),"Invalid source kind/origin")
    _template_require(m["source_ids"]==[_template_lineage_record(x) for x in template.lineage],"Original source IDs differ")
    _template_require(all(x->string(x.source_id.origin)==m["source_origin"] &&
        x.source_id.source_sha256==m["source_content_sha256"],template.lineage),"Source identity namespace differs")
    src=_template_read_bytes(a["source.tim"],bits);source_entries=toa_entries(src)
    final_entries=toa_entries(template.document)
    _template_require(length(source_entries)==length(final_entries) && 1<=length(final_entries)<=128,"Source count differs")
    _template_require([string(e.epoch_mjd) for e in source_entries]==m["source_epochs_mjd"],"Requested epochs differ")
    _template_require(m["diagnostic_ids"]==["T"*lpad(i,6,'0') for i in eachindex(source_entries)],"Diagnostic IDs differ")
    for (i,(s,f)) in enumerate(zip(source_entries,final_entries))
        _template_require(s.tokens[[1,2,4,5]]==f.tokens[[1,2,4,5]] &&
            f.flags==vcat(s.flags,[_TEMPLATE_ID_FLAG=>m["diagnostic_ids"][i]]),"Observing metadata changed")
    end
    _tim_precision(256) do
        displacements=[string((f.epoch_mjd-s.epoch_mjd)*86400) for (s,f) in zip(source_entries,final_entries)]
        _template_require(displacements==m["epoch_displacements_s"],"Epoch displacement provenance differs")
    end
    model=_prepare_template_par(String(copy(a["original.par"])))
    _template_require(a["reference.par"]==Vector{UInt8}(codeunits(model.text)) && m["transformations"]==model.transformations,"Controlled model differs from original preparation")
    b=m["bootstrap"];history=m["history"]
    _template_require(history isa AbstractVector && !isempty(history),"Missing required query history")
    for h in history
        _template_require_fields(h,("query","mode","stage","corrections","input_tim","executed_input_tim",
            "fit_weight_definition","pulse_numbers","bat_mjd","engine_row_by_output","max_abs_s",
            "max_weighted_centered_s","weighted_mean_inferred_s","fitted_offset_printed_s"),"template history")
    end
    _template_require(length(unique(h["query"] for h in history))==length(history),"Invalid query history")
    _template_require(generated==(b!==nothing) && generated==haskey(a,"bootstrap.par"),"Bootstrap origin/model mismatch")
    if generated
        _template_require_fields(b,("anchor_output_index","anchor_source_id","anchor_epoch_mjd",
            "initial_pulse_numbers","initial_margin_cycles","entrance_limit_s","model_transformation",
            "tim_transformation","mode_switch_max_difference_s","mode0_query","mode1_query"),"template bootstrap")
        _template_require(a["bootstrap.par"]==Vector{UInt8}(codeunits(_template_reference_mode(model.text,0))),"Bootstrap differs beyond MODE")
        _template_require(b["anchor_output_index"]==1 && b["anchor_source_id"]==first(m["source_ids"]) &&
            b["anchor_epoch_mjd"]==first(m["source_epochs_mjd"]),"Invalid bootstrap anchor")
        _template_require(b["model_transformation"]=="Set MODE = 0 only; unchanged controlled physical/JUMP model" &&
            b["tim_transformation"]=="Remove TIM MODE directive from generated evaluator copies; preserve source bytes","Invalid bootstrap transformations")
    end
    previous=nothing;previous_bytes=nothing;initial_pulses=nothing;nboot=0;nref=0;first_weighted=nothing
    for (index,h) in enumerate(history)
        mode=h["mode"];stage=h["stage"];q=_template_safe_relative(h["query"])
        _template_require(q=="queries/q"*lpad(index,3,'0'),"Unexpected query sequence")
        if mode==0
            nboot+=1
            _template_require(generated && nref==0 && stage=="bootstrap" &&
                h["corrections"]==nboot-1 && nboot<=settings.max_bootstrap_corrections+1,"Invalid bootstrap bounds/stage")
        else
            nref+=1
            _template_require(mode==1 && stage=="refinement" && h["corrections"]==nref-1 &&
                nref<=settings.max_corrections+1,"Invalid refinement bounds/stage")
        end
        bytes=a[q*"/candidate.tim"]
        executed=q*"/tmp/"*basename(h["input_tim"])
        _template_require(h["executed_input_tim"]==executed,"Executed input pathname differs from query")
        _template_require(bytes==a[_template_safe_relative(executed)],"Executed candidate differs from archive")
        _template_require(bytes==a[_template_safe_relative(h["input_tim"])],"Archived/evaluated candidate differs")
        _template_require(a[q*"/tmp/model_init.par"]==a[mode==0 ? "bootstrap.par" : "reference.par"],"Query init model differs")
        doc=_template_read_bytes(bytes,bits);entries=toa_entries(doc)
        _template_require(length(entries)==length(final_entries),"Stage observation count differs")
        for (e,f) in zip(entries,final_entries)
            _template_require(e.tokens[[1,2,4,5]]==f.tokens[[1,2,4,5]] && e.flags==f.flags,"Stage observing metadata changed")
        end
        t=_template_residual_table(String(copy(a[q*"/tmp/rows.tsv"])),entries;previous_pulses=initial_pulses)
        printed=_template_fit_audit(String(copy(a[q*"/tmp/param.labels"])),String(copy(a[q*"/tmp/design.matrix"])),
            String(copy(a[q*"/tmp/constraints.matrix"])),String(copy(a[q*"/tmp/param.vals"])),t;mode)
        expected_weight=mode==0 ? "unit weights; offset audited against arithmetic mean" : "inverse nominal variance; offset audited against weighted mean"
        _template_require(h["fit_weight_definition"]==expected_weight && h["pulse_numbers"]==t.pulses &&
            h["bat_mjd"]==t.bat_mjd && h["engine_row_by_output"]==t.engine_row_by_output &&
            h["max_abs_s"]==t.max_abs_s && h["max_weighted_centered_s"]==t.max_weighted_centered_s &&
            h["weighted_mean_inferred_s"]==t.weighted_mean_s && h["fitted_offset_printed_s"]==printed,"Stage diagnostics differ")
        if previous===nothing
            initial_pulses=t.pulses
            _template_require(bytes==Vector{UInt8}(codeunits(_template_evaluator_text(src;generated))),"Initial evaluator differs from source")
            if generated
                margin=_template_margin(t,model.f0_hz)
                _template_require(isfinite(margin) && b["initial_pulse_numbers"]==t.pulses &&
                    b["initial_margin_cycles"]==string(margin),"Invalid initial phase margin/pulses")
                _template_require(b["entrance_limit_s"]==min(0.01/model.f0_hz,minimum(t.errors_s)/10),"Bootstrap entrance domain differs")
            end
        elseif mode==1 && nref==1 && generated
            first_weighted=h
            _template_require(bytes==previous_bytes && max(previous.max_abs_s,previous.max_weighted_centered_s)<=b["entrance_limit_s"],"Invalid MODE handoff candidate/domain")
            parity=_tim_precision(256) do;maximum(abs,t.precise_residuals_s.-previous.precise_residuals_s);end
            _template_require(b["mode_switch_max_difference_s"]==string(parity) && parity<=settings.mode_switch_tolerance_s &&
                b["mode0_query"]==history[index-1]["query"] && b["mode1_query"]==q,"Invalid MODE parity/provenance")
        else
            prev_entries=toa_entries(_template_read_bytes(previous_bytes,bits))
            _tim_precision(256) do
                r=previous.precise_residuals_s;delay=mode==0 ? r.-r[1] : r
                for (s,f,d) in zip(prev_entries,entries,delay)
                    _template_require(abs((f.epoch_mjd-s.epoch_mjd)*86400+d)<=BigFloat("1e-18"),"Epoch correction differs from recorded residual")
                end
            end
        end
        mode==1 && _template_require(t.max_abs_s*model.f0_hz<=settings.max_step_phase,"Refinement phase guard differs")
        previous=t;previous_bytes=bytes
    end
    _template_require(nref>0 && (!generated || nboot>0 && first_weighted!==nothing),"Missing required stage")
    _template_require(last(history)["mode"]==1 && m["final_query"]==last(history)["query"] &&
        previous_bytes==a["final.tim"],"Final artifact is not strict verified candidate")
    return nothing
end

function _validate_noiseless_template(template::NoiselessTemplate)
    m=template.metadata;a=template.artifacts
    _template_require_fields(m,("format_version",),"template metadata")
    _template_require(m["format_version"] isa Integer && !(m["format_version"] isa Bool) &&
        m["format_version"]==_NOISELESS_TEMPLATE_FORMAT,"Unsupported template format: only complete format 3 is supported")
    _template_require_fields(m,("verified","initial_margin_policy","settings","source_origin","source_settings",
        "source_ids","source_content_sha256","source_epochs_mjd","final_epochs_mjd","epoch_displacements_s",
        "diagnostic_ids","bootstrap","transformations","history","final_query","artifact_sha256"),"template metadata")
    _template_require(m["verified"]===true,"Unverified template")
    _template_settings_from_record(m)
    _template_require_fields(m["source_settings"],("origin","kind"),"template source settings")
    _template_require(Set(keys(a))==Set(keys(m["artifact_sha256"])),"Artifact inventory differs")
    for (name,bytes) in a
        _template_safe_relative(name)
        _template_require(_template_sha(bytes)==m["artifact_sha256"][name],"Artifact checksum mismatch: $name")
    end
    _template_require(Vector{UInt8}(codeunits(_tim_render(template.document)))==a["final.tim"],"Template document changed after verification")
    entries=toa_entries(template.document)
    _template_require(length(entries)==length(template.lineage)==length(m["diagnostic_ids"]),"Invalid template lineage count")
    _template_require([x.output_index for x in template.lineage]==collect(eachindex(entries)) &&
        length(unique(x.source_id for x in template.lineage))==length(entries),"Invalid template lineage IDs/order")
    _template_require([e.tokens[3] for e in entries]==m["final_epochs_mjd"],"Final epoch strings differ")
    _validate_template_stages(template)
    q=_template_safe_relative(m["final_query"])
    _template_require(a[q*"/tmp/model_init.par"]==a["reference.par"],"Final engine init differs from reference")
    h=last(m["history"])
    table=_template_residual_table(String(copy(a[q*"/tmp/rows.tsv"])),entries;previous_pulses=h["pulse_numbers"])
    printed=_template_fit_audit(String(copy(a[q*"/tmp/param.labels"])),String(copy(a[q*"/tmp/design.matrix"])),
        String(copy(a[q*"/tmp/constraints.matrix"])),String(copy(a[q*"/tmp/param.vals"])),table)
    tolerance=m["settings"]["tolerance_s"]
    _template_require(isfinite(tolerance) && tolerance>0 && table.max_abs_s<=tolerance && table.max_weighted_centered_s<=tolerance,"Template does not meet tolerance")
    _template_require(h["max_abs_s"]==table.max_abs_s && h["max_weighted_centered_s"]==table.max_weighted_centered_s &&
        h["weighted_mean_inferred_s"]==table.weighted_mean_s && h["fitted_offset_printed_s"]==printed,"Final metrics inconsistent")
    return template
end

"""
    save_noiseless_template(template, target_directory) -> String

Save a self-contained create-only bundle, including exact TIM/PARs, lineage,
provenance and query artifacts. No engine is called. The parent must exist;
existing targets and symlinks are rejected. Partial writes remain incomplete and
cannot be loaded. JLD2 stores the current primitive record (format marker 3); checksums cover every
file. Caller-owned template bytes/document are validated before saving.
"""
function save_noiseless_template(template::NoiselessTemplate,target::AbstractString)
    _validate_noiseless_template(template)
    root=_template_new_directory(target)
    record=Dict{String,Any}("version"=>_NOISELESS_TEMPLATE_FORMAT,"metadata"=>deepcopy(template.metadata),
        "lineage"=>[_template_lineage_record(x) for x in template.lineage])
    hashes=Dict{String,String}()
    for (name,bytes) in template.artifacts
        _template_safe_relative(name)
        p=joinpath(root,name);mkpath(dirname(p));write(p,bytes);hashes[name]=_template_sha(bytes)
    end
    JLD2.jldsave(joinpath(root,"template.jld2");record)
    hashes["template.jld2"]=_template_sha(read(joinpath(root,"template.jld2")))
    write(joinpath(root,"checksums.tsv"),join((hashes[name]*"\t"*name for name in sort(collect(keys(hashes)))),"\n")*"\n")
    return root
end

"""
    load_noiseless_template(directory) -> NoiselessTemplate

Verify and load a complete relocated bundle without TEMPO2 or its original
workspace. Reject missing/extra/corrupt files, symlinks, unsafe paths and
inconsistent final scientific metadata. Only the complete current format (marker
3) is supported; older, unknown or incomplete records are rejected without
migration or inferred fields. Original observation identities remain
in `lineage`; re-reading the saved TIM supplies new document IDs.
"""
function load_noiseless_template(directory::AbstractString)
    root=abspath(directory)
    isdir(root) && !islink(root) || throw(ArgumentError("Expected a regular bundle directory"))
    actual=Set{String}()
    for (dir,dirs,names) in walkdir(root;follow_symlinks=false)
        for name in vcat(dirs,names)
            islink(joinpath(dir,name)) && throw(ArgumentError("Bundle contains a symlink"))
        end
        for name in names;push!(actual,relpath(joinpath(dir,name),root));end
    end
    manifest=joinpath(root,"checksums.tsv")
    isfile(manifest) || throw(ArgumentError("Incomplete bundle: checksums.tsv missing"))
    hashes=Dict{String,String}()
    for line in eachline(manifest)
        columns=split(line,'\t')
        _template_require(length(columns)==2 && occursin(r"^[0-9a-f]{64}$",columns[1]),"Invalid checksum manifest")
        name=_template_safe_relative(columns[2])
        _template_require(name!="checksums.tsv" && !haskey(hashes,name),"Duplicate checksum entry")
        hashes[name]=columns[1]
    end
    _template_require(actual==union(Set(keys(hashes)),Set(["checksums.tsv"])),"Bundle files missing or unexpected")
    for (name,hash) in hashes
        _template_require(_template_sha(read(joinpath(root,name)))==hash,"Bundle checksum mismatch: $name")
    end
    try
        record=JLD2.load(joinpath(root,"template.jld2"),"record")
        _template_require_fields(record,("version",),"bundle record")
        _template_require(record["version"] isa Integer && !(record["version"] isa Bool) &&
            record["version"]==_NOISELESS_TEMPLATE_FORMAT,"Unsupported bundle format: only complete format 3 is supported")
        _template_require_fields(record,("metadata","lineage"),"bundle record")
        _template_require_fields(record["metadata"],("format_version",),"template metadata")
        m=Dict{String,Any}(record["metadata"])
        _template_require(m["format_version"]==record["version"],"Bundle record/metadata versions differ")
        artifacts=Dict(name=>read(joinpath(root,name)) for name in keys(hashes) if name!="template.jld2")
        settings=_template_settings_from_record(m)
        doc=read_tim_document(joinpath(root,"final.tim");precision_bits=settings.precision_bits)
        lineage=[TOALineage(x["output_index"],TimTOAID(x["source_sha256"],x["source_line"],Symbol(x["source_origin"])),x["source_ordinal"]) for x in record["lineage"]]
        return _validate_noiseless_template(NoiselessTemplate(doc,lineage,m,artifacts))
    catch err
        throw(ArgumentError("Invalid template bundle: $(sprint(showerror,err))"))
    end
end
