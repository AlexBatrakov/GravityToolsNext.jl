const _TEMPLATE_FORMAT=" {FORMAT %.36Lg}{sat}|{bat}|{post}|{err}|{pulseN}|{del}|{tnrn}|{tndm}|{tnchrom}|{flags}\\n"

_template_require(ok,message)=ok || throw(ArgumentError(message))

function _template_tagged_text(doc::TimDocument)
    # Render current edited epochs before appending diagnostic metadata; original
    # comments, headers, metadata tokens and row order otherwise remain intact.
    text=_tim_render(doc); count=0
    lines=split(text,'\n';keepempty=true)
    for i in eachindex(lines)
        line=lines[i]
        (isempty(strip(line)) || startswith(line,"#") || occursin(r"^[Cc](?:\s|$)",line) ||
            occursin(r"^\s*(FORMAT|MODE)\s",line)) && continue
        count+=1
        lines[i]=rstrip(line,['\r'])*" $(_TEMPLATE_ID_FLAG) T"*lpad(count,6,'0')*(endswith(line,"\r") ? "\r" : "")
    end
    _template_require(count==length(toa_entries(doc)),"Diagnostic tagging changed the TOA count")
    return join(lines,'\n')
end

function _template_residual_table(text::String,entries;previous_pulses=nothing)
    _tim_precision(256) do
        ids=["T"*lpad(i,6,'0') for i in eachindex(entries)]
        rows=Dict{String,Any}(); engine_ids=String[]
        for line in split(chomp(text),'\n')
            cols=strip.(split(line,'|';keepempty=true))
            _template_require(length(cols)==10,"Invalid general2 residual columns")
            tokens=split(cols[10])
            _template_require(iseven(length(tokens)),"Malformed returned flags")
            flags=[String(tokens[i])=>String(tokens[i+1]) for i in 1:2:length(tokens)]
            tags=[v for (k,v) in flags if k==_TEMPLATE_ID_FLAG]
            _template_require(length(tags)==1,"Missing/duplicate diagnostic ID")
            id=only(tags)
            _template_require(id in ids && !haskey(rows,id),"Unexpected/duplicate returned ID $id")
            vals=[parse(BigFloat,cols[i]) for i in [1,2,3,4,7,8,9]]
            _template_require(all(isfinite,vals),"Nonfinite residual/coordinate/error/signal")
            _template_require(cols[6]=="0" && all(iszero,vals[5:7]),"Deleted TOA or stochastic reconstruction")
            pulse=string(parse(Int64,cols[5]))
            rows[id]=(sat=vals[1],bat=vals[2],residual=Float64(vals[3]),precise=vals[3],error_s=Float64(vals[4]/1000000),pulse=pulse,flags=flags)
            push!(engine_ids,id)
        end
        _template_require(length(rows)==length(entries),"Missing observations in engine output")
        residuals=Float64[]; precise=BigFloat[]; errors=Float64[]; pulses=String[]; bats=String[]
        for (id,e) in zip(ids,entries)
            row=rows[id]
            _template_require(row.flags==e.flags,"Returned flags differ for $id")
            _template_require(abs(row.sat-BigFloat(e.epoch_mjd))*86400<=BigFloat("1e-18"),"Returned site epoch differs for $id")
            _template_require(isfinite(row.residual) && isfinite(row.error_s) && row.error_s>0,"Invalid residual/error for $id")
            _template_require(isapprox(row.error_s,e.uncertainty_us*1e-6;rtol=1e-12,atol=0),"Effective error changed for $id")
            push!(residuals,row.residual);push!(precise,row.precise);push!(errors,row.error_s);push!(pulses,row.pulse);push!(bats,string(row.bat))
        end
        previous_pulses===nothing || _template_require(pulses==previous_pulses,"Pulse assignment changed between queries")
        _template_require(abs(sum(precise)/length(precise))<=BigFloat("1e-18"),"Residual table is not arithmetic-centered")
        weights=1 ./ errors.^2
        _template_require(all(isfinite,weights) && isfinite(sum(weights)),"Nonfinite residual weights")
        inferred=sum(weights.*residuals)/sum(weights)
        return (residuals_s=residuals,precise_residuals_s=precise,errors_s=errors,pulses=pulses,bat_mjd=bats,
            engine_row_by_output=[findfirst(==(id),engine_ids) for id in ids],
            weighted_mean_s=inferred,max_abs_s=maximum(abs,residuals),
            max_weighted_centered_s=maximum(abs,residuals.-inferred))
    end
end

function _template_fit_audit(labels::String,design::String,constraints::String,offset::String,table;mode=1)
    _template_require(mode in (0,1),"Unsupported evaluator mode")
    _template_require(strip(labels)=="0 param_ZERO 0","Fit contains parameters other than the nuisance offset")
    _template_require(isempty(strip(constraints)),"Unexpected fit constraints")
    lines=filter(!isempty,strip.(split(design,'\n')))
    _template_require(length(lines)==length(table.residuals_s),"Fit count differs from observation count")
    indices=Int[]
    for line in lines
        c=split(line)
        _template_require(length(c)==4 && parse(Int,c[2])==0 && parse(Float64,c[3])==1,"Invalid offset design matrix")
        i=parse(Int,c[1])+1
        _template_require(i in eachindex(table.errors_s),"Invalid fit row index")
        push!(indices,i)
        _template_require(isapprox(parse(Float64,c[4]),mode==0 ? 1.0 : 1/table.errors_s[i];rtol=1e-5),"Unexpected fit weights")
    end
    _template_require(sort(indices)==collect(eachindex(table.errors_s)),"Fit rows missing or duplicated")
    values=split(strip(offset))
    _template_require(length(values)==1,"Invalid offset parameter output")
    printed=parse(Float64,only(values))
    _template_require(isfinite(printed) && isapprox(printed,mode==0 ? 0.0 : table.weighted_mean_s;rtol=1e-5,atol=1e-15),"Fitted offset disagrees with the table convention")
    return String(only(values))
end

function _template_reference_mode(text::String,mode::Int)
    _tim_precision(256) do
        ps=_template_read_par(text)
        pf=TempoParFile("reference.par","",Dict(p.name_symbol=>p for p in ps),[p.name_symbol for p in ps])
        upsert_param!(pf,TP("MODE",mode))
        return join((get_par_file_representation(p) for p in params_as_vector(pf)),"\n")*"\n"
    end
end

_template_evaluator_text(source::ObservationSource)=_template_evaluator_text(source.document)
function _template_evaluator_text(doc::TimDocument;generated=doc.origin==:generated)
    text=_template_tagged_text(doc)
    # A TIM directive overrides PAR weighting. Keep caller/source bytes intact,
    # but remove that directive from generated evaluator copies in both stages.
    if generated
        text=join(filter(l->!occursin(r"^\s*MODE\s",l),split(text,'\n';keepempty=true)),'\n')
    end
    return text
end

function _template_history(result,corrections,mode,stage)
    t=result.table
    return Dict{String,Any}("corrections"=>corrections,"mode"=>mode,"stage"=>stage,
        "query"=>"queries/"*basename(result.job.job_root),
        "input_tim"=>joinpath("queries",basename(result.job.job_root),relpath(result.job.tim_in_path,result.job.job_root)),
        "executed_input_tim"=>joinpath("queries",basename(result.job.job_root),relpath(result.prepared.tim_path,result.job.job_root)),
        "max_abs_s"=>t.max_abs_s,"max_weighted_centered_s"=>t.max_weighted_centered_s,
        "weighted_mean_inferred_s"=>t.weighted_mean_s,"fitted_offset_printed_s"=>result.printed_offset,
        "fitted_offset_source"=>"param.vals: fitted param_ZERO, seconds, printed %g",
        "fit_weight_definition"=>mode==0 ? "unit weights; offset audited against arithmetic mean" : "inverse nominal variance; offset audited against weighted mean",
        "pulse_numbers"=>t.pulses,"bat_mjd"=>t.bat_mjd,"engine_row_by_output"=>t.engine_row_by_output)
end

function _template_correct(candidate,result,next,settings;anchor_fixed=false)
    doc=read_tim_document(candidate;precision_bits=settings.precision_bits)
    before=toa_entries(doc)
    _tim_precision(256) do
        r=result.table.precise_residuals_s
        delays=anchor_fixed ? r.-r[1] : r
        for (e,delay) in zip(before,delays)
            shift_toa_epoch!(doc,e.id,-delay)
        end
    end
    write_tim_document(doc,next)
    reread=read_tim_document(next;precision_bits=settings.precision_bits)
    _tim_precision(256) do
        for (a,b) in zip(toa_entries(doc),toa_entries(reread))
            _template_require(abs(BigFloat(a.epoch_mjd)-BigFloat(b.epoch_mjd))*86400<=BigFloat("1e-18"),"Serialized epoch lost precision")
        end
    end
    anchor_fixed && _template_require(toa_entries(reread)[1].epoch_mjd==before[1].epoch_mjd,"Bootstrap anchor moved")
    return next
end

_template_margin(table,f0)=_tim_precision(256) do
    BigFloat("0.5")-maximum(abs,(table.precise_residuals_s.-table.precise_residuals_s[1]).*BigFloat(f0))
end

function _template_bootstrap(root,candidate,model,source,settings,tempo,history)
    reference=_template_reference_mode(model.text,0)
    write(joinpath(root,"bootstrap.par"),reference)
    pulses=nothing;info=Dict{String,Any}();ready=nothing
    for correction in 0:settings.max_bootstrap_corrections
        result=_template_query(root,candidate,reference,settings,tempo,length(history)+1,pulses;mode=0,stage="bootstrap")
        t=result.table
        push!(history,_template_history(result,correction,0,"bootstrap"))
        if pulses===nothing
            pulses=copy(t.pulses)
            margin=_template_margin(t,model.f0_hz)
            _template_require(isfinite(margin),"Nonfinite initial phase margin")
            info=Dict{String,Any}("anchor_output_index"=>1,"anchor_source_id"=>_template_lineage_record(first(source.lineage)),
                "anchor_epoch_mjd"=>string(first(toa_entries(source.document)).epoch_mjd),
                "initial_pulse_numbers"=>pulses,"initial_margin_cycles"=>string(margin),
                "entrance_limit_s"=>min(0.01/model.f0_hz,minimum(t.errors_s)/10),
                "model_transformation"=>"Set MODE = 0 only; unchanged controlled physical/JUMP model",
                "tim_transformation"=>"Remove TIM MODE directive from generated evaluator copies; preserve source bytes")
        end
        if max(t.max_abs_s,t.max_weighted_centered_s)<=info["entrance_limit_s"]
            ready=result;break
        end
        correction<settings.max_bootstrap_corrections ||
            throw(IdealizationError(:bootstrap_nonconvergence,"Bootstrap correction budget exhausted",root))
        candidate=_template_correct(candidate,result,joinpath(root,"candidate"*lpad(length(history),3,'0')*".tim"),settings;anchor_fixed=true)
    end
    return (candidate=candidate,result=ready,info=info,pulses=pulses)
end

"""
    create_noiseless_template(source, reference_par; workspace,
        settings=IdealizationSettings(), tempo=Tempo2(), resource_paths=String[])

Create a reusable DDSTG GR template from a file-origin `tim_source` or genuine
`explicit_schedule`/`regular_schedule`, up to 128 TOAs. `workspace` is create-only.
Inputs are never edited. Every query uses the existing raw runner, a fresh
process, fixed physical/JUMP parameters and only the common nuisance offset fit.

Generated schedules first use bounded unit-weight MODE 0 initialization relative
to the first selected TOA, preserving its epoch and the initial engine pulse
vector. TIM MODE directives are removed only from generated evaluator copies,
with the transformation recorded. MODE 0 never supplies final acceptance: the
same serialized candidate must pass MODE 1 parity before strict refinement.
The initial half-cycle margin is finite, recorded diagnostics, not a rejection
threshold. The initial pulse vector must still remain exactly fixed afterward.
Both final residual maxima must meet tolerance in MODE 1 on the exact written
TIM. This is a relative phase convention, not absolute TZR pulse locking.

The final TIM retains collision-checked harmless `-gtnid` flags and observing
metadata/errors. Active TRACK, unsupported flags/models/noise controls and
pulse reassignment or nonconvergence are rejected. Failures throw `IdealizationError` and retain
the workspace. Resource provenance is explicitly partial, not an access trace.
No output-PAR chaining, noise generation or implicit cache is used.
"""
function create_noiseless_template(source::ObservationSource,reference_par::AbstractString;
        workspace::AbstractString,settings::IdealizationSettings=IdealizationSettings(),
        tempo::Tempo2=Tempo2(),resource_paths::AbstractVector{<:AbstractString}=String[])
    original=read(reference_par,String)
    model=_prepare_template_par(original)
    _validate_template_source(source,settings,model.selectors)
    root=_template_new_directory(workspace)
    try
        write(joinpath(root,"original.par"),original)
        write(joinpath(root,"reference.par"),model.text)
        write(joinpath(root,"source.tim"),_tim_render(source.document))
        candidate=joinpath(root,"candidate000.tim")
        write(candidate,_template_evaluator_text(source))
        history=Any[]; previous_pulses=nothing; final=nothing;bootstrap=nothing;seed=nothing
        if source.document.origin==:generated
            seed=_template_bootstrap(root,candidate,model,source,settings,tempo,history)
            candidate=seed.candidate;previous_pulses=seed.pulses;bootstrap=seed.info
        end
        for correction in 0:settings.max_corrections
            result=_template_query(root,candidate,model.text,settings,tempo,length(history)+1,previous_pulses)
            t=result.table;previous_pulses=t.pulses
            push!(history,_template_history(result,correction,1,"refinement"))
            if seed!==nothing && correction==0
                parity=_tim_precision(256) do
                    maximum(abs,t.precise_residuals_s.-seed.result.table.precise_residuals_s)
                end
                parity<=settings.mode_switch_tolerance_s ||
                    throw(IdealizationError(:mode_switch_failed,"Deterministic residuals changed at MODE handoff",root))
                bootstrap["mode_switch_max_difference_s"]=string(parity)
                bootstrap["mode0_query"]=history[end-1]["query"]
                bootstrap["mode1_query"]=last(history)["query"]
            end
            maximum(abs,t.residuals_s)*model.f0_hz<=settings.max_step_phase ||
                throw(IdealizationError(:phase_guard,"Residual step exceeds the configured local phase domain",root))
            if t.max_abs_s<=settings.tolerance_s && t.max_weighted_centered_s<=settings.tolerance_s
                final=result;break
            end
            correction<settings.max_corrections || throw(IdealizationError(:nonconvergence,"Correction budget exhausted",root))
            candidate=_template_correct(candidate,result,joinpath(root,"candidate"*lpad(length(history),3,'0')*".tim"),settings)
        end
        doc=read_tim_document(candidate;precision_bits=settings.precision_bits)
        _validate_prepared_tim(final.prepared)
        _template_require(read(final.prepared.tim_path)==read(candidate),"Final serialized candidate differs from evaluated TIM")
        artifacts=Dict{String,Vector{UInt8}}("final.tim"=>read(candidate),"source.tim"=>read(joinpath(root,"source.tim")),
            "original.par"=>Vector{UInt8}(codeunits(original)),"reference.par"=>Vector{UInt8}(codeunits(model.text)))
        bootstrap===nothing || (artifacts["bootstrap.par"]=read(joinpath(root,"bootstrap.par")))
        for h in history
            q=basename(h["query"]);dir=joinpath(root,q)
            for (d,_,names) in walkdir(dir), name in names
                p=joinpath(d,name);artifacts[joinpath("queries",q,relpath(p,dir))]=read(p)
            end
        end
        metadata=_template_provenance(source,reference_par,settings,tempo,resource_paths,model,history,artifacts;bootstrap)
        template=NoiselessTemplate(doc,deepcopy(source.lineage),metadata,artifacts)
        _validate_noiseless_template(template)
        return template
    catch err
        write(joinpath(root,"failure.txt"),sprint(showerror,err))
        err isa IdealizationError && rethrow()
        throw(IdealizationError(:validation_failed,sprint(showerror,err),root))
    end
end

function _template_query(root,candidate,reference,settings,tempo,number,previous_pulses;mode=1,stage="refinement")
    name="q"*lpad(number,3,'0')
    s=TempoRunSettings(work_dir=root,par_input=mode==0 ? "bootstrap.par" : "reference.par",tim_input=basename(candidate),
        par_output="model_out.par",tempo_version=tempo,nits=1,gain=1.0,
        flags="-nofit -writeres -output general2 -file fmt.txt -outfile rows.tsv",
        work_mode=:jobdir,job_name=name,overwrite=:error,layout=:split,link_tim=false,
        cleanup_before_run=false,keep_tmp_on_success=true,keep_tmp_on_error=true,
        write_output=true,write_residuals=true,timeout_s=settings.timeout_s)
    job=materialize_job(s)
    write(joinpath(job.run_cwd,"fmt.txt"),_TEMPLATE_FORMAT)
    prepared=_tim_precision(256) do;_prepare_tempo_run(job);end
    _template_require(read(prepared.init_path,String)==reference,"Actual runner init differs from controlled reference")
    _template_require(read(prepared.tim_path)==read(candidate),"Actual runner TIM differs from candidate")
    _write_json(joinpath(job.job_root,"invocation.json"),Dict("argv"=>collect(prepared.cmd.exec),
        "engine_env"=>tempo_env(tempo),"cwd"=>job.run_cwd,"timeout_s"=>settings.timeout_s,
        "kill_grace_s"=>1,"mode"=>mode,"stage"=>stage,"toa_count"=>length(toa_entries(read_tim_document(candidate;precision_bits=settings.precision_bits))),
        "init_sha256"=>prepared.init_sha256,"tim_sha256"=>prepared.tim_sha256,
        "executed_tim_path"=>prepared.tim_path,"executed_tim_sha256"=>prepared.tim_sha256,
        "archive_tim_path"=>job.tim_in_path,"archive_tim_sha256"=>_template_sha(read(job.tim_in_path))))
    started=time()
    out,err,code=run_tempo_raw(prepared)
    write(joinpath(job.job_root,"stdout.txt"),out);write(joinpath(job.job_root,"stderr.txt"),err)
    diagnostics=Dict{String,Any}(string(k)=>v for (k,v) in prepared.diagnostics)
    diagnostics["duration_s"]=time()-started;diagnostics["exit_code"]=code
    _write_json(joinpath(job.job_root,"execution.json"),diagnostics)
    if code!=0
        status=get(prepared.diagnostics,:timed_out,false) ? :timeout : :engine_failed
        throw(IdealizationError(status,"Engine query failed; partial artifacts are diagnostic only",root))
    end
    _template_require(!occursin(r"(?i)(unrecognized|unrecognised|unknown parameter|not recognised)",out*err),"Engine reports unrecognized model controls")
    _template_require(_template_sha(read(prepared.init_path))==prepared.init_sha256,"Engine input PAR was modified")
    entries=toa_entries(read_tim_document(candidate;precision_bits=settings.precision_bits))
    table=_template_residual_table(read(joinpath(job.run_cwd,"rows.tsv"),String),entries;previous_pulses)
    printed=_template_fit_audit(read(joinpath(job.run_cwd,"param.labels"),String),
        read(joinpath(job.run_cwd,"design.matrix"),String),read(joinpath(job.run_cwd,"constraints.matrix"),String),
        read(joinpath(job.run_cwd,"param.vals"),String),table;mode)
    _validate_prepared_tim(prepared)
    _template_require(read(prepared.tim_path)==read(candidate),"Evaluated TIM differs from candidate")
    write(joinpath(job.job_root,"candidate.tim"),read(candidate))
    return (job=job,table=table,printed_offset=printed,prepared=prepared)
end
