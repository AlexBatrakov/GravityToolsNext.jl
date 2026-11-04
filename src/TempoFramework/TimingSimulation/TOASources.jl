"""Origin-to-output mapping; source ordinal and ID survive subsequent selections."""
struct TOALineage
    output_index::Int
    source_id::TimTOAID
    source_ordinal::Int
end

"""
    ObservationSource

Owned `document::TimDocument`, ordered `lineage::Vector{TOALineage}`, and saved
`settings::NamedTuple`. Construct via `tim_source`, `explicit_schedule`, or
`regular_schedule`; select with `select_toas`. None of these operations writes
files or evaluates a timing model. Persist lineage separately from a written TIM.
"""
struct ObservationSource
    document::TimDocument
    lineage::Vector{TOALineage}
    settings::NamedTuple
end

function _observation_source(doc, settings)
    entries = toa_entries(doc)
    isempty(entries) && throw(ArgumentError("observation source must contain at least one TOA"))
    lineage = [TOALineage(i,e.id,e.ordinal) for (i,e) in enumerate(entries)]
    return ObservationSource(doc, lineage, deepcopy(settings))
end

"""
    tim_source(document::TimDocument)
    tim_source(path; precision_bits=128)

Create an independent, nonempty existing-TIM source. Preserve current edited
epochs, metadata, duplicate observations, raw document lines and original IDs.
Writing its document retains original-file alias protection.
"""
function tim_source(doc::TimDocument)
    owned = subset_tim_document(doc, [e.id for e in toa_entries(doc)])
    return _observation_source(owned, (kind=:tim, origin=doc.origin,
        source_path=doc.source_path, source_sha256=doc.source_sha256,
        precision_bits=doc.precision_bits, selections=NamedTuple[]))
end
tim_source(path::AbstractString; precision_bits::Integer=128) =
    tim_source(read_tim_document(path; precision_bits=precision_bits))

function _source_expand(value, n, label)
    if value isa AbstractVector
        length(value) == n || throw(ArgumentError("$label requires exactly $n values"))
        return Any[x for x in value]
    end
    return Any[value for _ in 1:n]
end

function _source_flags(value, n)
    value isa AbstractVector || throw(ArgumentError("flags must be a vector of pairs or per-TOA pair vectors"))
    if all(x -> x isa Pair, value)
        return [deepcopy(value) for _ in 1:n]
    end
    length(value) == n && all(x -> x isa AbstractVector && all(p -> p isa Pair,x), value) ||
        throw(ArgumentError("per-TOA flags require exactly $n vectors of pairs"))
    return deepcopy(value)
end

function _source_rows(epochs; observatory=nothing, frequency_mhz=nothing, uncertainty_us=nothing,
                      backend=nothing, flags=Pair{String,String}[], filenames=nothing, systems=nothing)
    n=length(epochs)
    n>0 || throw(ArgumentError("schedule must contain at least one epoch"))
    backends=_source_expand(backend,n,"backend")
    extra_flags=_source_flags(flags,n)
    names=filenames === nothing ? ["obs_"*lpad(i,6,'0') for i in 1:n] : _source_expand(filenames,n,"filenames")
    sites=_source_expand(observatory,n,"observatory")
    freqs=_source_expand(frequency_mhz,n,"frequency_mhz")
    errors=_source_expand(uncertainty_us,n,"uncertainty_us")
    if systems !== nothing
        systems isa AbstractDict || throw(ArgumentError("systems must map backend strings to observing metadata"))
        all(x -> x === nothing,(observatory,frequency_mhz,uncertainty_us)) ||
            throw(ArgumentError("systems cannot be mixed with direct site/frequency/uncertainty inputs"))
    end
    rows=NamedTuple[]
    for i in 1:n
        rowflags=Pair{String,String}[]
        if systems !== nothing
            be=backends[i]
            be isa AbstractString && haskey(systems,be) || throw(ArgumentError("missing or unknown backend system at observation $i"))
            system=systems[be]
            all(k -> hasproperty(system,k),(:observatory,:frequency_mhz,:uncertainty_us)) ||
                throw(ArgumentError("system $be requires observatory, frequency_mhz, uncertainty_us"))
            sites[i]=system.observatory; freqs[i]=system.frequency_mhz; errors[i]=system.uncertainty_us
            for p in get(system,:flags,Pair{String,String}[])
                p isa Pair || throw(ArgumentError("system flags must be pairs"))
                push!(rowflags,string(first(p))=>string(last(p)))
            end
        end
        all(x -> x !== nothing,(sites[i],freqs[i],errors[i])) ||
            throw(ArgumentError("observatory, frequency_mhz and uncertainty_us are required"))
        append!(rowflags,[string(first(p))=>string(last(p)) for p in extra_flags[i]])
        if backends[i] !== nothing
            backends[i] isa AbstractString || throw(ArgumentError("backend must be a string"))
            any(p -> first(p)=="-be",rowflags) && throw(ArgumentError("backend and an explicit -be flag overlap"))
            push!(rowflags,"-be"=>String(backends[i]))
        end
        push!(rows,(epoch_mjd=epochs[i],filename=names[i],observatory=sites[i],
            frequency_mhz=freqs[i],uncertainty_us=errors[i],flags=rowflags))
    end
    return rows
end

function _schedule_source(epochs; precision_bits=128, mode=nothing, kwargs...)
    rows=_source_rows(collect(epochs);kwargs...)
    doc=create_tim_document(rows;precision_bits=precision_bits,mode=mode)
    saved=[(epoch_mjd=e.tokens[3],filename=e.filename,frequency_mhz=e.tokens[2],
        uncertainty_us=e.tokens[4],observatory=e.observatory,flags=copy(e.flags)) for e in toa_entries(doc)]
    return _observation_source(doc,(kind=:explicit_schedule,origin=:generated,
        schedule_id=doc.source_sha256,precision_bits=doc.precision_bits,mode=mode,
        observations=saved,systems=deepcopy(get(kwargs,:systems,nothing)),selections=NamedTuple[]))
end

"""
    explicit_schedule(epochs; observatory, frequency_mhz, uncertainty_us,
                      backend=nothing, flags=[], filenames=nothing,
                      systems=nothing, precision_bits=128, mode=nothing)

Create a deterministic generated-origin source in supplied epoch order, retaining
duplicates. Epochs are plain decimal strings, >=128-bit BigFloat or exact integer
days. Metadata accepts scalars or vectors of exactly matching length. `flags` is
a shared ordered pair vector or one pair vector per TOA. Optional `systems` maps
backend strings to named tuples containing site/frequency/uncertainty and optional
flags; it replaces direct site/frequency/uncertainty inputs. A supplied backend
adds one `-be` pair and must not overlap explicit `-be` flags. No RNG, hidden write,
pulse numbering or model evaluation is used. Saved settings contain resolved
per-observation metadata and decimal epoch strings.
"""
explicit_schedule(epochs::AbstractVector; kwargs...) = _schedule_source(epochs;kwargs...)

"""
    regular_schedule(start_mjd; count, cadence_s, precision_bits=128, kwargs...)

Generate exactly `count > 0` observations at `start + (i-1)*cadence_s/86400`.
Cadence is positive finite site seconds; absolute MJD never passes through
Float64. Compute each epoch from the anchor at explicit precision, avoiding
accumulated cadence rounding. Other metadata keywords match `explicit_schedule`.
This constructs a schedule, not timing-model-aligned or noiseless arrivals.
"""
function regular_schedule(start_mjd; count::Integer, cadence_s::Real, precision_bits::Integer=128, kwargs...)
    count>0 && !(count isa Bool) || throw(ArgumentError("count must be a positive integer"))
    bits=_tim_bits(precision_bits)
    start=_tim_generated_epoch(start_mjd,bits)
    cadence,epochs=_tim_precision(bits) do
        step=BigFloat(cadence_s)
        isfinite(step) && step>0 && !(cadence_s isa Bool) || throw(ArgumentError("cadence_s must be finite and positive"))
        step,[start+BigFloat(i-1)*step/BigFloat(86400) for i in 1:Int(count)]
    end
    source=_schedule_source(epochs;precision_bits=bits,kwargs...)
    settings=merge(source.settings,(kind=:regular_schedule,start_mjd=_tim_generated_token(start,bits),
        count=Int(count),cadence_s=string(cadence)))
    return ObservationSource(source.document,source.lineage,settings)
end
