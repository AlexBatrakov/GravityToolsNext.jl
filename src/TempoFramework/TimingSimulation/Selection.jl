function _selection_count(n, count, fraction)
    (count === nothing) != (fraction === nothing) ||
        throw(ArgumentError("provide exactly one of count or fraction"))
    k = if count !== nothing
        count isa Integer && !(count isa Bool) || throw(ArgumentError("count must be an integer"))
        count
    else
        fraction isa Real && !(fraction isa Bool) && isfinite(fraction) && 0 < fraction <= 1 ||
            throw(ArgumentError("fraction must be in (0,1]"))
        if fraction isa BigFloat
            _tim_precision(max(256,precision(fraction)+64)) do
                floor(Int,BigFloat(n)*fraction)
            end
        else
            floor(Int,n*fraction)
        end
    end
    1 <= k <= n || throw(ArgumentError("selection count must be between 1 and $n (fraction uses floor)"))
    return Int(k)
end

function _selection_groups(entries, flag)
    flag isa AbstractString || throw(ArgumentError("group_flag must be a string"))
    groups=Dict{String,Vector{Int}}()
    for (i,e) in enumerate(entries)
        values=[last(p) for p in e.flags if first(p)==flag]
        length(values)==1 || throw(ArgumentError("TOA $(e.id) requires exactly one $flag flag for stratification"))
        push!(get!(groups,only(values),Int[]),i)
    end
    return groups
end

function _selection_mjd(value,bits)
    value === nothing && return nothing
    if value isa AbstractString
        occursin(r"^[+-]?[0-9]+(?:\.[0-9]+)?$",value) || throw(ArgumentError("time bounds must be plain decimal MJD"))
        return _tim_precision(bits) do
            parse(BigFloat,value)
        end
    end
    value isa BigFloat || (value isa Integer && !(value isa Bool)) ||
        throw(ArgumentError("time bounds must be decimal strings, BigFloat or exact integer MJD"))
    value isa BigFloat && precision(value)<128 && throw(ArgumentError("MJD bound BigFloat must have at least 128 bits"))
    x=BigFloat(value,RoundNearest;precision=bits)
    isfinite(x) || throw(ArgumentError("time bounds must be finite"))
    return x
end

"""
    select_toas(source; method=:all, count=nothing, fraction=nothing,
                rng=nothing, seed=nothing, stride=1, offset=1,
                time_start=nothing, time_finish=nothing,
                group_flag="-be", group_counts=nothing)

Select a nonempty independent source in input order, without replacement.
Methods: `:all`; `:random` (exact count or floor(N*fraction)); `:stratified`
(proportional largest-remainder counts with lexical group ties, or explicit
complete group_counts); `:stride` (offset:stride:N); `:time_window` (inclusive
MJD bounds, at least one bound). Do not combine method-specific settings.

Random methods require exactly one explicit AbstractRNG or nonnegative integer
seed (private Xoshiro). Supplied RNGs are consumed; settings save a copy of the
initial state, type, seed if provided and Julia version. Replay starts from the
same state, not a reused consumed RNG. No default/global RNG is used. Stratified
TOAs require exactly one grouping flag, even if repeated values agree. Explicit
group counts must name every observed group, allow zero per group, fit capacities
and have positive total. Non-random methods consume no RNG.

All supported raw document lines are retained. Original IDs/current epochs and
lineage persist; writing and re-reading gives new file IDs, so persist the mapping
separately. Selection is not a temporal train/test split or model idealization.
"""
function select_toas(source::ObservationSource; method::Symbol=:all,count=nothing,fraction=nothing,
                     rng=nothing,seed=nothing,stride::Integer=1,offset::Integer=1,
                     time_start=nothing,time_finish=nothing,group_flag="-be",group_counts=nothing)
    method in (:all,:random,:stratified,:stride,:time_window) || throw(ArgumentError("unknown selection method $method"))
    entries=toa_entries(source.document)
    n=length(entries)
    n>0 || throw(ArgumentError("cannot select from an empty source"))
    stochastic=method in (:random,:stratified)
    if !stochastic
        all(x -> x === nothing,(count,fraction,rng,seed,group_counts)) ||
            throw(ArgumentError("count/fraction/RNG/group_counts apply only to random selection"))
    end
    method == :stratified || (group_counts === nothing && group_flag=="-be") ||
        throw(ArgumentError("group settings apply only to stratified selection"))
    method == :stride || (stride==1 && offset==1) || throw(ArgumentError("stride/offset require method=:stride"))
    method == :time_window || (time_start === nothing && time_finish === nothing) ||
        throw(ArgumentError("time bounds require method=:time_window"))
    random=nothing
    rng_metadata=nothing
    if stochastic
        (rng === nothing) != (seed === nothing) || throw(ArgumentError("provide exactly one of rng or seed"))
        seed === nothing || (seed isa Integer && !(seed isa Bool) && seed>=0) ||
            throw(ArgumentError("seed must be a nonnegative integer"))
        random=seed === nothing ? rng : Random.Xoshiro(seed)
        random isa AbstractRNG || throw(ArgumentError("rng must be an AbstractRNG"))
        rng_metadata=(type=string(typeof(random)),seed=seed,state_before=copy(random),julia_version=string(VERSION))
    end
    allocation=nothing
    lower=nothing;upper=nothing
    indices = if method == :all
        collect(1:n)
    elseif method == :random
        k=_selection_count(n,count,fraction)
        sort!(randperm(random,n)[1:k])
    elseif method == :stratified
        groups=_selection_groups(entries,group_flag)
        labels=sort!(collect(keys(groups)))
        if group_counts === nothing
            k=_selection_count(n,count,fraction)
            allocation=Dict(g=>div(k*length(groups[g]),n) for g in labels)
            remaining=k-sum(values(allocation))
            ranking=sort(labels;by=g -> (-mod(k*length(groups[g]),n),g))
            for g in ranking[1:remaining]; allocation[g]+=1; end
        else
            count === nothing && fraction === nothing || throw(ArgumentError("group_counts cannot be combined with count/fraction"))
            group_counts isa AbstractDict && Set(keys(group_counts))==Set(labels) ||
                throw(ArgumentError("group_counts must name every observed group, without unknown groups"))
            all(g -> group_counts[g] isa Integer && !(group_counts[g] isa Bool) &&
                0<=group_counts[g]<=length(groups[g]),labels) || throw(ArgumentError("invalid or impossible group count"))
            sum(values(group_counts))>0 || throw(ArgumentError("group counts select no TOAs"))
            allocation=Dict(g=>Int(group_counts[g]) for g in labels)
        end
        chosen=Int[]
        for g in labels
            k=allocation[g]
            k==0 && continue
            append!(chosen,groups[g][randperm(random,length(groups[g]))[1:k]])
        end
        sort!(chosen)
    elseif method == :stride
        stride>0 && offset>0 && offset<=n && !(stride isa Bool) && !(offset isa Bool) ||
            throw(ArgumentError("stride must be positive and offset must be in 1:$n"))
        collect(Int(offset):Int(stride):n)
    else
        (time_start !== nothing || time_finish !== nothing) || throw(ArgumentError("time window needs a bound"))
        bits=source.document.precision_bits
        lower=_selection_mjd(time_start,bits); upper=_selection_mjd(time_finish,bits)
        lower === nothing || upper === nothing || lower<=upper || throw(ArgumentError("reversed time window"))
        findall(e -> (lower === nothing || e.epoch_mjd>=lower) && (upper === nothing || e.epoch_mjd<=upper),entries)
    end
    isempty(indices) && throw(ArgumentError("selection contains no TOAs"))
    document=subset_tim_document(source.document,[entries[i].id for i in indices])
    lineage=[TOALineage(j,source.lineage[i].source_id,source.lineage[i].source_ordinal) for (j,i) in enumerate(indices)]
    saved=(method=method,count=count,fraction=fraction,stride=stride,offset=offset,
        time_start=lower === nothing ? nothing : string(lower),time_finish=upper === nothing ? nothing : string(upper),
        group_flag=method==:stratified ? group_flag : nothing,group_counts=allocation,
        selected_input_indices=copy(indices),selected_count=length(indices),rng=rng_metadata)
    history=deepcopy(source.settings.selections)
    push!(history,saved)
    return ObservationSource(document,lineage,merge(deepcopy(source.settings),(selections=history,)))
end
