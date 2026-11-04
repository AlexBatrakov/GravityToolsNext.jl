"""Stable TOA origin: content SHA256, original line, and :file/:generated namespace.
`TimTOAID(sha, line)` retains the original two-argument file-origin constructor.
"""
struct TimTOAID
    source_sha256::String
    source_line::Int
    origin::Symbol
end
TimTOAID(sha::String, line::Int) = TimTOAID(sha, line, :file)

Base.:(==)(a::TimTOAID, b::TimTOAID) =
    a.source_sha256 == b.source_sha256 && a.source_line == b.source_line && a.origin == b.origin
Base.isequal(a::TimTOAID, b::TimTOAID) = a == b
Base.hash(id::TimTOAID, h::UInt) = hash(id.origin, hash(id.source_line, hash(id.source_sha256, h)))

"""A strict FORMAT 1 error, with source path and physical line number."""
struct TimFormatError <: Exception
    path::String
    line::Int
    message::String
end
Base.showerror(io::IO, e::TimFormatError) = print(io, e.path, ':', e.line, ": ", e.message)

struct TimRawLine
    kind::Symbol
    raw::String
end

"""
    TimTOA

Read-only TOA snapshot returned by [`toa_entries`](@ref). `epoch_mjd` is a
BigFloat; `freq_mhz` and `uncertainty_us` are numeric metadata views. `tokens`,
ordered `flags`, `raw`, origin `id`, `source_line`, and `ordinal` retain input
information. Use the document editing functions to change epochs.
"""
struct TimTOA
    id::TimTOAID
    source_line::Int
    ordinal::Int
    raw::String
    epoch_span::UnitRange{Int}
    tokens::Vector{String}
    filename::String
    freq_mhz::Float64
    uncertainty_us::Float64
    observatory::String
    flags::Vector{Pair{String,String}}
    epoch_mjd::BigFloat
    edited::Bool
end

"""
    TimDocument

A strict FORMAT 1 document with immutable source provenance and an explicit
`precision_bits` policy. Construct with [`read_tim_document`](@ref) or
[`create_tim_document`](@ref). Edit via
[`set_toa_epoch!`](@ref) / [`shift_toa_epoch!`](@ref), not by modifying fields.
"""
struct TimDocument
    source_path::Union{Nothing,String}
    source_realpath::Union{Nothing,String}
    source_file_id::Union{Nothing,Tuple{UInt64,UInt64}}
    source_sha256::String
    precision_bits::Int
    nodes::Vector{Union{TimRawLine,TimTOA}}
    origin::Symbol
end

const _TIM_DIRECTIVES = Set(split("FORMAT MODE INCLUDE END SKIP NOSKIP TIME PHASE JUMP " *
    "EFAC EQUAD GLOBAL_EFAC T2EFAC T2EQUAD SIGMA ESET EMAX EFLOOR FMIN FMAX INFO " *
    "PROFILE_DIR HEAD TOAS TRACK DITHER"))
const _TIM_EPOCH = r"^[1-9][0-9]{4}\.[0-9]+$"
const _TIM_FLAG = r"^-[A-Za-z_][A-Za-z0-9_.:+-]*$"
const _TIM_INTEGER = r"^[+-]?[0-9]+$"
const _TIM_PRECISION_LOCK = ReentrantLock()

_tim_error(path, line, message) = throw(TimFormatError(path, line, message))
_tim_file_id(s) = (UInt64(s.device), UInt64(s.inode))

function _tim_bits(bits::Integer)
    bits >= 128 || throw(ArgumentError("TIM precision_bits must be at least 128"))
    return Int(bits)
end

# Scope both settings, including error paths. Serialize our own precision scopes;
# unrelated concurrent code changing global MPFR settings must be synchronized too.
function _tim_precision(f::Function, bits::Int)
    lock(_TIM_PRECISION_LOCK) do
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, RoundNearest) do
                f()
            end
        end
    end
end

function _tim_epoch(token::String, bits::Int, path::String, line::Int)
    occursin(_TIM_EPOCH, token) || _tim_error(path, line,
        "epoch must be plain five-digit decimal MJD in [10000, 100000)")
    epoch = _tim_precision(bits) do
        parse(BigFloat, token)
    end
    10000 <= epoch < 100000 || _tim_error(path, line, "epoch outside [10000, 100000)")
    return epoch
end

function _tim_number(token::String, path::String, line::Int, label::String; positive=false)
    x = tryparse(Float64, token)
    (x !== nothing && isfinite(x) && (positive ? x > 0 : x >= 0)) ||
        _tim_error(path, line, "$label must be finite and $(positive ? "positive" : "nonnegative")")
    return x::Float64
end

function _parse_tim_text(text::String, path::String, canonical::Union{Nothing,String},
                         file_id::Union{Nothing,Tuple{UInt64,UInt64}}, bits::Int; origin::Symbol=:file)
    isvalid(text) || _tim_error(path, 1, "invalid UTF-8 input")
    digest = bytes2hex(sha256(codeunits(text)))
    nodes = Union{TimRawLine,TimTOA}[]
    chunks = split(text, '\n'; keepempty=true)
    format_seen = false
    mode_seen = false
    ordinal = 0
    for (line, chunk) in enumerate(chunks)
        # A trailing LF terminates the previous line, rather than creating a node.
        line == length(chunks) && isempty(chunk) && break
        raw = String(chunk) * (line < length(chunks) ? "\n" : "")
        body = endswith(chunk, '\r') && line < length(chunks) ? chop(chunk) : chunk
        any(c -> iscntrl(c) && c != '\t', body) &&
            _tim_error(path, line, "unsupported control character or bare CR line ending")
        if all(c -> c in (' ', '\t'), body)
            push!(nodes, TimRawLine(:blank, raw))
            continue
        elseif startswith(body, '#') || occursin(r"^C(?:[ \t]|$)", body)
            push!(nodes, TimRawLine(:comment, raw))
            continue
        end
        matches = collect(eachmatch(r"[^ \t]+", body))
        tokens = String[m.match for m in matches]
        key = uppercase(tokens[1])
        if key == "FORMAT"
            !format_seen && ordinal == 0 && tokens == ["FORMAT", "1"] ||
                _tim_error(path, line, "expected exactly one leading FORMAT 1")
            format_seen = true
            push!(nodes, TimRawLine(:format, raw))
            continue
        end
        format_seen || _tim_error(path, line, "FORMAT 1 is required before MODE or TOAs")
        if key == "MODE"
            !mode_seen && ordinal == 0 && length(tokens) == 2 &&
                tokens[1] == "MODE" && tokens[2] in ("0", "1") ||
                _tim_error(path, line, "only one MODE 0/1 before TOAs is supported")
            mode_seen = true
            push!(nodes, TimRawLine(:mode, raw))
            continue
        end
        key in _TIM_DIRECTIVES && _tim_error(path, line, "unsupported TIM directive $(tokens[1])")
        tokens[1] == "I" && _tim_error(path, line, "inactive/deleted I rows are unsupported")
        first(body) in ('C', 'I') && _tim_error(path, line,
            "column-1 filenames starting C/I are ignored or deleted by TEMPO2")
        length(tokens) >= 5 || _tim_error(path, line, "expected five TOA fields; unsupported directive or malformed row")
        any(t -> any(c -> c in ('\"', '\''), t), tokens) &&
            _tim_error(path, line, "quoted tokens are unsupported")
        any(t -> startswith(t, '#'), tokens) && _tim_error(path, line, "inline comments are unsupported")
        iseven(length(tokens)-5) || _tim_error(path, line, "flags require one value each")
        freq = _tim_number(tokens[2], path, line, "frequency (MHz)")
        epoch = _tim_epoch(tokens[3], bits, path, line)
        error_us = _tim_number(tokens[4], path, line, "nominal uncertainty (us)"; positive=true)
        flags = Pair{String,String}[]
        for i in 6:2:length(tokens)
            flag, value = tokens[i], tokens[i+1]
            occursin(_TIM_FLAG, flag) || _tim_error(path, line, "invalid flag $flag")
            if startswith(value, '-') && !occursin(_TIM_INTEGER, value)
                number = tryparse(Float64, value)
                (number !== nothing && isfinite(number)) ||
                    _tim_error(path, line, "missing or ambiguous value for flag $flag")
            end
            flag == "-pn" && !occursin(_TIM_INTEGER, value) &&
                _tim_error(path, line, "-pn requires an exact signed integer token")
            push!(flags, flag => value)
        end
        ordinal += 1
        span = matches[3].offset:(matches[3].offset+ncodeunits(tokens[3])-1)
        push!(nodes, TimTOA(TimTOAID(digest, line, origin), line, ordinal, raw, span, tokens,
            tokens[1], freq, error_us, tokens[5], flags, epoch, false))
    end
    format_seen || _tim_error(path, 1, "missing FORMAT 1")
    return TimDocument(origin == :file ? path : nothing, canonical, file_id, digest, bits, nodes, origin)
end

_tim_label(doc::TimDocument) = something(doc.source_path, "<generated TIM>")

"""
    read_tim_document(path; precision_bits=128) -> TimDocument

Read the strict, flat FORMAT 1 subset without Float64 epoch conversion. Retain
untouched bytes, ordered flag pairs, line endings and TOA origin IDs. Unsupported
directives or malformed rows throw `TimFormatError`; no observations are skipped.
Precision is explicit (at least 128 significand bits), independent of ambient MPFR
precision. Does not follow TIM INCLUDE or open TOA filename tokens.
"""
function read_tim_document(path::AbstractString; precision_bits::Integer=128)
    bits = _tim_bits(precision_bits)
    source = abspath(path)
    canonical = realpath(source)
    text, file_id = open(source, "r") do io
        read(io, String), _tim_file_id(stat(io))
    end
    return _parse_tim_text(text, source, canonical, file_id, bits)
end

"""
    toa_entries(document) -> Vector{TimTOA}

Return independent TOA snapshots in input order, including duplicates. IDs encode
the original input SHA256 and physical line; edits preserve them. Re-reading an
edited file gives new source IDs; retain these snapshots as an explicit lineage
mapping when persisting derived files.
"""
toa_entries(doc::TimDocument) = TimTOA[deepcopy(n) for n in doc.nodes if n isa TimTOA]

function _tim_node_index(doc::TimDocument, id::TimTOAID)
    i = findfirst(n -> n isa TimTOA && n.id == id, doc.nodes)
    i === nothing && throw(ArgumentError("TOA ID does not belong to this TIM document"))
    return i
end

function _tim_replace_epoch!(doc::TimDocument, i::Int, value::BigFloat)
    old = doc.nodes[i]::TimTOA
    isfinite(value) && 10000 <= value < 100000 ||
        _tim_error(_tim_label(doc), old.source_line, "epoch outside [10000, 100000) or nonfinite")
    doc.nodes[i] = TimTOA(old.id, old.source_line, old.ordinal, old.raw,
        old.epoch_span, old.tokens, old.filename, old.freq_mhz, old.uncertainty_us,
        old.observatory, old.flags, value, true)
    return doc
end

"""
    set_toa_epoch!(document, id, epoch::BigFloat) -> document

Set one site MJD, rounded at the document precision. Require a BigFloat with at
least 128 bits; parse decimal epoch strings explicitly before calling. Metadata
and source identity stay unchanged. Invalid edits leave the document unchanged.
"""
function set_toa_epoch!(doc::TimDocument, id::TimTOAID, epoch::BigFloat)
    precision(epoch) >= 128 || throw(ArgumentError("epoch BigFloat must have at least 128 bits"))
    i = _tim_node_index(doc, id)
    value = BigFloat(epoch, RoundNearest; precision=doc.precision_bits)
    return _tim_replace_epoch!(doc, i, value)
end

set_toa_epoch!(::TimDocument, ::TimTOAID, ::Real) =
    throw(ArgumentError("absolute epoch must be a BigFloat parsed from a decimal string"))

"""
    shift_toa_epoch!(document, id, delay_s::Real) -> document

Add a site delay in seconds: `epoch_mjd += delay_s / 86400`, using the document
precision. A negative delay advances the TOA. This performs no barycentric or
residual conversion. Float64 seconds are allowed; absolute MJD is never downcast.
"""
function shift_toa_epoch!(doc::TimDocument, id::TimTOAID, delay_s::Real)
    i = _tim_node_index(doc, id)
    value = _tim_precision(doc.precision_bits) do
        delay = BigFloat(delay_s)
        isfinite(delay) || throw(ArgumentError("site delay must be finite (seconds)"))
        (doc.nodes[i]::TimTOA).epoch_mjd + delay / BigFloat(86400)
    end
    return _tim_replace_epoch!(doc, i, value)
end

function _tim_render(doc::TimDocument)
    _tim_bits(doc.precision_bits)
    io = IOBuffer()
    for node in doc.nodes
        if node isa TimTOA && node.edited
            token = _tim_precision(doc.precision_bits) do
                isfinite(node.epoch_mjd) && 10000 <= node.epoch_mjd < 100000 ||
                    _tim_error(_tim_label(doc), node.source_line, "invalid edited epoch")
                digits = max(24, ncodeunits(node.tokens[3])-6)
                token = Printf.format(Printf.Format("%.$(digits)f"), node.epoch_mjd)
                # A valid value immediately below 100000 must not round out of
                # the supported domain merely because its input token was short.
                if startswith(token, "100000.")
                    digits = max(digits, ceil(Int, doc.precision_bits*log10(2))+2)
                    token = Printf.format(Printf.Format("%.$(digits)f"), node.epoch_mjd)
                end
                token
            end
            print(io, node.raw[1:first(node.epoch_span)-1], token,
                  node.raw[last(node.epoch_span)+1:end])
        else
            print(io, node.raw)
        end
    end
    text = String(take!(io))
    # Validate the entire serialized candidate before touching the filesystem.
    _parse_tim_text(text, _tim_label(doc), doc.source_realpath, doc.source_file_id, doc.precision_bits; origin=doc.origin)
    return text
end

function _tim_output_check(doc::TimDocument, target::String, overwrite::Bool)
    target in (doc.source_path, doc.source_realpath) &&
        throw(ArgumentError("TIM output may not overwrite its source"))
    if ispath(target)
        _tim_file_id(stat(target)) == doc.source_file_id &&
            throw(ArgumentError("TIM output aliases its source (symlink/hardlink)"))
        realpath(target) == doc.source_realpath &&
            throw(ArgumentError("TIM output aliases its source"))
    end
    islink(target) && throw(ArgumentError("TIM output symlinks are not supported"))
    if ispath(target)
        overwrite || throw(ArgumentError("TIM output already exists: $target"))
        isfile(target) || throw(ArgumentError("TIM output must be a regular file"))
    end
    # Resolve parent-directory symlinks even for a nonexistent output file.
    joinpath(realpath(dirname(target)), basename(target)) == doc.source_realpath &&
        throw(ArgumentError("TIM output aliases its source through its parent directory"))
    return nothing
end

"""
    subset_tim_document(document, ids) -> TimDocument

Copy selected TOAs in original document order, retaining every supported raw
comment/header/blank node in its original relative position. Reject repeated or
foreign IDs. Preserve current edited epochs, original IDs/ordinals, precision,
origin and source-file overwrite protection. Empty subsets are representable.
The returned document is independent of the input. This does not rebase pulses.
"""
function subset_tim_document(doc::TimDocument, ids::AbstractVector{TimTOAID})
    wanted = Set(ids)
    length(wanted) == length(ids) || throw(ArgumentError("duplicate TOA IDs in subset request"))
    available = Set(n.id for n in doc.nodes if n isa TimTOA)
    issubset(wanted, available) || throw(ArgumentError("subset contains a foreign TOA ID"))
    nodes = Union{TimRawLine,TimTOA}[deepcopy(n) for n in doc.nodes if !(n isa TimTOA) || n.id in wanted]
    return TimDocument(doc.source_path, doc.source_realpath, doc.source_file_id,
        doc.source_sha256, doc.precision_bits, nodes, doc.origin)
end

function _tim_input_token(value, label)
    token = string(value)
    if isempty(token) || any(isspace, token) || any(iscntrl, token)
        throw(ArgumentError("$label must be a single nonempty token"))
    end
    return token
end

function _tim_generated_epoch(value, bits)
    epoch = if value isa AbstractString
        _tim_epoch(String(value), bits, "<generated TIM>", 1)
    elseif value isa BigFloat
        precision(value) >= 128 || throw(ArgumentError("epoch BigFloat must have at least 128 bits"))
        BigFloat(value, RoundNearest; precision=bits)
    elseif value isa Integer && !(value isa Bool)
        BigFloat(value, RoundNearest; precision=bits)
    else
        throw(ArgumentError("absolute MJD must be a decimal string, BigFloat, or exact integer"))
    end
    isfinite(epoch) && 10000 <= epoch < 100000 || throw(ArgumentError("epoch outside [10000, 100000)"))
    return epoch
end

# Enough fractional digits to round-trip a generated binary epoch, including
# values close to the domain boundary. No absolute-MJD Float64 conversion.
function _tim_generated_token(epoch::BigFloat, bits::Int)
    _tim_precision(bits) do
        digits = max(24, ceil(Int, bits*log10(2))+2)
        Printf.format(Printf.Format("%.$(digits)f"), epoch)
    end
end

"""
    create_tim_document(observations; precision_bits=128, mode=nothing)

Construct a generated-origin FORMAT 1 document without files or hidden writes.
Each observation supplies `epoch_mjd`, `observatory`, `frequency_mhz`, and
`uncertainty_us`; optional `filename` defaults to `obs_000001`, etc., and `flags`
is an ordered vector of pairs. Epochs accept plain decimal strings, >=128-bit
BigFloat or exact integer days, never Float64. Validate the S1 subset. Generated
IDs hash the constructed bytes in the `:generated` namespace; source paths and
inodes are `nothing`. No pulse numbers or phase connection are invented.
"""
function create_tim_document(observations::AbstractVector; precision_bits::Integer=128, mode=nothing)
    bits = _tim_bits(precision_bits)
    mode === nothing || (mode isa Integer && !(mode isa Bool) && mode in (0,1)) ||
        throw(ArgumentError("mode must be nothing, 0 or 1"))
    io = IOBuffer()
    print(io, "FORMAT 1\n")
    mode === nothing || print(io, "MODE ", mode, '\n')
    for (i, row) in enumerate(observations)
        all(k -> hasproperty(row,k),(:epoch_mjd,:observatory,:frequency_mhz,:uncertainty_us)) ||
            throw(ArgumentError("observation $i requires epoch_mjd, observatory, frequency_mhz, uncertainty_us"))
        epoch = _tim_generated_epoch(row.epoch_mjd, bits)
        name = _tim_input_token(get(row, :filename, "obs_"*lpad(i,6,'0')), "filename")
        first(name) in ('C','I','#') && throw(ArgumentError("generated filename must not start with C/I/# (comment/inactive syntax)"))
        fields = [name, _tim_input_token(row.frequency_mhz, "frequency"),
            _tim_generated_token(epoch,bits), _tim_input_token(row.uncertainty_us,"uncertainty"),
            _tim_input_token(row.observatory,"observatory")]
        row_flags = get(row,:flags,Pair{String,String}[])
        row_flags isa AbstractVector || throw(ArgumentError("flags must be an ordered vector of pairs"))
        for pair in row_flags
            pair isa Pair || throw(ArgumentError("flags must contain ordered flag/value pairs"))
            append!(fields, [_tim_input_token(first(pair),"flag"), _tim_input_token(last(pair),"flag value")])
        end
        print(io, join(fields,' '), '\n')
    end
    return _parse_tim_text(String(take!(io)), "<generated TIM>", nothing, nothing, bits; origin=:generated)
end

"""
    write_tim_document(document, new_path; overwrite=false) -> String

Write unchanged bytes verbatim; for edited TOAs replace only the epoch token
with at least 24 fractional decimal digits (retain a greater original width).
Validate before creating output. Never overwrite the source, its saved inode,
symlink/hardlink aliases, or output symlinks, even with `overwrite=true`.
The parent directory must exist. Publish a completed temporary sibling; existing
regular outputs require explicit overwrite. Filesystem failures clean up the
temporary file. No engine execution or noise processing is performed.
"""
function write_tim_document(doc::TimDocument, path::AbstractString; overwrite::Bool=false)
    text = _tim_render(doc)
    target = abspath(path)
    _tim_output_check(doc, target, overwrite)
    temporary, io = mktemp(dirname(target))
    try
        write(io, text)
        close(io)
        _tim_output_check(doc, target, overwrite)
        if overwrite
            # rename replaces a directory entry; it does not truncate other hardlinks.
            Base.Filesystem.rename(temporary, target)
        else
            # An atomic no-replace publish, including a target created after validation.
            hardlink(temporary, target)
        end
    finally
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
    end
    return target
end
