# TempoDataFiles.jl
# Reading TEMPO/TEMPO2 .tim and residual files, and combining them into a single per-TOA view.

# ------------------------------------------------------------
# TIM (.tim) entries
# ------------------------------------------------------------

"""
    TimTOAEntry

Single TOA from a `.tim` file, after removing comments and headers.

Fields:
- `index::Int`          : sequential index after removing comments/headers (1-based)
- `timfile_line::Int`   : original line number in the `.tim` file (1-based)
- `toa::Float64`        : MJD
- `freq::Float64`       : observing frequency (MHz)
- `uncertainty::Float64`: TOA uncertainty (microseconds)
- `backend::String`     : backend label (extracted from `-be <name>` if present, otherwise `"unknown"`)
"""
struct TimTOAEntry
    index::Int
    timfile_line::Int
    toa::Float64
    freq::Float64
    uncertainty::Float64
    backend::String
end

# Internal helper: detect comment/blank lines
_is_comment_or_blank(s::AbstractString) = begin
    t = strip(s)
    isempty(t) || startswith(t, "C") || startswith(t, "c") || startswith(t, "#")
end

# Internal helper: naive header length.
# By default we keep backward-compat with "first 2 non-comment lines are header".
function _header_length(clean_lines::Vector{Tuple{Int,String}})
    # clean_lines: (original_line_number, line_text) after removing comments/blank
    # If there are fewer than 2 lines, use whatever is available.
    return min(2, length(clean_lines))
end

"""
    read_tim_file(path::String) -> Vector{TimTOAEntry}

Reads a `.tim` file, skipping comment/blank lines, preserving original
line numbers, and assuming the first two **non-comment** lines are header/meta.

Notes:
- Expects data lines with (at least) 4 tokens where the 2nd/3rd/4th tokens are
  `freq`, `MJD`, `uncertainty` respectively (compatible with your current format).
- Backend is taken from the `-be <name>` flag if present; otherwise `"unknown"`.
- Unparseable lines are skipped with a warning.
"""
function read_tim_file(path::String)::Vector{TimTOAEntry}
    raw = open(path, "r") do io
        readlines(io)
    end

    # Keep non-comment, non-blank lines with their original indices
    clean = Tuple{Int,String}[]
    for (ln, line) in enumerate(raw)
        _is_comment_or_blank(line) && continue
        push!(clean, (ln, line))
    end

    # Determine header length (kept for backward-compat: 2 non-comment lines)
    hlen = _header_length(clean)
    n_toas = max(length(clean) - hlen, 0)

    entries = Vector{TimTOAEntry}(undef, n_toas)
    out_i = 0

    for i in (hlen+1):length(clean)
        ln, line = clean[i]
        toks = split(strip(line))
        if length(toks) < 4
            @warn "Skipping short .tim line $ln: expected ≥ 4 tokens, got $(length(toks))" line=line
            continue
        end

        # Backend detection: "-be <name>"
        be_idx = findfirst(==("-be"), toks)
        backend = (be_idx === nothing || be_idx == length(toks)) ? "unknown" : String(toks[be_idx + 1])

        # Your current convention: toks[2]=freq, toks[3]=MJD, toks[4]=uncertainty
        # (If your files differ, adjust here or make this branch heuristic.)
        freq        = try parse(Float64, toks[2]) catch; NaN end
        toa         = try parse(Float64, toks[3]) catch; NaN end
        uncertainty = try parse(Float64, toks[4]) catch; NaN end

        if !(isfinite(freq) && isfinite(toa) && isfinite(uncertainty))
            @warn "Skipping unparsable .tim line $ln (freq/mjd/unc not finite)" line=line
            continue
        end

        out_i += 1
        entries[out_i] = TimTOAEntry(out_i, ln, toa, freq, uncertainty, backend)
    end

    # If we skipped any broken lines, shrink the vector
    if out_i != n_toas
        resize!(entries, out_i)
    end

    return entries
end

# ------------------------------------------------------------
# Residuals file (residuals.dat) entries
# ------------------------------------------------------------

"""
    TempoResidualEntry

One line from `residuals.dat` (TEMPO/TEMPO2). Values are converted to microseconds
where appropriate (signal, residual, residual_tn, uncertainty).
"""
struct TempoResidualEntry
    time::Float64
    signal::Float64
    residual::Float64
    residual_tn::Float64
    uncertainty::Float64
end

"""
    read_residual_file(path::String) -> Vector{TempoResidualEntry}

Reads `residuals.dat` (or similar) with 5 columns:
`time  signal  residual  residual_tn  uncertainty`.
Signal, residuals, and uncertainty are converted to microseconds (×1e6).
"""
function read_residual_file(path::String)::Vector{TempoResidualEntry}
    lines = open(path, "r") do io
        readlines(io)
    end

    entries = TempoResidualEntry[]
    for (i, raw) in enumerate(lines)
        s = strip(raw)
        isempty(s) && continue
        fields = split(s)
        if length(fields) != 5
            @warn "Skipping line $i in residuals: expected 5 fields, got $(length(fields))" line=raw
            continue
        end

        try
            time         = parse(Float64, fields[1])
            signal       = parse(Float64, fields[2]) * 1e6
            residual     = parse(Float64, fields[3]) * 1e6
            residual_tn  = parse(Float64, fields[4]) * 1e6
            uncertainty  = parse(Float64, fields[5]) * 1e6
            push!(entries, TempoResidualEntry(time, signal, residual, residual_tn, uncertainty))
        catch err
            @warn "Failed to parse residuals line $i" exception=err line=raw
        end
    end
    return entries
end

# ------------------------------------------------------------
# Combined entries
# ------------------------------------------------------------

"""
    CombinedTOAEntry

Merged view of one TOA with corresponding residual entries.

Fields:
- `index::Int`               : sequential index after header removal (from `.tim`)
- `timfile_line::Int`        : original line number in the `.tim` file
- `in_fit::Bool`             : whether this TOA is within the requested time window
- `toa::Float64`             : MJD
- `freq::Float64`            : frequency (MHz)
- `signal::Float64`          : residuals `signal` column (µs)
- `red_noise::Float64`       : `residual - residual_tn` (µs)
- `residual::Float64`        : residual with red noise (µs)
- `residual_tn::Float64`     : residual with TN removed (µs)
- `uncertainty_orig::Float64`: original TOA uncertainty from `.tim` (µs)
- `uncertainty::Float64`     : transformed uncertainty from residuals file (µs)
- `backend::String`          : backend label
- `weight::Float64`          : `1 / uncertainty^2` (or `0` if uncertainty ≤ 0 or non-finite)
"""
struct CombinedTOAEntry
    index::Int
    timfile_line::Int
    in_fit::Bool
    toa::Float64
    freq::Float64
    signal::Float64
    red_noise::Float64
    residual::Float64
    residual_tn::Float64
    uncertainty_orig::Float64
    uncertainty::Float64
    backend::String
    weight::Float64
end

"""
    combine_tim_and_residuals(tim_entries, residuals; time_start=nothing, time_finish=nothing)
        -> Vector{CombinedTOAEntry}

Combine `.tim` entries and residuals 1:1 (by index). Optional time window (`MJD`) marks
`in_fit` for each entry. Length must match.

Notes:
- `weight = 1 / uncertainty^2` if `uncertainty` is positive and finite; otherwise `0.0`.
- Time window is inclusive: `time_start ≤ toa ≤ time_finish` when both are provided.
"""
function combine_tim_and_residuals(
    tim_entries::Vector{TimTOAEntry},
    residuals::Vector{TempoResidualEntry};
    time_start::Union{Nothing,Float64}=nothing,
    time_finish::Union{Nothing,Float64}=nothing
)::Vector{CombinedTOAEntry}

    N = length(tim_entries)
    length(residuals) == N || error("Mismatched lengths: tim_entries = $N, residuals = $(length(residuals))")

    entries = Vector{CombinedTOAEntry}(undef, N)

    @inbounds for i in 1:N
        toa = tim_entries[i]
        res = residuals[i]

        in_fit = true
        if time_start !== nothing && toa.toa < time_start
            in_fit = false
        end
        if time_finish !== nothing && toa.toa > time_finish
            in_fit = false
        end

        red_noise = res.residual - res.residual_tn

        w = (isfinite(res.uncertainty) && res.uncertainty > 0.0) ? (1.0 / (res.uncertainty^2)) : 0.0

        entries[i] = CombinedTOAEntry(
            toa.index,
            toa.timfile_line,
            in_fit,
            toa.toa,
            toa.freq,
            res.signal,
            red_noise,
            res.residual,
            res.residual_tn,
            toa.uncertainty,
            res.uncertainty,
            toa.backend,
            w
        )
    end

    return entries
end