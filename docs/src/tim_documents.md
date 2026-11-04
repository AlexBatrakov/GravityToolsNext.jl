# Precision-preserving TIM documents

`read_tim_document` reads a strict, flat FORMAT 1 document for epoch editing.
It preserves untouched bytes and uses explicit BigFloat precision (128 binary
significand bits by default). The existing `read_tim_file` / `TimTOAEntry`
analysis API remains separate and retains its Float64 behavior.

```julia
using GravityToolsNext

doc = read_tim_document("observations.tim"; precision_bits=128)
toa = first(toa_entries(doc))
shift_toa_epoch!(doc, toa.id, 1e-9) # add a 1 ns site delay

# Absolute epochs must be BigFloat, parsed directly from a decimal string.
epoch = BigFloat("55825.123456789012345678901234", RoundNearest; precision=128)
set_toa_epoch!(doc, toa.id, epoch)
write_tim_document(doc, "edited.tim") # output must not already exist
```

TOA snapshots expose `id`, `source_line`, `ordinal`, `epoch_mjd`, `filename`,
`freq_mhz`, `uncertainty_us`, `observatory`, original `tokens` and ordered
`flags::Vector{Pair{String,String}}`. Duplicate flags and identical TOAs remain
distinct. Snapshots are independent copies; use editing functions to update
the document. Positive site delays increase MJD; no residual, barycentric,
noise or pulse-number adjustment is implicit.

Source IDs combine the input content SHA256 and original physical line number.
Edits preserve those IDs in memory. Re-reading a changed file produces new IDs;
retain an explicit old-ID/output-row mapping when persisting derived datasets.
Unchanged content read again has the same IDs regardless of filename.

Generated documents and safe subsetting are supported through
[`create_tim_document`](@ref) and [`subset_tim_document`](@ref). A generated
document has `origin == :generated`, no source path or inode, and content-based
IDs in a separate namespace. Existing file-origin fields/IDs keep their meaning.
See [Observation sources and selection](observation_sources.md).

Supported syntax is one leading `FORMAT 1`, optionally one `MODE 0` or `MODE 1`
before TOAs, blank lines, column-1 `#` comments and column-1 `C` followed by
whitespace/end. TOA rows contain filename, frequency in MHz (finite, >=0),
plain five-digit decimal MJD in `[10000,100000)`, nominal uncertainty in us
(finite, >0), observatory, and zero or more flag/value pairs. Delimiters are
spaces/tabs. LF, CRLF and a missing final newline are preserved. `-pn` values
are exact signed integer tokens, including integers beyond Float64 precision.

Unsupported directives (including INCLUDE, END, SKIP, TIME, PHASE, TIM JUMP
and error/filter directives), inactive `I` rows, legacy formats, quoted values,
inline comments and malformed rows throw `TimFormatError` with path and physical
line. Nothing is silently skipped or followed. Column-1 filenames beginning
`C`/`I` are rejected because TEMPO2 interprets them as comments/inactive rows.
Flag names start with `-` plus an ASCII letter/underscore; subsequent characters
may include letters, digits, `_ . : + -`. Values beginning with `-` must be
numeric to avoid ambiguous missing values. Flag preservation does not validate
their scientific effects; active timing flags need a separate engine policy.

Editing replaces only the epoch token, with at least 24 fractional decimal
digits (retaining a greater original width). Precision is independent of ambient
MPFR precision/rounding and is restored after success or failure. The I/O
criterion is <=1e-18 s against a 256-bit requested epoch, and <=1e-15 s after
1,000 write/read cycles. Binary BigFloat is not decimal-exact or IEEE binary128
storage. Precision scopes are serialized within this API; unrelated concurrent
code changing MPFR global settings must synchronize its own access.

The writer validates the complete candidate before creating output and publishes
a completed temporary sibling. Its parent directory must exist. Existing
regular output files require `overwrite=true`; source paths, saved source
inodes, aliases and output symlinks remain protected even with that option.
Overwriting replaces the output directory entry rather than truncating other
hardlinks. Non-overwriting publication requires filesystem hardlink support;
unsupported publication fails and cleans up the temporary file. The API does
not provide crash-durability guarantees or protection against adversarial
concurrent directory changes. Modify documents only through the supported API.

```@docs
TimDocument
TimTOA
TimTOAID
TimFormatError
read_tim_document
toa_entries
set_toa_epoch!
shift_toa_epoch!
write_tim_document
```
