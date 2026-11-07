# Observation sources and selection

These operations prepare observations entirely in Julia. Epochs have not been
aligned to a timing model. They perform no engine execution, noise draw, fitting
or hidden file write. `ObservationSource` owns a `document`, ordered `lineage`
and saved `settings`; write its document explicitly with `write_tim_document`.

```julia
using GravityToolsNext, Random

source = tim_source("observations.tim"; precision_bits=128)
selected = select_toas(source; method=:random, count=32, rng=Xoshiro(42))
write_tim_document(selected.document, "selected.tim")
mapping = selected.lineage # retain separately when persisting the output

schedule = regular_schedule("55825.123456789012345678901234";
    count=12, cadence_s=600, observatory="pks",
    frequency_mhz=1400, uncertainty_us=2, backend="PDFB1")
```

`explicit_schedule(epochs; ...)` retains supplied order and duplicates.
`regular_schedule(start_mjd; count, cadence_s, ...)` uses exactly `count > 0`
epochs, each computed from the anchor as `start+(i-1)*cadence_s/86400`, with
positive finite cadence in site seconds. Both use explicit >=128-bit precision.
Absolute MJDs accept S1 decimal strings, BigFloats with >=128 bits, or exact
integer days; Float64 absolute MJDs are rejected. Regular cadence can be Float64
seconds. Generation preserves the binary epoch through enough decimal digits
to round-trip it at document precision. It does not invent `-pn` or validate
the physical effects of explicitly supplied flags.

Site/frequency/nominal uncertainty, backend and filenames accept shared scalars
or vectors matching the epoch count. `flags` accepts an ordered shared pair
vector or one pair vector per observation; duplicates are retained. Backend
adds one `-be` pair; combining `backend` with an explicit `-be` is rejected.
Alternatively, pass a backend sequence and a simple system dictionary:

```julia
systems = Dict(
    "A" => (observatory="pks", frequency_mhz=1400, uncertainty_us=2,
            flags=["-fe"=>"receiver1"]),
    "B" => (observatory="@", frequency_mhz=0, uncertainty_us=3),
)
schedule = explicit_schedule(["55825.0", "55825.1", "55825.1"];
    backend=["A", "B", "A"], systems=systems)
```

`systems` replaces direct site/frequency/uncertainty arguments; mixing them is
an error. Unknown/missing systems fail. Additional supplied flags follow system
flags, then the backend pair. Unused system definitions are permitted and saved.
Resolved per-observation metadata, exact decimal epochs and schedule parameters
remain in settings. No cadence-distribution or scheduling language is implied.
Generated filenames beginning `C`, `I` or `#` are rejected explicitly so no
requested observation can turn into a comment or inactive row.

All selection methods retain input order, not chronological order. There is no
replacement or duplicate-ID selection; identical observations remain distinct.

| Method | Count and order policy |
| --- | --- |
| `:all` | Every TOA, unchanged order |
| `:random` | Exactly `count` or `floor(N*fraction)`, using `randperm` then restoring input order |
| `:stratified` | Exactly a total `count`/fraction, proportional largest-remainder group allocation; ties use lexical group labels. Alternatively a complete `group_counts` dictionary names every observed group, including zero-count groups |
| `:stride` | `offset:stride:N`, default offset=1; positive stride and offset within input bounds |
| `:time_window` | Inclusive `time_start`/`time_finish` MJD, at least one bound; current edited epochs are used |

Counts must be positive and not exceed capacity. Fractions lie in `(0,1]`; a
floored zero result fails rather than being rounded up. Fraction multiplication
uses its numeric type; BigFloat products use explicit sufficient precision,
independent of ambient settings. Use an exact rational for boundary-sensitive
count fractions. No implicit count is supplied to random methods. Time bounds
accept decimal strings, >=128-bit BigFloat or integer days, not Float64. Empty
results, reversed windows and method-inapplicable keywords fail clearly.

Stratification defaults to `group_flag="-be"`. Every TOA must have exactly one
such flag: missing or repeated flags, even equal repeated values, are errors.
Unknown, omitted, negative, noninteger or over-capacity explicit group counts
are errors; their total must be positive. Groups are processed in lexical order.
Selection retains all supported header/comment/blank nodes in original relative
order, including comments adjacent to removed TOAs; those comments are treated
as document text, not automatically reassigned observation annotations.

Random methods require exactly one explicit `rng` or nonnegative integer `seed`.
Seed convenience creates a private `Xoshiro`. Supplied RNGs are consumed; repeating
a request from the same initial RNG state reproduces it. Reusing a consumed RNG
does not promise repetition. Settings retain RNG type, initial state copy, seed
when supplied, Julia version, allocation, selected input indices and counts.
Algorithms use `randperm`; results are not promised identical across Julia/RNG
versions. Nonrandom selection and schedules consume no RNG. Keep other random
processes on independent streams. Selection is not train/test splitting.

`subset_tim_document(doc, ids)` is the supported low-level composition API;
it copies in input order and preserves current epochs, raw lines, original
ordinals/IDs and source-file protection. `create_tim_document(rows)` constructs
generated observations without inventing original paths or inodes. Generated
documents have `origin == :generated` and `source_path`, `source_realpath`,
`source_file_id` equal to `nothing`. Their content hash and TOA IDs belong to
the generated namespace. File reads use `:file`; the original two-argument
`TimTOAID(hash, line)` constructor still means file origin.

Lineage gives output index, original source ID and original source ordinal,
including through successive selections. A changed output read from disk has
new file IDs; keep the returned mapping separately. Documents/snapshots/settings
are independently owned. Mutate epochs via the accepted TIM editing API.

See `examples/observation_sources.jl` for an example with explicit output paths.
Existing file-origin sources and genuine explicit/regular schedules can feed the
bounded [noiseless template API](noiseless_templates.md). Generated inputs retain
their origin and selected IDs through phase bootstrap and strict final MODE 1
verification; source construction alone makes no timing-model claim.

```@docs
ObservationSource
TOALineage
tim_source
explicit_schedule
regular_schedule
select_toas
subset_tim_document
create_tim_document
```
