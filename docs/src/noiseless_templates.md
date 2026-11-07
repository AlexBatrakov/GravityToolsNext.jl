# Reusable noiseless templates

`create_noiseless_template` aligns existing TIM observations or a genuine schedule to a fixed deterministic
reference PAR using ordinary installed TEMPO2 and Julia epoch corrections.
`save_noiseless_template` and `load_noiseless_template` explicitly persist and
reuse the verified result; neither operation runs an engine.

```julia
using GravityToolsNext

source = tim_source("observations.tim"; precision_bits=128)
selected = select_toas(source; method=:stratified, count=25, seed=42)
template = create_noiseless_template(selected, "reference.par";
    workspace="new-idealization-job",
    settings=IdealizationSettings(tolerance_s=1e-9, max_corrections=5))
save_noiseless_template(template, "new-template-bundle")
reused = load_noiseless_template("new-template-bundle")

# Genuine generated input uses the same create/save/load operations.
schedule = regular_schedule("53590.0"; count=16, cadence_s=600,
    observatory="pks", frequency_mhz=1400, uncertainty_us=2, backend="PDFB1")
generated = create_noiseless_template(schedule, "reference.par";
    workspace="new-schedule-job")
save_noiseless_template(generated, "new-schedule-bundle")
```

The workspace and bundle target must not exist, and their parents must exist.
Inputs and existing jobs are never overwritten. Failed idealization throws
`IdealizationError` with a status and retained workspace; an unconverged or
failed engine output never becomes a template. No noise is drawn and no
estimator or analysis task is called.

## Supported local domain

Engine validation covers **DDSTG GR**, file-origin `tim_source` and genuine
`explicit_schedule`/`regular_schedule` inputs,
and 1–128 TOAs. The development reference is Linux ARM64 Julia 1.12.4 with the
prepared TEMPO2 installation. Installed stock **general2** is an output plugin;
no generator or custom compiled adapter is required. The adapter stages short
filenames and validates its output as a residual query, independently of the
existing fit stdout/RMS parser.

Generated schedules use bounded initial phase alignment before strict weighted
refinement. Requested epochs are approximate observing times that move to a
nearby permitted phase branch. Origin and original selected IDs remain generated,
including after relocation. Selection occurs before idealization; its first
retained row becomes the anchor. Arbitrary cadence, absolute pulse labels and
global root selection are not guaranteed. A small initial half-cycle margin is
diagnostic; subsequent pulse reassignment or nonconvergence still fails without
changing physical parameters or nominal uncertainties.

The current policy rejects TRACK, active START/FINISH selection, unvalidated
phase/timing flags, unknown noise/whitening controls, and unsupported PAR syntax.
JUMPs must use `JUMP -flag string_selector value [fit [uncertainty]]` syntax;
numeric selectors and other JUMP forms are not supported. Ordinary fixed
DDSTG parameters are preserved even when their names are unfamiliar. PAR
serialization must retain their numeric values exactly; unsupported precision
or syntax fails explicitly. Repeated parameter names are rejected.

Supported metadata flags are `-bw`, `-gof`, `-tobs`, `-nbin`, `-fe`, `-tmplt`,
`-be`, `-nch`, `-epoch`, `-f`, and inactive integer `-pn`, plus custom flags used
by deterministic JUMP selectors. Unverified active flags such as `-pnadd`,
`-padd`, `-addsat`, and `-to` are rejected regardless of selectors. The bounded
stock-output subset permits at most 32 pairs per observation and tokens up to
80 bytes without the output delimiter `|`. Original comments, metadata tokens,
nominal errors, source order and duplicates are retained.

The reserved `-gtnid` flag is **deliberately retained** in the evaluated final TIM
and bundle. It maps each output row to `template.lineage` and original source
IDs; it does not change timing. Input collisions and use as a JUMP selector
are rejected. Removing it changes the verified artifact bytes and requires a
new explicit verification policy. Document IDs assigned on bundle reload differ
from original input IDs; original identities remain in `lineage`.

## Fixed model, phase, and convergence

The original PAR is kept separately. A controlled reference clears all physical
and JUMP fit flags, including special -1 flags, sets REFPHS MEAN and MODE 1,
and fixes NITS/GAIN to 1. The actual runner-written init is audited before every
launch. The logical TIM pathname in the engine CWD is checked against the
expected candidate and staged archive before launch and after capture; a stale,
replaced or retargeted working file cannot be called the evaluated candidate.
Every MODE 1 query starts a fresh process from that same reference; output PARs
are never chained. Unflagged START/FINISH values are retained in the init, while
the ordinary engine may update extrema within a query. No observation may be
deleted or omitted. DDSTG computed/bookkeeping state is distinguished from
fitted independent physical parameters.

The adapter removes these explicitly classified evaluator noise settings:

- Selected white-noise terms: TNEF/TNEQ/T2EFAC/T2EQUAD and
  TNGlobalEF/TNGlobalEQ/addTNGlobalEQ.
- The standard scalar TNRedAmp/TNRedGam/TNRedC/TNRedFLow/TNRedCorner family and
  TNsubtractRed. Removing only fit flags or TNsubtractRed is insufficient.

It records each transformation, rejects other TN/ECORR, WAVE, IFUNC, RNAMP/RNIDX,
DM/CM reconstruction and constraint families, checks that only `param_ZERO` fits,
and requires zero reconstructed red/DM/chromatic signals. This validates
preparation of the listed red-enabled PAR family, **not** red-noise generation
or spectral normalization.

The correction residual is general2 `{post}` in **seconds**, equivalent to column
3 of ordinary `residuals_1.dat`. It is recalculated after fitting and has its
**unweighted arithmetic mean removed**. The fitted weighted offset is stored by
the engine but not subtracted from this recalculated table. Column 4 / residualtn
has different centering and must not be used for correction.

For residuals `r` and effective errors `sigma`, success requires both
`maximum(abs, r)` and `maximum(abs, r .- mu_w)` to meet the tolerance, where
`mu_w = sum(r ./ sigma.^2) / sum(1 ./ sigma.^2)`. The default is 1 ns.
The printed fitted offset from `param.vals` and the inferred weighted mean are
saved separately with their definitions. This convention does not establish
absolute TZR phase equality; later noise/recovery work must account for common
component removal and resulting residual dependence.

Corrections add `-r` site seconds through the precision-preserving TIM API.
Site and residual delays are approximately, not exactly, equal. Every candidate
is written/re-read before its engine evaluation, including the successful final
candidate. At most five corrections (six queries) are allowed. The default local
step bound is 0.01 cycles, scaled by F0, and engine pulse numbers must remain
unchanged across queries. This conservative guard is not a universal guarantee
for large phase shifts or arbitrary timing models. Exact zero chi-square is not
required.

## Generated schedule initialization

The generated-only bootstrap uses a second reference differing only by MODE 0.
Fit weights are one while nominal TIM errors remain unchanged. This avoids
inverse-variance amplification in the installed engine's initial fit-statistic
guard. The numeric MODE 0 statistic is not a probabilistic reduced chi-square and
does not establish phase connection. Physics/JUMPs remain fixed; only param_ZERO
may fit. Its offset is audited against arithmetic centering, separately from the
MODE 1 weighted-mean audit. Actual design weights establish the effective mode.
A generated TIM MODE directive is removed only from disposable evaluator copies
in both stages; original bytes and the transformation are saved. Existing-TIM
MODE 0 is rejected before launch.

Choose the nearest model pulse branch relative to the **first selected TOA** in
source order. Record its original ID/epoch, initial engine pulse vector and
half-cycle margin in cycles. The finite, recomputable margin is **diagnostic
only**, including values below 0.01; it is not an acceptance threshold. Hold the
anchor epoch fixed during bootstrap, and shift epochs by `-(r_i-r_anchor)` site
seconds. At most `max_bootstrap_corrections=3` corrections/four queries are allowed.
A changed pulse integer aborts: no TRACK, wrapping, rebasing or silent anchor
changes. This does not infer absolute TZR pulse locking.

Enter MODE 1 only when both residual maxima are at most
`min(0.01/F0, minimum(nominal_errors_s)/10)`. Evaluate the **same serialized
candidate** in MODE 1, require unchanged pulses and deterministic residual parity
within `mode_switch_tolerance_s=1e-14`, then apply the ordinary small corrections
and final criteria. Bootstrap/refinement query bounds and MODE parity limits
remain enforced. The separate 0.01-cycle local correction bound and entrance
domain remain unchanged. Existing-TIM failures never silently fall back to
bootstrap. A correction across a native glitch that changes an engine pulse
integer remains an error; no automatic root selection or exclusion window is
provided.

Small measured cases cover a regular equal-error schedule and an irregular
multi-week schedule with unequal errors, two backends and a fixed JUMP. The
formerly margin-screened 16-TOA schedule now converges through the public API
with an initial margin of about 0.0034426 cycles and final strict residual maxima
below 0.073 ns, confirmed by independent MODE 1 replay. This does not guarantee
unlimited cadence, arbitrary models or convergence on every near-boundary input.

## Persistence and provenance

A bundle contains full-precision final/source TIMs, original/controlled PARs,
actual init files, query stdout/stderr and diagnostic artifacts, a versioned
`template.jld2` record, and per-file SHA256 checksums. Decimal epochs are saved
as strings. Settings, actual source/selection IDs, RNG implementation and
version-bound state representation, transformations, residual/offset/pulse
policy, query mapping, convergence history, argv, engine identities and dirty
package-source hashes are retained.

The API supports **one complete current format, with marker 3**, for both
existing TIMs and generated schedules. Its explicit
`initial_margin_policy="diagnostic_only"` records the finite initial margin and
pulse vector without an initial threshold. Origin/source IDs, requested/final
epochs and displacements, anchor, both models, mode-specific histories/weights,
stage limits and MODE parity are retained and validated.

Query records distinguish executed TIM and archive paths/hashes. Every history
entry must name its `executed_input_tim` snapshot; loading checks that snapshot,
the archived input and candidate bytes agree, alongside model relationships,
epoch corrections and final scientific criteria. Required settings and history
fields are never reconstructed from defaults or inferred paths.

Older/unknown markers, record/metadata marker mismatches, unsupported policies
and incomplete records are rejected with an explicit error. This includes older
experimental records that share marker 3 but lack current required fields.
The marker identifies the supported record, not a promise of backward
compatibility. Historical TIM/PAR/provenance files remain evidence; the API
provides no automatic migration or legacy load/resave path. Saving/loading does
not regenerate a schedule or perform bootstrap.

A bundle can move and reload after the old job disappears. Loading checks the
complete file inventory, hashes, lineage, epoch strings and final scientific
metrics, and rejects incomplete bundles, unsafe paths, symlinks and corruption.
Saving/loading never executes an engine. Treat a template's fields as read-only;
saving rejects changes to verified document/artifact bytes.

Resource provenance is explicitly **partial**. `resource_paths` can add hashes
of caller-declared resources, but this is not an engine access trace and never
claims exhaustive clock/EOP/ephemeris/DDSTG resource coverage. Full reproduction
also requires the recorded compatible engine, runtime and external resources.

## Deadlines and verification

Templates require a finite 1–60 second engine deadline through the existing raw
runner. On POSIX the runner owns a new process group, retains partial output,
sends TERM at the deadline and KILL after a bounded one-second grace, and reaps
its direct child. Inherited descendants are terminated; orphan zombie collection
belongs to the OS reaper. Descendants that deliberately escape into new sessions
are outside this process-group contract. Other platforms reject finite deadline
requests rather than silently running without a limit. Default runner calls
with `timeout_s=nothing` retain their existing behavior and return contract.

Default tests use pure code and fake executables. A separate opt-in test requires
`GTN_ENABLE_NOISELESS_INTEGRATION=1` and an existing persistent
`GTN_NOISELESS_BUDGET_DIR`. It reserves the maximum possible queries before
launch, records actual invocations, and refuses reruns that could exceed 24
real calls in that evidence directory. Unresolved reservations stay charged.
`GTN_ENABLE_S3_CLOSURE=1` adds regular/irregular generated cases, the formerly
margin-screened case (still named `ambiguous` for case selection), and independent
final MODE 1 verification. All share the same
budget. `GTN_S3_CASES` selects named cases (`existing,regular,ambiguous,explicit`)
for a bounded continuation; it never resets recorded attempts or changes the
scientific API. A separate `GTN_ENABLE_MARGIN_INTEGRATION=1` gate runs the exact
16-TOA small-margin public create/save/load, independent MODE 1 replay and a
25-TOA existing-TIM regression. It requires a new persistent budget directory and
uses a 16-call cap; leave the other real gates disabled for that packet. Default
tests never launch real engines. Keep the older full-fixture integration gate
disabled.

See `examples/noiseless_template.jl` for a lightweight explicit create/save/load
workflow. Existing task, fit-parser and analysis-result contracts are preserved.

```@docs
IdealizationSettings
IdealizationError
NoiselessTemplate
create_noiseless_template
save_noiseless_template
load_noiseless_template
```
