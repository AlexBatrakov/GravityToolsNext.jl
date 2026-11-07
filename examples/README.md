# Examples

This directory contains lightweight supported examples for the current public
package surface.

Notes:
- Examples are repo-local and avoid machine-specific paths.
- TEMPO/TEMPO2-backed examples are intentionally construction-first; uncomment
  the actual `run_task(...)` lines only when local binaries and runtime data are available.
- Examples use repository-local fixtures and avoid machine-specific data paths.
- For CI-friendly unit tests, see the package test suite under `test/`.

## Observation preparation

`observation_sources.jl OUTPUT_DIR` creates a regular schedule and a stratified
TIM subset with explicit output paths. It runs entirely in Julia and does not
align epochs to a timing model or call TEMPO2. Activate the package environment
when running from an external experiment directory. Existing output files are
not overwritten.

`noiseless_template.jl` creates, saves and reloads a bounded DDSTG GR noiseless
template from an existing TIM/PAR, or a literal 16-TOA schedule with
`--regular reference.par new-job new-bundle`, in a prepared TEMPO2 environment. See
[noiseless templates](../docs/src/noiseless_templates.md) for support boundaries.
