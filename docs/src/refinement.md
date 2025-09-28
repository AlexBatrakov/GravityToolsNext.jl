# Refinement

Control how the grid refines based on result metrics.

## Units
- `LocalMinimaUnit(:chi2_marginalized)` — focus on local minima of a metric
- `FullUnit()` — dense refinement everywhere
- `DiffUnit()` / `RelDiffUnit()` — refine where absolute/relative differences are large
- `ContourUnit(levels)` / `DiffContourUnit(levels)` — contour-based refinement

Combine units in `RefinementSettings((units...,))`. The grid engine can save chosen metric names in `params_to_save` for later analysis.

## Example
```julia
ref = RefinementSettings((LocalMinimaUnit(:chi2_marginalized), ContourUnit([2.3, 6.17])))
```
