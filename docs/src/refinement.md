# Refinement

Control how the grid refines based on result metrics.

## Units
- `LocalMinimaUnit(:chi2_marginalized; ...)` — focus on local minima of a metric
- `FullUnit(:chi2_marginalized; ...)` — dense refinement everywhere for a named metric
- `DiffUnit(:chi2_marginalized; diff=..., ...)` / `RelDiffUnit(:chi2_marginalized; rel_diff=..., ...)`
- `ContourUnit(:chi2_marginalized; contours=[...], ...)`
- `DiffContourUnit(:chi2_marginalized; diffs=[...], contours=[...], ...)`

Combine units as positional arguments to `RefinementSettings(...)`. The grid
engine can save chosen metric names in `params_to_save` for later analysis.

## Example
```julia
ref = RefinementSettings(
    LocalMinimaUnit(:chi2_marginalized; max=20.0, max_diff=0.1),
    ContourUnit(:chi2_marginalized; contours=[2.3, 6.17]),
    desired_refinement_level = 1,
    params_to_save = (:chi2_marginalized,)
)
```
