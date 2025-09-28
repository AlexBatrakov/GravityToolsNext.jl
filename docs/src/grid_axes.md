# Grid Axes

Define parameter axes using grid rules.

## Types
- `LinAxis(name::Symbol, min, max, n)` — linear spacing
- `LogAxis(name::Symbol, min, max, n)` — logarithmic spacing
- `ExplicitAxis(name::Symbol, values::AbstractVector)` — explicit points

Underlying rules: `LinRule`, `LogRule`, `ExplicitRule` with helper functions
like `linspace`, `axisvalues`, and `refine`.

## Examples
```julia
x = LinAxis(:PX, 0.0, 5.0, 51)
y = LogAxis(:PY, 1e-4, 1e0, 41)
z = ExplicitAxis(:PZ, [0.1, 1.0, 10.0])

xs = axisvalues(x)
ys = axisvalues(y)
```
