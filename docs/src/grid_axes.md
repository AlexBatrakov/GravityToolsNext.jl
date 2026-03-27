# Grid Axes

Define parameter axes with `GridAxis` and a sampling rule.

## Types
- `GridAxis(name; min, max, N, rule=LinRule())` — general constructor
- `LinRule()` — linear spacing
- `LogRule()` / `LogRule(sign)` — logarithmic spacing in `log10(abs(x))`
- `ExplicitRule(vals)` or `GridAxis(name, vals)` — explicit points

`LinAxis`, `LogAxis`, and `ExplicitAxis` are exported type aliases for typed
axes, but the supported constructor surface is `GridAxis(...)`.

Helpers: `linspace`, `axisvalues`, and `refine`.

## Examples
```julia
x = GridAxis(:PX; min=0.0, max=5.0, N=51, rule=LinRule())
y = GridAxis(:PY; min=1e-4, max=1e0, N=41, rule=LogRule(+1))
z = GridAxis(:PZ, [0.1, 1.0, 10.0])

xs = axisvalues(x)
ys = axisvalues(y)
```
