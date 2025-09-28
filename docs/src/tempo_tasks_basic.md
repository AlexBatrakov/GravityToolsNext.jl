# BasicTempoTask — Single run

The simplest way to execute a TEMPO run is via `BasicTempoTask`, which wraps a `TempoRunSettings`.

## Minimal example
```julia
using GravityToolsNext

s = TempoRunSettings(
  work_dir="/abs/work", par_input="a.par", tim_input="a.tim", par_output="a_out.par", tempo_version=Tempo2(),
  write_output=true, write_residuals=true,
  work_mode=:jobdir, layout=:split
)

basic = BasicTempoTask(s)
res = run_task(basic)
```

## Staging inputs
You can stage `.par`/`.tim` into a directory (e.g., for batch workflows):
```julia
mkpath("/abs/work/staging")
task_stage_inputs!(basic, "/abs/work/staging")
```

## Overrides
Use `copy_with` to derive settings for a variant run:
```julia
s2 = copy_with(s; par_output="variant_out.par", override_params_upsert=[TP("DDOT", 5e-19)])
basic2 = BasicTempoTask(s2)
res2 = run_task(basic2)
```
