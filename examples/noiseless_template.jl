using GravityToolsNext

# Run in an environment with ordinary TEMPO2 and its installed general2 output.
# Supply an existing pulse-local TIM or --regular for a literal small schedule.
length(ARGS)==4 || error("usage: noiseless_template.jl input.tim|--regular reference.par new-job new-bundle")
tim,par,job,bundle=ARGS
source=tim=="--regular" ? regular_schedule("53590.0";count=16,cadence_s=600,
    observatory="pks",frequency_mhz=1400,uncertainty_us=2,backend="PDFB1") : tim_source(tim;precision_bits=128)
selected=tim=="--regular" ? source : select_toas(source;method=:stratified,count=min(25,length(source.lineage)),seed=42)
template=create_noiseless_template(selected,par;workspace=job,settings=IdealizationSettings(max_corrections=3))
save_noiseless_template(template,bundle)
reloaded=load_noiseless_template(bundle) # no engine call
println("Maximum verified residual (ns): ",last(reloaded.metadata["history"])["max_abs_s"]*1e9)
println("Observations: ",length(reloaded.lineage))
println("Source origin: ",reloaded.metadata["source_origin"])
