# A binary-free residual protocol fixture. This is not a physical timing model.
using Printf
setprecision(BigFloat,256)
mode=popfirst!(ARGS)
mode=="timeout" && (println("partial-output");flush(stdout);sleep(20))
mode=="nonzero" && (println("partial-output");flush(stdout);exit(9))
par=ARGS[findfirst(==("-f"),ARGS)+1]
tim=ARGS[findfirst(==("-f"),ARGS)+2]
weightmode=parse(Int,split(only(filter(l->occursin(r"^MODE\s",l),readlines(par))))[2])
rows=[split(l) for l in readlines(tim) if !isempty(strip(l)) && !startswith(l,"FORMAT") && !startswith(l,"MODE") && !startswith(l,"#") && !startswith(l,"C ")]
n=length(rows)
target=[BigFloat("55825")+BigFloat(i-1)/86400 for i in 1:n]
res=[(parse(BigFloat,r[3])-t)*86400 for (r,t) in zip(rows,target)]
if mode in ["nonconvergence","phase_guard","ambiguity"]
    step=mode=="phase_guard" ? BigFloat("0.1") : mode=="ambiguity" ? BigFloat("0.1") : BigFloat("0.0001")
    res=[BigFloat(i-2)*step for i in 1:n]
end
res.-=sum(res)/n
mode=="parity" && weightmode==1 && (res[1]+=BigFloat("1e-12");res[end]-=BigFloat("1e-12"))
errors=[parse(Float64,r[4])*1e-6 for r in rows]
mu=sum(Float64.(res)./errors.^2)/sum(1 ./ errors.^2)
open("rows.tsv","w") do io
    for (i,r) in enumerate(rows)
        mode=="missing" && i==n && continue
        flags=join(r[6:end]," ")
        mode=="duplicate" && i==n && (flags=replace(flags,r"T\d{6}"=>"T000001"))
        pulse=100+i+(mode=="pulse" && !occursin("candidate000",tim) ? 1 : 0)
        sat=parse(BigFloat,r[3])+(mode=="sat" ? BigFloat("0.0001") : 0)
        err=mode=="error" ? 2*parse(BigFloat,r[4]) : parse(BigFloat,r[4])
        value=mode=="nonfinite" && i==1 ? "NaN" : string(res[i])
        println(io,join([string(sat),string(sat+BigFloat("0.001")),value,string(err),string(pulse),
            mode=="deleted" && i==n ? "1" : "0",mode=="red" ? "0.1" : "0","0","0",flags],'|'))
    end
end
write("param.labels",mode=="fit" ? "0 param_ZERO 0\n1 param_f 0\n" : "0 param_ZERO 0\n")
open("design.matrix","w") do io
    for i in 1:n
        println(io,"$(i-1) 0 1 ",mode=="wrongweight" ? 7.0 : weightmode==0 ? 1.0 : 1/errors[i])
    end
end
write("constraints.matrix",mode=="constraint" ? "0 0 1\n" : "")
write("param.vals",@sprintf("%.6g\n",mode=="offset" ? 0.1 : weightmode==0 ? 0.0 : mu))
mode=="init" && write(par,read(par,String)*"F3 1 0\n")
if mode in ["tim_replace","tim_retarget","archive_edit","candidate_edit"]
    bytes=read(tim,String)
    if mode=="candidate_edit"
        candidate=joinpath(dirname(dirname(pwd())),basename(tim))
        write(candidate,bytes*"C candidate changed by fake\n")
    else
        rm(tim)
        if mode=="tim_retarget"
            write("other.tim",bytes*"C retargeted by fake\n");symlink("other.tim",tim)
        elseif mode=="archive_edit"
            write(tim,bytes)
            write(joinpath(dirname(pwd()),"input",basename(tim)),bytes*"C archive changed by fake\n")
        else
            write(tim,bytes*"C execution input changed by fake\n")
        end
    end
end
println("fake-residual-query")
