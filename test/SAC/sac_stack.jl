# Script to compare the stacking reoutines in SAC and those
# in Beamforming.jl
#
# Requires Plots in addition to the packages in this project

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Revise
import Beamforming, Seis, Plots, Seis.Plot
using Plots: plot

t = Seis.sample_data(:array)
t.meta.file = joinpath.(@__DIR__, "stack." .* Seis.channel_code.(t))
Seis.write_sac.(t, t.meta.file)

sacdata = Dict{String,Seis.Trace{Float32,Vector{Float32},String}}()
seisdata = deepcopy(sacdata)

t1, t2 = -5, 10

for stack_type in ("linear", "nthroot", "phaseweight"), n in 1:2:5
    stack_file = joinpath(@__DIR__, "$(stack_type)")
    if stack_type == "linear"
        n > 1 && continue
    else
        stack_file *= "-$n"
    end
    stack_file *= ".SHZ"

    open(pipeline(`sac`, stdin), "w") do f
        println(f, "sss")
        for tt in t
            println(f, "as $(tt.meta.file) delay $(-tt.picks.A.time)")
        end
        println(f, """
            tw $t1 $t2
            sumstack type $stack_type $n
            ws $(stack_file)
            """)
    end

    sacdata[stack_file] = Seis.read_sac(stack_file)
    s = Beamforming.stack(t, (t1, t2), :A,
        method=Symbol(stack_type), n=n)
    s.meta.file = stack_file
    seisdata[stack_file] = s
end

sacdata = sort(collect(values(sacdata)), by=x->x.meta.file)
seisdata = sort(collect(values(seisdata)), by=x->x.meta.file)

p = plot(sacdata .|> Seis.normalise, label=last.(splitpath.(sacdata.meta.file)))
for i in 1:length(seisdata)
    s = seisdata[i] |> Seis.normalise
    plot!(p[i], Seis.times(s), Seis.trace(s), line=(:dash, :red))
end
p
