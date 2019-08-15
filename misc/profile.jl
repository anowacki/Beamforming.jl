import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
import Juno
using Revise
using Beamforming, Seis
using Plots

Beamforming.verbose(false)

function f()
    t = Seis.sample_data(:array)
    Beamforming.beamform(t, 1110, 1140, 6, 0.25)
end

Juno.@profiler f()
@time fk = f()

heatmap(fk.sx, fk.sy, fk.power', aspect_ratio=:equal,
    xlim=extrema(fk.sx), ylim=extrema(fk.sy),
    c=Plots.cgrad([:white,:blue,:green,:yellow,:orange,:red]))
