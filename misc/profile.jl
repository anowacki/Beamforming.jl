import Juno
using SAC
using Revise
using Beamforming
using Plots

Beamforming.verbose(false)

function f()
    s = sample(:array)
    t = cut(s, 1110, 1140)
    sx, sy, fk = Beamforming.fk(t, 1110, 1140, 6, 0.5, 0.1, 2, true)
    sx, sy, fk
end

Juno.@profiler f()
sx, sy, fk = f()

heatmap(sx, sy, fk')
