import Juno
using SAC
using Revise
using Beamforming
using Plots

Beamforming.verbose(false)

function f()
    t = s = sample(:array)
    # t = cut(s, 1110-100, 1140+100)
    Beamforming.beamform(t, 1110, 1140, 6, 0.1)
end

Juno.@profiler f()
