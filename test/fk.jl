using Base.Test

using SAC
import SphericalGeom
using Beamforming

"""
    synthetic_arrival(sx, sy) -> ::Vector{SACtr}

Calculate synthetic traces with an arrival at a set slowness (`sx`, `sy`)
s/°
"""
function synthetic_arrival(sx, sy)
    # Synthetic arrival
    data = exp.(-(-50:0.01:50).^2/0.01)
    # Grid of stations
    delta = 0.01
    nsta = 20
    s = [SACtr(data, delta) for _ in 1:nsta]
    s[:stlo] = lon = 5rand(nsta) - 2.5
    s[:stla] = lat = 5rand(nsta) - 2.5
    mlon, mlat = mean(lon), mean(lat)
    # Shift arrivals to correct time
    for (i, ss) in enumerate(s)
        dist = SphericalGeom.delta(ss[:stlo], ss[:stla], mlon, mlat, true)
        az = SphericalGeom.azimuth(mlon, mlat, ss[:stlo], ss[:stla], true)
        dx = dist*sind(az) # °
        dy = dist*cosd(az)
        dt = sx*dx + sy*dy # s
        ss.t = circshift(ss.t, -round(Int, dt/delta))
    end
    s
end

@testset "Form beam" begin
    let sxi = 2rand()-1, syi = 2rand()-1
        s = synthetic_arrival(sxi, syi)
        bf = Beamforming.beamform(s, 40, 60, -2, 2, -2, 2, 0.1)
        @test bf.sx_max[1] ≈ sxi atol=0.05
        @test bf.sy_max[1] ≈ syi atol=0.05
    end
end
