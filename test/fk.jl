using Test
using Statistics: mean

using Seis
using Beamforming


"""
    synthetic_arrival(sx, sy) -> ::Vector{SACtr}

Calculate synthetic traces with an arrival at a set slowness (`sx`, `sy`)
s/°
"""
function synthetic_arrival(sx, sy)
    # Synthetic arrival
    data = exp.(-(-50:0.01:50).^2 ./ 0.01)
    # Grid of stations
    delta = 0.01
    nsta = 20
    s = [Trace(0, delta, data) for _ in 1:nsta]
    s.sta.lon = lon = 5rand(nsta) .- 2.5
    s.sta.lat = lat = 5rand(nsta) .- 2.5
    mlon, mlat = mean(lon), mean(lat)
    # Shift arrivals to correct time
    for (i, ss) in enumerate(s)
        dist = Beamforming.delta(ss.sta.lon, ss.sta.lat, mlon, mlat, true)
        az = Beamforming.azimuth(mlon, mlat, ss.sta.lon, ss.sta.lat, true)
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
        nsta = length(s)
        bf = Beamforming.beamform(s, 40, 60, -2, 2, -2, 2, 0.1)
        @test bf isa Beamforming.BeamformGrid
        @test length(bf.x) == length(bf.y) == nsta
        @test bf.sx == -2:0.1:2
        @test bf.sy == -2:0.1:2
        @test size(bf.power) == (length(bf.sx), length(bf.sy))
        @test bf.nmax == 1
        @test length(bf.bazmax) == length(bf.sx_max) ==
            length(bf.sy_max) == length(bf.pmax)
        @test bf.sx_max[1] ≈ sxi atol=0.05
        @test bf.sy_max[1] ≈ syi atol=0.05
    end
end
