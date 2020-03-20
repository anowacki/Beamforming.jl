"""
# Seis.Synth

Create synthetic array data for testing.
"""
module Synth

using Statistics: mean

import Seis
import Geodesics
import ..Beamforming

export synthetic_arrival

"""
    synthetic_arrival(sx, sy; delta=0.01, nsta=20, lon, lat, distance, noise) -> ::Vector{Seis.Trace}

Calculate synthetic traces with an arrival at a set slowness (`sx`, `sy`) s/°.

Optionally specify the number of stations `nsta`, the trace sampling interval
`delta` s and the station coordinates `lon`° and `lat`°.

If `distance` is not `nothing`, then populate the traces' event info with
the coordinates of the event.  Specify the level of Gaussian `noise`
to add to the trace as a value greater than 0.
"""
function synthetic_arrival(sx, sy;
        delta = 0.01,
        nsta = 20,
        lon = 5rand(nsta) .- 2.5,
        lat = 5rand(nsta) .- 2.5,
        distance = nothing,
        noise = 0.0
        )
    # Synthetic arrival
    data = exp.(-(-50:0.01:50).^2 ./ 0.01)
    maxamp = maximum(data)
    # Add noise
    if noise > 0
        data .+= 2maxamp*noise.*randn.()
    end
    # Grid of stations
    s = [Seis.Trace(0, delta, data) for _ in 1:nsta]
    s.sta.lon = lon
    s.sta.lat = lat
    mlon, mlat = mean(lon), mean(lat)

    # Add event if requested
    evt = Seis.Event{eltype(s[1])}()
    if !isnothing(distance)
        β = atand(sx, sy)
        evt.lon, evt.lat, _ = Geodesics.angular_step(mlon, mlat, β,
            distance, true)
        evt.dep = 0
        evt.id = "Imaginary event"
    end
    s.evt = evt

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

end # module
