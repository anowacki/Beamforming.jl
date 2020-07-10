
"""
    beamform(s::AbstractVector{Seis.Trace}, t1, t2, sx1, sx2, sy1, sy2, ds; maxima=1, method=:linear, n=1) -> ::BeamformGrid

Compute normalised beam power at each horizontal slowness point
`sx1:ds:sx2` × `sy1:ds:sy2` s/°, between times `t1` and `t2` s.  The slowness
grid is spaced by `ds` s/°.

The `maxima` maximum beam power points are set in the returned `BeamformGrid`.

    beamform(s, t1, t2, smax, ds; maxima=1, method=:linear, n=1) -> ::BeamformGrid

Use a slowness grid of `-smax:ds:smax` in both dimensions.

# Keyword arguments

- `maxima=1`: The top `maxima` points with the highest beam power are saved to the
  `BeamformGrid` returned.
- `method=:linear`: Set the stacking method to use.  See [`stack`](@ref) for available
  stacking methods.
- `n`: Power of the nth-root or phase-weighted stack.
- `wavefront=:plane`: Geometry of problem to consider.  If `:plane`, consider
  a plane wave moving across the array.  If `circle`, assumes a circular wavefront
  with backazimuth and absolute slowness determined by the slowness grid points.
"""
function beamform(s::TraceArray{T}, t1, t2, sx1, sx2, sy1, sy2, ds;
        maxima=1, method=:linear, n=nothing, wavefront=:plane) where T
    sx = s_per_km.(sx1:ds:sx2)
    sy = s_per_km.(sy1:ds:sy2)
    nx, ny = length.((sx, sy))
    beam_power = Array{T}(undef, nx, ny)
    lon, lat = s.sta.lon, s.sta.lat
    evt = first(s).evt
    x, y, z, mean_lon, mean_lat = array_geometry(lon, lat, zeros(length(s)))
    Δ̄ = wavefront === :circle ? delta(mean_lon, mean_lat, evt.lon, evt.lat) : nothing
    for (iy, isy) in enumerate(sy), (ix, isx) in enumerate(sx)
        align = shift_times(x, y, isx, isy, wavefront;
            lon=lon, lat=lat, mean_lon=mean_lon, mean_lat=mean_lat, mean_dist=Δ̄)
        S = stack(s, (t1,t2), align; method=method, n=n)
        beam_power[ix,iy] = sum(x->x^2, Seis.trace(S))
    # end
    end
    sx, sy = s_per_degree.(sx), s_per_degree.(sy)
    sx_max, sy_max, p, β = find_maxima(sx, sy, beam_power, maxima)
    BeamformGrid{T,eltype(s)}(mean_lon, mean_lat, x, y,
        collect(sx1:ds:sx2), collect(sy1:ds:sy2),
        beam_power./maximum(beam_power), maxima,
        sx_max, sy_max, p, β, s, method, n, t1, t2)
end
beamform(s::TraceArray, t1, t2, smax, ds; kwargs...) =
    beamform(s, t1, t2, -smax, smax, -smax, smax, ds; kwargs...)

"""
    vespagram(traces::AbstractVector{<:Seis.Trace}, t1, t2, s1, s2, ds; kwargs...) -> vesp::VespaGrid

Create a vespagram for a set of `traces`, between times `t1` and `t2` seconds,
and for the range of slownesses from `s1` to `s2` s/° in steps of `ds`
s/°.

# Keyword arguments
- `maxima=1`: The top `maxima` points with the highest beam power are saved to the
  `VespaGrid` returned.
- `method=:linear`: Set the stacking method to use.  See [`stack`](@ref) for available
  stacking methods.
- `n`: Power of the nth-root or phase-weighted stack.
- `wavefront=:plane`: Geometry of problem to consider.  If `:plane`, consider
  a plane wave moving across the array.  If `circle`, assumes a circular wavefront
  with backazimuth and absolute slowness determined by the slowness grid points.
- `envelope=false`: Whether or not to stack the envelopes of traces
"""
function vespagram(s::TraceArray{T}, t1, t2, s1, s2, ds;
        maxima=1, method=:linear, n=nothing, envelope=false, wavefront=:plane) where T
    slow = s_per_km.(s1:ds:s2)
    t = t1:first(s).delta:t2
    nslow = length(slow)
    nt = length(t)
    vespa = Array{T}(undef, nt, nslow)
    lon, lat = s.sta.lon, s.sta.lat
    evt = first(s).evt
    x, y, z, mean_lon, mean_lat = array_geometry(lon, lat, zeros(length(s)))
    Δ̄ = wavefront === :circle ? delta(mean_lon, mean_lat, evt.lon, evt.lat) : nothing
    ᾱ = Geodesics.azimuth(mean_lon, mean_lat, evt.lon, evt.lat) #+ 180
    for (i, islow) in enumerate(slow)
        sx, sy = islow.*sincos(deg2rad(ᾱ))
        align = shift_times(x, y, sx, sy, wavefront; lon=lon, lat=lat,
            mean_lon=mean_lon, mean_lat=mean_lat, mean_dist=Δ̄)
        S = stack(s, (t1, t2), align; method=method, n=n)
        envelope && (Seis.envelope!(S))
        vespa[:,i] .= Seis.trace(S)
    end
    slow = s_per_degree(slow)
    t_max, slow_max, _, _ = find_maxima(t, slow, vespa, maxima)
    VespaGrid{T,eltype(s)}(mean_lon, mean_lat, x, y, t, slow, vespa,
        maxima, t_max, slow_max, s, method, n, wavefront, envelope)
end

"""
    shift_times(x, y, sx, sy, wavefront=:plane; kwargs...) -> shifts::Vector

Compute a set of shift times in s, `shifts`, which represent the relative
time at which a wavefront travelling at x-slowness `sx`  and y-slowness `sy` s/km
occurs at each point in space `x` and `y` km.  By default, a plane wave is assumed;
set `wavefront=:circle` to use a circular wavefront.

# Keyword arguments
## Required if `wavefront == :circle`:
- `lon, lat`: The longitude and latitude (°) of each station.
- `mean_lon, mean_lat`: Reference longitude and latitude (°) of the array
- `mean_dist`: Distance from event to reference point
## Optional if `wavefront == :circle`:
- `earth_radius = $(R_EARTH_KM_DEFAULT)`: Planetary radius in km
"""
function shift_times(x, y, sx, sy, wavefront=:plane;
                     lon=nothing, lat=nothing, mean_lon=nothing, mean_lat=nothing,
                     mean_dist=nothing, earth_radius=R_EARTH_KM_DEFAULT)
    if wavefront === :plane
        -(x.*sx .+ y.*sy)
    elseif wavefront === :circle
        any(x -> x===nothing, (lon, lat, mean_lon, mean_lat, mean_dist)) &&
            throw(ArgumentError("lon, lat, mean_lon, mean_lat and mean_dist " *
                                "must be given when `wavefront` is `:circle`"))
        s = sqrt(sx^2 + sy^2) # s/km
        β = rad2deg(atan(sx, sy))
        # Location of new 'event' at this slowness and backazimuth
        new_evt_lon, new_evt_lat = geodesic_endpoint(mean_lon, mean_lat, β, mean_dist)
        # Distances to new 'event' from each station, in °
        Δs = delta.(new_evt_lon, new_evt_lat, lon, lat)
        # Distance from array reference to new 'event', in °
        Δ = delta(new_evt_lon, new_evt_lat, mean_lon, mean_lat)
        # Station delays are relative times compared to array centre
        (Δs .- Δ).*s_per_degree(s, earth_radius)
    else
        throw(ArgumentError("unknown `wavefront` type"))
    end
end
