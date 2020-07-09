# Functions for correlation beamforming

"""
    crosscorrelation_beamform(traces::AbstractVector{<:Seis.Trace}, t1, t2, sx1, sx2, sy1, sy2, ds; t=(t2-t1)/2, maxima=1, method=:linear, n=1)
    crosscorrelation_beamform(traces, t1, t2, s_max, ds; kwargs...)

Perform cross-correlation beamforming as in Ruigrok et al. (2017).

In this method, the cross-correlations between all pairs of stations in
`traces` are computed for the data between `t1` and `t2` s, then standard beamforming
is applied to these correlation traces on a horizontal slowness grid.

Stacking is performed over a window `t` seconds either side of the zero-lag
time of the cross-correlations.  Hence the total window length is `2t` s.
If `t` is not set, then it is given a default value of `(t2-t1)/2` s.

Normalised **squared** beam power is computed at each horizontal slowness point
`sx1:ds:sx2` × `sy1:ds:sy2` s/°.  The slowness grid is spaced by `ds` s/°.
If only `s_max` and `ds` are specified, then horizontal slownesses between
`-s_max` and `s_max` s/° are searched over in both directions.

If you already have cross-correlations computed, then use
[`crosscorrelation_beamform_corrs`](@ref) to perform the beamforming.

!!! note
    Only a plane wave is currently supported for use with cross-correlation
    beamforming.

See also: [`beamform`](@ref), [`crosscorrelation_beamform_corrs`](@ref).

# Keyword arguments
- `t`: Cross-correlation stacking window (s).  This is applied symmetrically
  about the zero-lag time, so the window length is really `2t` s long.
  (Default `(t2 - t1)/2` s.)
- `maxima`: Store the `maxima` (default `1`) highest-power peaks in the returned
  `BeamformGrid` object.
- `method`: Stacking method to use (default `:linear`).
- `n`: Power of phaseweighted or nth-root stack if used.
- `wavefront`: Geometry of wave used for stacking.  The default (`:plane`) is
  the only geometry currently supported.

# Example
```
julia> t = Beamforming.Synth.synthetic_arrival(-2, 3, nsta=6, lon=0.1rand(6), lat=0.2rand(6))
6-element Array{Seis.Trace{Float64,Array{Float64,1},Seis.Geographic{Float64}},1}:
 Seis.Trace(.1..: delta=0.01, b=0.0, nsamples=10001)
 Seis.Trace(.2..: delta=0.01, b=0.0, nsamples=10001)
 Seis.Trace(.3..: delta=0.01, b=0.0, nsamples=10001)
 Seis.Trace(.4..: delta=0.01, b=0.0, nsamples=10001)
 Seis.Trace(.5..: delta=0.01, b=0.0, nsamples=10001)
 Seis.Trace(.6..: delta=0.01, b=0.0, nsamples=10001)

julia> bf = crosscorrelation_beamform(t, 40, 60, 4, 0.02, t=1);

julia> bf.sx_max[1], bf.sy_max[1] # Maximum power slowness points in x and y
(-2.0, 3.02)
```

# References
- Ruigrok, E., Gibbons, S., Wapenaar, K., 2017.
  Cross-correlation beamforming.
  Journal of Seismology 21, 495–508.
  https://doi.org/10.1007/s10950-016-9612-6
"""
function crosscorrelation_beamform(s::TraceArray{T}, t1, t2, sx1, sx2, sy1, sy2, ds;
        t=(t2-t1)/2, kwargs...) where T
    x, y, z, mean_lon, mean_lat = array_geometry(s.sta.lon, s.sta.lat,
        coalesce.(s.sta.elev, zero(T)))
    corrs = compute_crosscorrelations(s, t1, t2)
    _do_crosscorrelation_beamforming(s, corrs, t, sx1, sx2, sy1, sy2, ds, x, y, z,
        mean_lon, mean_lat; t1=t1, t2=t2, kwargs...)
end

crosscorrelation_beamform(s, t1, t2, s_max, ds; kwargs...) =
    crosscorrelation_beamform(s, t1, t2, -s_max, s_max, -s_max, s_max, ds; kwargs...)

"""
    crosscorrelation_beamform_corrs(corrs, t, sx1, sx2, sy1, sy2, ds; maxima=1, method=:linear, n=nothing, wavefront=:plane)
    crosscorrelation_beamform_corrs(corrs, t, s_max, ds; kwargs...)

Compute the beam power on a horizontal slowness grid for a set of inter-station
cross-correlations `corrs`.  The zero-lag time should occur at 0 s in the
traces, but they need not be the same length.  In the cross-correlations,
peaks at **positive** lags should correspond to energy arriving at the virtual
receiver **later** than at the virtual source.  For virtual source station A and
virtual receiver station B, this is calculated in Julia via `DSP.xcorr(B, A)`,
and this convention is the same in MATLAB and NumPy.

!!! note
    If your cross-correlations have the opposite sense to the convention
    used here, you can reverse them (e.g., with `reverse!.(trace.(corrs))`)
    before calling this function.

Slowness ranges are specified as for [`crosscorrelation_beamform`](@ref).

Specify `t`, the cross-correlation stacking window length (s).
This is applied symmetrically about the zero-lag time, so the window length is
really `2t` s long.

Unique real stations are identified by the coordinates in each trace's
`sta` and `evt` fields.  It is assumed that no two real stations can be at
exactly the same position in space, and the array reference is calculated
relative to the spherical mean position of all real stations.

Each trace in `corrs` must have the following fields set:
- `sta` must contain the station coordinates of the first of the cross-correlation
  pair (the virtual receiver).
- `evt` must contain the station coordinates of the second of the pair
  (the virtual source).

# Keyword arguments
- `maxima`: Store the `maxima` (default `1`) highest-power peaks in the returned
  `BeamformGrid` object.
- `method`: Stacking method to use (default `:linear`).
- `n`: Power of phaseweighted or nth-root stack if used.
- `wavefront`: Geometry of wave used for stacking.  The default (`:plane`) is
  the only geometry currently supported.

See also [`crosscorrelation_beamform`](@ref).
"""
function crosscorrelation_beamform_corrs(corrs::TraceArray{T}, t, sx1, sx2, sy1, sy2, ds;
        kwargs...) where T
    # Find unique event and receiver coordinates to create array reference.
    # This assumes stations have unique positions, and ignores elevation
    stas_lonlatelev = [(s.lon, s.lat, s.elev) for s in unique(corrs.sta)]
    evts_lonlatelev = [(e.lon, e.lat, -1000e.dep) for e in unique(corrs.evt)]
    lonlatelevs = unique(x -> (x[1], x[2]), vcat(stas_lonlatelev, evts_lonlatelev))
    lons = first.(lonlatelevs)
    lats = getindex.(lonlatelevs, 2)
    elevs = coalesce.(last.(lonlatelevs), zero(T))
    x, y, z, mean_lon, mean_lat = array_geometry(lons, lats, elevs)
    _do_crosscorrelation_beamforming([], corrs, t, sx1, sx2, sy1, sy2, ds, x, y, z,
        mean_lon, mean_lat; kwargs...)
end

crosscorrelation_beamform_corrs(corrs, t, s_max, ds; kwargs...) =
    crosscorrelation_beamform_corrs(corrs, t, -s_max, s_max, -s_max, s_max, ds; kwargs...)

"""
Do the actual work of performing the cross-correlation beamforming;
used by `crosscorrelation_beamform` and `crosscorrelation_beamform_corrs`.

Arguments other than `corrs` are passed into the `BeamformGrid` constructor.

Keyword arguments `t1` and `t2` are stored in the `BeamformGrid`, and
others are used to determine the stacking method and parameters.
"""
function _do_crosscorrelation_beamforming(traces, corrs::TraceArray{T},
        t, sx1, sx2, sy1, sy2, ds, x, y, z, mean_lon, mean_lat;
        t1=nothing, t2=nothing,
        maxima=1, method=:linear, n=nothing, wavefront=:plane) where T
    wavefront === :plane ||
        throw(ArgumentError("wavefront type `$wavefront` not currently supported"))
    all(x -> x.delta == first(corrs).delta, corrs) ||
        throw(ArgumentError("all cross-correlations must have the same sampling interval"))
    sx = s_per_km.(sx1:ds:sx2)
    sy = s_per_km.(sy1:ds:sy2)
    nx, ny = length.((sx, sy))
    beam_power = Array{T}(undef, nx, ny)
    for (iy, isy) in enumerate(sy), (ix, isx) in enumerate(sx)
        align = shift_times_crosscorr(corrs, isx, isy)
        S = stack(corrs, (-t,t), align; method=method, n=n)
        beam_power[ix,iy] = sum(abs2, Seis.trace(S))
    end
    sx, sy = s_per_degree.(sx), s_per_degree.(sy)
    sx_max, sy_max, p, β = find_maxima(sx, sy, beam_power, maxima)
    BeamformGrid{T,eltype(corrs)}(mean_lon, mean_lat, x, y,
        collect(sx1:ds:sx2), collect(sy1:ds:sy2),
        beam_power./maximum(beam_power), maxima,
        sx_max, sy_max, p, β, traces, method, n, t1, t2, :crosscorrelation, corrs, t)
end

"""
    compute_crosscorrelations(s, t1, t2, x, y, z) -> corrs

Create a set of traces `corrs`, where each trace contains the cross-correlation
between a pair of traces in `s`, between `t1` and `t2` s of each trace.
The zero time of the trace corresponds to the zero-lag correlation of the traces.

For each trace in `corrs`, the `evt` field has the station location of the first
trace in the cross-correlation.  Its `evt.id` field is the channel code of the
first station.  The `sta` field is the station of the second station.

An error is thrown if any trace starts after `t1` or ends before `t2` s.
"""
function compute_crosscorrelations(s::TraceArray{T}, t1, t2) where T
    delta = first(s).delta
    all(x -> x.delta == delta, s) ||
        throw(ArgumentError("all traces must have the same sampling interval"))
    maximum(Seis.starttime, s) <= t1 && minimum(Seis.endtime, s) >= t2 ||
        throw(ArgumentError("not all traces have data between $t1 and $t2 s"))
    n = length(s)
    # Trace type passed in
    trace_t = eltype(s)
    # Number of correlation traces
    ncorrs = n*(n - 1)÷2
    corrs = Vector{trace_t}(undef, ncorrs)
    k = 1
    for i in 1:n
        # s1 is the 'virtual source' station
        s1 = Seis.cut(s[i], t1, t2)
        # Transform coordinates into m
        evt = Seis.Event{T}(lon=s1.sta.lon, lat=s1.sta.lat, dep=-s1.sta.elev/1000,
                            id=Seis.channel_code(s1))
        for j in (i+1):n
            i == j && continue
            # s2 is the 'virtual receiver' station
            s2 = Seis.cut(s[j], t1, t2)
            corr = DSP.xcorr(Seis.trace(s2), Seis.trace(s1), padmode=:longest)
            npts_corr = length(corr)
            b = -delta*(npts_corr - 1)/2
            t = trace_t(b, delta, corr)
            t.sta = s2.sta
            t.evt = evt
            corrs[k] = t
            k += 1
        end
    end
    corrs
end

"""
    shift_times_crosscorr(corrs, isx, isy) -> shifts

Return a set of `shifts`, which are delays in seconds for each inter-station
cross-correlation in `corrs`, at an east-west slowness of `isx` s/km and
north-south slowness `isy` s/km.
"""
function shift_times_crosscorr(corrs::TraceArray{T}, isx, isy) where T
    α = atan(isx, isy) # Azimuth of wavevector (radians)
    s = sqrt(isx^2 + isy^2) # Absolute horizontal slowness (s/km)
    shifts = Vector{T}(undef, length(corrs))
    for (i, corr) in enumerate(corrs)
        azimuth = deg2rad(Seis.azimuth(corr))
        distance = Seis.distance_km(corr)
        shifts[i] = -s*distance*cos(azimuth - α)
    end
    shifts
end

"""
    crosscorrelation_array_response(α, Δ, sx1, sx2, sy1, sy2, ds, f1, f2, df, s_deg=false; powf=abs2) -> arf
    crosscorrelation_array_response(α, Δ, s_max, ds, f1, f2, df, s_deg; powf=abs2) -> arf

Return the cross-correlation array response function for a set of inter-station
cross-correlations with inter-station azimuths `α` (°) and inter-station
distances `Δ` (km).

The horizontal slowness grid is determined by setting the lower and upper
east slownesses, `sx1` and `sx2, and similarly for the north slowness with
`sy1` and `sy2`, and the grid spacing `ds`.  All values are in s/km, unless
`s_deg` is `true`.

It is also possible in the common case of a square grid centred on zero
slowness to specify the maximum horizontal slowness `s_max` and the grid
spacing `ds`.

Specify the function to apply to each frequency before summation, `powf`.
This function must take a `Complex` number and return a `Real` number.
The default is `abs2`; `abs` is also common.

Note that this method does not assume that cross-correlations are available
in both directions, so there should be two values of α for each station pair,
180° apart.

The array response object returned `arf` is a
[`Beamforming.ArrayResponseCrosscorrelation`](@ref), which can be plotted with
`Plots.plot(arf)`.  See [`Beamforming.plot`](@ref) for more details.
"""
function crosscorrelation_array_response(α, Δ, sx1, sx2, sy1, sy2, ds, f1, f2, df,
        s_deg=false; powf=abs2
    )
    if s_deg
        sx1, sx2, sy1, sy2, ds = s_per_km.((sx1, sx2, sy1, sy2, ds))
    end
    arf = ArrayResponseCrosscorrelation{float(eltype(α))}(α, Δ, sx1, sx2, sy1, sy2, ds, f1, f2, df)
    _compute_arf_crosscorrelation!(arf, powf)
end

crosscorrelation_array_response(α, Δ, s_max, ds, f1, f2, df, s_deg=false; kwargs...) =
    crosscorrelation_array_response(α, Δ, -s_max, s_max, -s_max, s_max, ds, f1, f2, df, s_deg; kwargs...)

"""
    crosscorrelation_array_response(stations, sx1, sx2, sy1, st2, ds, f1, f2, df, s_deg=false; powf=abs2) -> arf
    crosscorrelation_array_response(stations, s_max, ds, f1, f2, df, s_deg=false; powf=abs2) -> arf

Compute the cross-correlation ARF from a set of `stations` (of type `Seis.Station`)
which define the locations of each array element.
"""
function crosscorrelation_array_response(stations::AbstractArray{<:Seis.GeogStation{T}},
        args...; kwargs...
    ) where T
    x, y, z, lon, lat = array_geometry(stations.lon, stations.lat,
        coalesce.(stations.elev, zero(T)))
    α, Δ = _azimuth_distance_ccbf(x, y)
    crosscorrelation_array_response(α, Δ, args...; kwargs...)
end

function crosscorrelation_array_response(stations::AbstractArray{<:Seis.CartStation{T}},
        args...; kwargs...
    ) where T
    to_km = T(0.001)
    α, Δ = _azimuth_distance_ccbf(to_km*stations.x, to_km*stations.y)
    crosscorrelation_array_response(α, Δ, args...; kwargs...)
end

"""
    crosscorrelation_array_response(traces, sx1, sx2, sy1, sy2, ds, f1, f2, df, s_deg=false; powf=abs2) -> arf
    crosscorrelation_beamform_corrs(traces, s_max, ds, f1, f2, df, s_deg=false; powf=abs2) -> arf

Compute the cross-correlation ARF from a set of `traces` which define the locations
of each array element.
"""
crosscorrelation_array_response(traces::TraceArray, args...; kwargs...) =
    crosscorrelation_array_response(traces.sta, args...; kwargs...)

"""
Convert a set of `x` and `y` Cartesian coordinates into azimuths and distances
for all possible pairs (apart from self-pairs).  Distances are in km, and
angles are in °.
"""
function _azimuth_distance_ccbf(x, y)
    n = length(x)
    length(y) == n ||
        throw(DimensionMismatch("x and y coordinate vectors not the same length"))
    α = [atand(x[j] - x[i], y[j] - y[i]) for i in 1:n for j in 1:n if i != j]
    Δ = [hypot(x[j] - x[i], y[j] - y[i]) for i in 1:n for j in 1:n if i != j]
    α, Δ
end

"""
Do the actual ARF computation.  `powf` is the power function to apply to the summed
beam power at each frequency, defaulting to `abs2`.  
"""
function _compute_arf_crosscorrelation!(arf::ArrayResponseCrosscorrelation{T}, powf) where T
    α = deg2rad.(arf.azimuth)
    @inbounds for (j, slowy) in enumerate(arf.sy), (i, slowx) in enumerate(arf.sx)
        α₀ = atan(slowx, slowy)
        s = sqrt(slowx^2 + slowy^2)
        for freq in arf.freqs
            pow = zero(T)*im
            ω = 2π*freq
            @simd for ipair in 1:length(α)
                pow += cis(ω*arf.distance[ipair]*s*cos(α[ipair] - α₀))
            end
            arf.power[i,j] += powf(pow)
        end
    end
    arf
end
