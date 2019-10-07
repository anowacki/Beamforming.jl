
"""
    beamform(s::AbstractVector{Seis.Trace}, t1, t2, sx1, sx2, sy1, sy2, ds; maxima=1, method=:linear, n=1) -> ::BeamformGrid

Compute normalised beam power at each horizontal slowness point
`sx1:ds:sx2` × `sy1:ds:sy2` s/°, between times `t1` and `t2` s.  The slowness
grid is spaced by `ds` s/°.

The `maxima` maximum beam power points are set in the returned `BeamformGrid`.

    beamform(s, t1, t2, smax, ds; maxima=1, method=:linear, n=1) -> ::BeamformGrid

Use a slowness grid of `-smax:ds:smax` in both dimensions.
"""
function beamform(s::TraceArray{T}, t1, t2, sx1, sx2, sy1, sy2, ds;
        maxima=1, method=:linear, n=nothing) where T
    sx = s_per_km.(sx1:ds:sx2)
    sy = s_per_km.(sy1:ds:sy2)
    nx, ny = length.((sx, sy))
    beam_power = Array{T}(undef, nx, ny)
    x, y, z = array_geometry(s.sta.lon, s.sta.lat, zeros(length(s)))
    for (iy, isy) in enumerate(sy), (ix, isx) in enumerate(sx)
        align = -(x.*isx + y.*isy)
        S = stack(s, (t1,t2), align; method=method, n=n)
        beam_power[ix,iy] = sum(Seis.trace(S).^2)
    end
    sx, sy = s_per_degree.(sx), s_per_degree.(sy)
    lon, lat = array_centre(s)
    sx_max, sy_max, p, β = find_maxima(sx, sy, beam_power, maxima)
    BeamformGrid{T}(lon, lat, x, y,
        collect(sx1:ds:sx2), collect(sy1:ds:sy2),
        beam_power./maximum(beam_power), maxima,
        sx_max, sy_max, p, β)
end
beamform(s::TraceArray, t1, t2, smax, ds; kwargs...) =
    beamform(s, t1, t2, -smax, smax, -smax, smax, ds; kwargs...)

"""
    array_geometry(lon, lat, alt) -> easting, northing, elevation
    array_geometry(s::AbstractArray{Seis.Trace}) -> easting, northing, elevation

Return the local `easting` and `northing` in km of each station in the array `s` in
a transerse Mercator transform, relative to the centre of the array.

`lon` and `lat` should be in ° and `alt` in m.  Or provide an array of `Seis.Trace`s.
"""
function array_geometry(lon, lat, alt)
    length(lon) == length(lat) == length(alt) ||
        throw(ArgumentError("`lon` and `lat` must be the same length"))
    # Find spherical mean of the array, ignoring altitude
    x, y, z = geog2cart(lon, lat, true)
    X, Y, Z = Float64.(sum.((x, y, z)))
    mean_lon, mean_lat, r = cart2geog(X, Y, Z, true)
    # Create transform into a local TM projection
    origin = Geodesy.LLA(mean_lat, mean_lon, zero(eltype(lat)))
    transform = Geodesy.ENUfromLLA(origin, Geodesy.wgs84)
    # Do the transform.  x and y are in metres
    x = Array{Float64}(undef, length(lon))
    y = similar(x)
    z = similar(x)
    for (i, (lo, la, al)) in enumerate(zip(lon, lat, alt))
        p = transform(Geodesy.LLA(la, lo, al))
        x[i] = p.e
        y[i] = p.n
        z[i] = p.u
    end
    x./1e3, y./1e3, z./1e3
end
array_geometry(s::TraceArray{T}) where T = array_geometry(s.sta.lon, s.sta.lat,
    [ismissing(ss.sta.elev) ? zero(T) : ss.sta.elev for ss in s])

"""
    array_centre(lon, lat) -> midlon, midlat

Return the mean location of the array in spherical coordinates.
`lon` and `lat` are vectors in °.
"""
function array_centre(lon, lat)
    length(lon) == length(lat) ||
        throw(ArgumentError("`lon` and `lat` must be the same length"))
    x, y, z = geog2cart(lon, lat, true)
    X, Y, Z = sum.((x, y, z))
    mlon, mlat, r = cart2geog(X, Y, Z, true)
    mlon, mlat
end
array_centre(s::TraceArray) = array_centre(s.sta.lon, s.sta.lat)

"""
    array_response(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df, s_deg=false) -> arf::ArrayResponse
    array_response(x, y, smax, ds, f1, f2, df, s_deg=false) -> arf::ArrayResponse
    array_response(traces::Array{Seis.Trace}, args...) -> arf::ArrayResponse

Return the array response function `arf.power` (an [`ArrayResponse`]@ref) at
slownesses with east and north components `arf.sx` and `arf.sy` respectively,
in s/km.  The ARF is computed between
frequencies `f1` and `f2` Hz, spaced by `df` Hz.

Provide either explicit limits on slowness `(sx1, sx2)` and `(sy1, sy2)` in the
east and north components respectively, or the maximum slowness `smax`,
plus the slowness spacing `ds`.

`x` and `y` are the station coordinates in km; or an array of `Seis.Trace` `traces`
can be given.

Default input units of slowness are s/km.  Give `s_deg` as `true`
to use s/°.
"""
function array_response(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df, s_deg=false)
    slowx = sx1:ds:sx2
    slowy = sy1:ds:sy2
    s_deg && ((slowx, slowy) = s_per_km.((slowx, slowy)))
    f = f1:df:f2
    arf = zeros(length(slowx), length(slowy))
    _compute_arf!(arf, x, y, slowx, slowy, f)
    ((slowx, slowy) = s_per_degree.((slowx, slowy)))
    ArrayResponse(f, x, y, slowx, slowy, arf)
end
array_response(x, y, smax, ds, f1, f2, df, s_deg::Bool=false) =
    array_response(x, y, -smax, smax, -smax, smax, ds, f1, f2, df, s_deg)
array_response(s::TraceArray, args...) =
    array_response(array_geometry(s)[1:2]..., args...)

"""    _compute_arf!(arf, x, y, sx, sy, f)

Compute the array response function `arf` at indices `[ix,iy]`, where `ix` and
`iy` refer to the slownesses `sx[ix]` and `sy[iy]` (s/km).  Integrate over the
frequencies (Hz) in `f`.  `x` and `y` are the station locations (km)."""
function _compute_arf!(arf::AbstractArray{T,2} where T, x, y, sx, sy, f)
    length(x) == length(y) && size(arf) == (length(sx), length(sy)) ||
        throw(ArgumentError("Array sizes do not match"))
    @inbounds for (j, slowy) in enumerate(sy), (i, slowx) in enumerate(sx)
        for k in 1:length(f)
            sum = 0.0im
            @simd for istat in 1:length(x)
                sum += exp((slowx.*x[istat] + slowy.*y[istat])*2π*f[k]*im)
            end
            arf[i,j] += abs(sum)^2
        end
    end
    nothing
end

"""
    delay_times(x, y, sx, sy) -> delays
    delay_times(s::AbstractArray{Seis.Trace}, sx, sy) -> delays

Return the delay times `delays` relative to the array centre for an incoming
plane wave with horizontal slownesses `sx` and `sy` s/km.  Stations have
coordinates `x` and `y` in the east and north directions respectively (km).
Or give an array of `Seis.Trace`s.
"""
function delay_times(x, y, sx, sy)
    length(x) == length(y) ||
        throw(ArgumentError("Length of `x` and `y` must be the same"))
    delays = Array{float(eltype(x))}(undef, length(x), length(sx), length(sy))
    @inbounds for (j, ssyy) in enumerate(sy), (i, ssxx) in enumerate(sx)
        @simd for k in eachindex(x)
            delays[k,i,j] = ssxx*x[k] + ssyy*y[k]
        end
    end
    delays
end
delay_times(s::TraceArray, sx, sy) = delay_times(array_geometry(s)[1:2]..., sx, sy)

"""
    find_maxima(sx, sy, power) -> sx_max, sy_max, slowness, backazimuth

Find `n` maxima in the beam power grid `power`, and return vectors
of `sx_max`, `sy_max`, `slowness` and `backazimuth` for each.
"""
function find_maxima(sx, sy, power::AbstractArray{T,2}, n=1) where T
    inds = sortperm(view(power, 1:length(power)), rev=true)
    sx_max, sy_max, p, β = T[], T[], T[], T[]
    for i in 1:n
        ix, iy = Tuple(CartesianIndices(power)[inds[i]])
        push!(sx_max, sx[ix])
        push!(sy_max, sy[iy])
        push!(p, sqrt(sx[ix]^2 + sy[iy]^2))
        push!(β, rad2deg(atan(sx[ix], sy[iy])))
    end
    sx_max, sy_max, p, β
end

"""
    fk(t::AbstractArray{T,2}, fs, x, y, starttime, endtime, smax, ds, f1, f2, s_deg=false) -> sx, sy, fk_stack
    fk(s::AbstractArray{Seis.Trace}, starttime, endtime, smax, ds, f1, f2, s_deg=false) -> sx, sy, fk_stack

Compute the frequency-wavenumber spectrum for the seismic traces in `s`.  Specify the
start and end times of the traces in s relative to the beginning of the trace in the
first form, or the B marker in the second.  The search is performed
up to a maximum slowness in s/km in steps of `slowness_step` s/km.  Integrate
the power between frequencies `f1` and `f2` Hz.

Slownesses are in units of s/km unless `s_deg` is `true`, in which case s/° are used.
"""
function fk(t, fs, x, y, t1, t2, sx1, sx2, sy1, sy2, ds, f1, f2, s_deg=false)
    sx = sx1:ds:sx2
    sy = sy1:ds:sy2
    if s_deg
        sx = s_per_km.(collect(sx))
        sy = s_per_km.(collect(sy))
    end
    nx, ny = length(sx), length(sy)
    fks = zeros(eltype(t), nx, ny)
    _calc_fk!(fks, t, fs, x, y, sx, sy, ds, f1, f2)
    sx1:ds:sx2, sy1:ds:sy2, fks
end
fk(t, fs, x, y, t1, t2, smax, ds, f1, f2, s_deg=false) =
    fk(t, fs, x, y, t1, t2, -smax, smax, -smax, smax, ds, f1, f2, s_deg)

function fk(s::TraceArray, args...)
    all(isapprox.(s[1].b, s.b, atol=s[1].delta/4)) &&
        all(Seis.nsamples(s[1]) .== Seis.nsamples.(s[:npts])) ||
        throw(ArgumentError("All traces must have same start time and length"))
    all(s[1].delta .== s.delta) ||
        throw(ArgumentError("All traces must have same sampling rate"))
    x, y, z = array_geometry(s)
    fk(hcat(s[:t]...), 1/s[1].delta, x, y, args...)
end

function _calc_fk!(fks::AbstractArray{T,2}, t::AbstractArray{Tt,2},
        fs, x, y, sx, sy, ds, f1, f2) where {T,Tt}
    size(fks) == (length(sx), length(sy)) ||
        throw(ArgumentError("Dimensions of `fks` do not match slownesses"))
    length(x) == length(y) || throw(ArgumentError("`x` and `y` must be same length"))
    # Create arrays for summation
    freqs, spec = _make_spectra(t, fs, f1, f2)
    shifted_spec = similar(spec)
    # Travel time delays
    delays = delay_times(x, y, sx, sy)
    nfreqs, ntraces = size(spec)
    local_delays = Array{T}(undef, ntraces)
    @inbounds for iy in eachindex(sy), ix in eachindex(sx)
        VERBOSE[] && @show ix, iy
        for i in 1:ntraces
            local_delays .= delays[:,ix,iy]
            _tshift!(shifted_spec, spec, local_delays, fs)
        end
        _abs!(shifted_spec)
        # fks[ix,iy] .= zero(T)
        # for j in 1:size(shifted_spec, 2), i in 1:size(shifted_spec, 1)
        #     fks[ix,iy] += shifted_spec[i,j]^2
        # end
        fks[ix,iy] = sum(abs(ss^2) for ss in shifted_spec)
    end
    fks ./= maximum(fks)
    fks
end

"""    _abs!(a::AbstractArray{T,2})

Replace each element of a with its absolute value."""
@inline function _abs!(a::AbstractArray)
    @inbounds @simd for i in 1:length(a)
        a[i] = abs(a[i])
    end
    a
end

"""   _make_spectra(t, fs, f1, f2) -> freqs, spec

Return an array `spec` to hold the Fourier-transformed traces between frequencies
`f1` and `f2` Hz, assuming the signal is real.  Traces are stored as columns of `t`
with sampling frequency `fs` Hz."""
function _make_spectra(t, fs, f1, f2)
    npts = size(t, 1)
    N = round(Int, npts/2) + 1
    df = fs/npts
    freqs = df*(1:N)
    inds = findall(f1 .<= freqs .<= f2)
    spec = fft(t, 1)
    VERBOSE[] && @show N, extrema(freqs), size(inds), size(spec)
    freqs[inds], spec[inds,:]
end

"""    _covariance(spec) -> c

Return the covariance matrix `c` for the spectra in `spec`.  These contain
a spectrum in each column.  `c` has dimensions `(n_freq, n_traces, n_traces)`,
thus to access the covariance between stations i and j, at frequency k, do
`c[k,i,j]`."""
function _covariance(spec)
    n_freq, n_traces = size(spec)
    c = Array{eltype(spec)}(undef, n_freq, n_traces, n_traces)
    @inbounds for j in 1:n_traces, i in j:n_traces
        c[:,i,j] .= spec[:,i].*conj.(spec[:,j])
        if i != j
            c[:,j,i] .= conj.(c[:,i,j])
        end
    end
    c
end

function _tshift!(shift_spec, spec, delay, fs)
    nf, nsta = size(spec)
    k = convert(eltype(shift_spec), -2π*im*fs/2/(nf-1))
    # shift_spec .= spec.*exp.(delay'.*(0:(nf-1)).*k)
    #
    pow = Float32.(-2π*fs*collect(0:(nf-1))/2/(nf-1))
    @inbounds for j in 1:nsta
        @simd for i in 1:nf
            shift_spec[i,j] = spec[i,j]*cis(delay[j]*pow[i])
        end
    end
    nothing
end

"""    s_per_degree(s_per_km, r_earth=$R_EARTH_KM_DEFAULT)

Convert slowness in s/km at the surface of the Earth to s/°."""
s_per_degree(s_km, r_earth=R_EARTH_KM_DEFAULT) = s_km*2π*r_earth/360

"""
    s_per_km(s_per_degree, r_earth=$R_EARTH_KM_DEFAULT)

Convert slowness in s/° at the surface of the Earth to s/km."""
s_per_km(s_degree, r_earth=R_EARTH_KM_DEFAULT) = s_degree*360/(2π*r_earth)

"List of available stacking routines"
const available_stack_methods = (:linear, :nthroot, :phaseweight)

"""
    stack(S::AbstractArray{Seis.Trace}, time_range, align=zeros(length(S)), weight; method=:linear) -> s::Seis.Trace

Return the linear stack of all traces in `S`, aligned in time on the
value in `align`, between `times[1]` and `times[end]`.

`times` can be a Tuple, Range or Array; in all cases, only the first and last values
are used.  **N.B. If using a non-integer Range, this may not include the end value
you expect.  E.g., `0.1:0.9` by default contains only the value `0.1` because the
default step is `1`.**

`align` can be a header symbol, or an array of equal length to `S`.  If it is a Symbol,
then the values in the headers are used.

All traces must have the same sampling frequency, and an error is thrown
if not all traces can be included in the stack.

If provided, traces are weighted by the values in `weight`, which are normalised
to sum to unity.

Available stacking methods are: $(available_stack_methods).
"""
function stack(S::AbstractArray{<:Seis.Trace{T}}, time_range,
               align::Union{Symbol,AbstractArray}=zeros(T, length(S));
               weights=ones(T, length(S))./length(S),
               method=:linear,
               n=nothing) where T
    # Convert to a vector of times if given as a symbol
    if typeof(align) == Symbol
        picks = getproperty(S.picks, align)
        any(ismissing, picks) && error("not all traces have a pick ")
        align = picks.time
    end
    all(x -> x == S[1].delta, S.delta) || error("all traces must have same delta")
    length(S) == length(weights) ||
        throw(ArgumentError("weights vector must be same length as traces"))
    delay = -align
    t1, t2 = time_range[1], time_range[end]
    weights ./= sum(weights) # Normalise weights
    if method == :linear
        stack_trace = _stack_linear([Seis.trace(s) for s in S],
            S[1].delta, S.b, delay, t1, t2, weights)
    elseif method == :nthroot
        n === nothing &&
            throw(ArgumentError("n must be provided for nthroot stacking"))
        stack_trace = _stack_nthroot([Seis.trace(s) for s in S],
            n, S[1].delta, S.b, delay, t1, t2, weights)
    elseif method == :phaseweight
        n === nothing &&
            throw(ArgumentError("n must be provided for phaseweighted stacking"))
        stack_trace = _stack_phaseweight([Seis.trace(s) for s in S],
            n, S[1].delta, S.b, delay, t1, t2, weights)
    else
        throw(ArgumentError("unrecognised stacking method"))
    end
    Seis.Trace{T}(t1, S[1].delta, stack_trace)
end

"""
    _stack_linear(a, delta, b, delay, t1, t2) -> stacked_trace::Vector{eltype(a)}

Return an array containing the linear stack of the vector of vectors `a`, whose
start times are given in `b`, and where the stack is returned between times `t1` and `t2`
(which can be negative if delay is also negative).  `delta` is the sampling interval
for all the traces.

The trace starts at time `t1` and has sample spacing `delta`.
"""
function _stack_linear(a::Vector{<:AbstractArray}, delta::Real, b::Vector,
        delay::Vector, t1::Real, t2::Real, weights)
    N, npts, b_shift, stack_trace = _stack_array(a, delta, b, delay, t1, t2)
    # stack_trace checks that we are always in array bounds
    @inbounds for j in 1:N
        ip1 = round(Int, (t1 - b_shift[j])/delta) + 1
        w = weights[j]
        @simd for i in ip1:(ip1+npts-1)
            stack_trace[i-ip1+1] += a[j][i]*w
        end
    end
    stack_trace
end

"""
    _stack_nthroot(a::Vector{<:AbstractArray}, n, delta, b, delay, t1, t2, weights)

Return an array containing the `n`throot stack of the vector of vectors `a`, whose
start times are given in `b`, and where the stack is returned between times `t1` and `t2`
(which can be negative if delay is also negative).  `delta` is the sampling interval
for all the traces.

The trace starts at time `t1` and has sample spacing `delta`.
"""
function _stack_nthroot(a::Vector{<:AbstractArray}, n, delta, b, delay, t1, t2, weights)
    N, npts, b_shift, stack_trace = _stack_array(a, delta, b, delay, t1, t2)
    @inbounds for j in 1:N
        ip1 = round(Int, (t1 - b_shift[j])/delta) + 1
        w = weights[j]
        @simd for i in ip1:(ip1+npts-1)
            stack_trace[i-ip1+1] += sign(a[j][i])*abs(a[j][i])^(1/n)*w
        end
    end
    stack_trace .= sign.(stack_trace).*abs.(stack_trace).^n
    stack_trace
end

"""
    _stack_phaseweight(a::Vector{<:AbstractArray}, delta, b, delay, t1, t2, weights)

Return an array containing the phaseweighted stack of the vector of vectors `a`, whose
start times are given in `b`, and where the stack is returned between times `t1` and `t2`
(which can be negative if delay is also negative).  `delta` is the sampling interval
for all the traces.  The weighting power is given by `n`.

The trace starts at time `t1` and has sample spacing `delta`.
"""
function _stack_phaseweight(a::Vector{<:AbstractArray{T}}, n, delta, b, delay, t1, t2, weights) where T
    N, npts, b_shift, stack_trace = _stack_array(a, delta, b, delay, t1, t2)
    # Collect shifted traces together
    traces = Array{T}(undef, npts, N)
    @inbounds for j in 1:N
        ip1 = round(Int, (t1 - b_shift[j])/delta) + 1
        w = weights[j]
        @simd for i in ip1:(ip1+npts-1)
            traces[i-ip1+1,j] = a[j][i]*w
        end
    end
    # Create instantaneous phases
    phases = DSP.hilbert(traces)
    phases .= phases./abs.(phases)
    # Stack
    @inbounds for i in 1:npts
        weight = abs(sum(phases[i,:])/N)^n
        @simd for j in 1:N
            stack_trace[i] += traces[i,j]*weight/N
        end
    end
    stack_trace
end

"""
    _stack_array(a::Vector{Vector{<:Real}}, delta, b::Vector, delay::Vector, t1, t2) -> N, npts, b_shifted, stack_trace

Return the number of traces `N`, the number of points in the output stack trace,
`npts`, an array of shifted trace start times `b_shifted` and a zeroed array
for the stack itself, `stack_trace`.  This routine also checks that the input
arrays are the required length and no points are missing from the desired stack.
"""
function _stack_array(a::Vector{A}, delta, b, delay, t1, t2) where {A<:AbstractArray}
    N = length(a)
    length(b) == length(delay) == N ||
        throw(ArgumentError("Lengths of `a`, `b` and `delay` must all be the same"))
    t1 < t2 || throw(ArgumentError("`t2` is before or the same as `t1`"))
    b_shift = b .+ delay
    e_shift = b_shift .+ (length.(a) .- 1).*delta
    all(b_shift .<= t1) && all(t2 .<= e_shift) ||
        error("Not all traces have data within time range $(t1):$(t2) s
        range in b: $(extrema(b_shift))
        range in e: $(extrema(e_shift))")
    npts = round(Int, (t2 - t1)/delta) + 1
    stack_trace = zeros(eltype(eltype(a)), npts)
    N, npts, b_shift, stack_trace
end
