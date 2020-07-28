"""
    array_geometry(lon, lat, alt) -> easting, northing, elevation, mean_lon, mean_lat
    array_geometry(traces::AbstractArray{Seis.Trace}) -> easting, northing, elevation, mean_lon, mean_lat
    array_geometry(stations::AbstractArray{Seis.Station}) -> easting, northing, elevation, missing, missing

Return the local `easting` and `northing` in km of each station in the array `s` in
a transverse Mercator transform, relative to the centre of the array.  The longitude
and latitude of the array centre are given by `mean_lon` and `mean_lat` (°) respectively.

`lon` and `lat` should be in ° and `alt` in m.  Or provide an array of `Seis.Trace`s
or `Seis.Station`s.  If `traces` or `stations` is an array of `CartTrace` or
`CartStation`s, then no mean longitude or latitude are calculated and instead
they are returned as `nothing`.
"""
function array_geometry(lon, lat, alt)
    T = eltype(lon)
    length(lon) == length(lat) == length(alt) ||
        throw(DimensionMismatch("`lon` and `lat` must be the same length"))
    # Find spherical mean of the array, ignoring altitude
    x, y, z = geog2cart(lon, lat, true)
    X, Y, Z = Float64.(sum.((x, y, z)))
    mean_lon, mean_lat, r = cart2geog(X, Y, Z, true)
    # Create transform into a local TM projection
    origin = Geodesy.LLA(mean_lat, mean_lon, zero(eltype(lat)))
    transform = Geodesy.ENUfromLLA(origin, Geodesy.wgs84)
    # Do the transform.  x and y are in metres
    x = Array{T}(undef, length(lon))
    y = similar(x)
    z = similar(x)
    for (i, (lo, la, al)) in enumerate(zip(lon, lat, alt))
        p = transform(Geodesy.LLA(la, lo, al))
        x[i] = p.e/1e3
        y[i] = p.n/1e3
        z[i] = p.u/1e3
    end
    x, y, z, mean_lon, mean_lat
end

array_geometry(s::TraceArray) = array_geometry(s.sta)

array_geometry(stas::AbstractArray{<:Seis.GeogStation{T}}) where T =
    array_geometry(stas.lon, stas.lat, coalesce.(stas.elev, zero(T)))

function array_geometry(stas::AbstractArray{<:Seis.CartStation{T}}) where T
    to_km = T(0.001)
    xx = stas.x
    yy = stas.y
    zz = stas.z
    any(x -> x[1] === missing || x[2] === missing, zip(xx, yy)) &&
        throw(ArgumentError("all stations must have both x and y coordinates"))
    x = to_km*(xx .- mean(xx))
    y = to_km*(yy .- mean(yy))
    z = to_km*coalesce.((zz .- mean(zz)), zero(T))
    x, y, z, nothing, nothing
end

"""
    array_centre(lon, lat) -> midlon, midlat

Return the mean location of the array in spherical coordinates.
`lon` and `lat` are vectors in °.
"""
function array_centre(lon, lat)
    length(lon) == length(lat) ||
        throw(DimensionMismatch("`lon` and `lat` must be the same length"))
    x, y, z = geog2cart(lon, lat, true)
    X, Y, Z = sum.((x, y, z))
    mlon, mlat, r = cart2geog(X, Y, Z, true)
    mlon, mlat
end
array_centre(s::TraceArray) = array_centre(s.sta.lon, s.sta.lat)

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
    inds = sortperm(vec(power), rev=true)
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

"List of available stacking routines"
const AVAILABLE_STACK_METHODS = (:linear, :nthroot, :phaseweight)

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

Available stacking methods are: $(AVAILABLE_STACK_METHODS).
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
