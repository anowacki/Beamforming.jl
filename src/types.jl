"""
    BeamformGrid{T, S<:Seis.Trace}

Struct containing the results of beamforming on a regular slowness grid and
array coordinates.

This can be plotted using [`Plots.plot`](@ref).
"""
struct BeamformGrid{T, S<:Seis.Trace}
    "Reference longitude in ° of array"
    lon::T
    "Reference latitude in ° array"
    lat::T
    "Cartesian x-coordinates of each station (km)"
    x::Vector{T}
    "Cartesian y-coordinates of each station"
    y::Vector{T}
    "Horizontal slownesses in s/° in x-direction of beamforming grid (first dimension)"
    sx::Vector{T}
    "Horizontal slownesses in s/° in y-direction of beamforming grid (second dimension)"
    sy::Vector{T}
    "Beam power at each (sx, sy) point"
    power::Array{T,2}
    "Number of maxima identified"
    nmax::Int
    "X-slownesses of nmax maxima (s/°)"
    sx_max::Vector{T}
    "Y-slownesses of nmax maxima (s/°)"
    sy_max::Vector{T}
    "Absolute slowness of nmax maxima (s/°)"
    pmax::Vector{T}
    "Backazimuth of nmax maxima (°)"
    bazmax::Vector{T}
    "Input traces.  May be empty"
    traces::Vector{S}
    "Stacking method (e.g., :linear, :phaseweight)"
    method::Symbol
    "Order of stack (e.g. in nthroot or phaseweighted stacking)"
    n::Union{Int,Nothing}
    "Start time of stack (`nothing` if not available)"
    t1::Union{T,Nothing}
    "End time of stack (`nothing` if not available)"
    t2::Union{T,Nothing}
    # Cross-correlation beamforming
    "Kind of beamforming applied.  Can be one of:
     - `:normal`: ordinary beamforming on single traces.
     - `:crosscorrelation`: cross-correlation beamforming on receiver pairs"
    kind::Symbol
    "Cross-correlation traces.  May be empty"
    corrs::Vector{S}
    "Cross-correlation stacking window length (s)"
    corr_t::Union{T,Nothing}
end

# Keyword argument constructor
function BeamformGrid{T,S}(; lon, lat, x, y, sx, sy, power,
        nmax, sx_max, sy_max, pmax, bazmax, traces=S[], method, n=nothing,
        t1, t2, kind=:normal, corrs=S[], corr_t=nothing) where {T,S}
    BeamformGrid{T,S}(lon, lat, x, y, sx, sy, power, nmax, sx_max, sy_max, pmax,
        bazmax, traces, method, n, t1, t2, kind, corrs, corr_t)
end

# Default-argument contructor
BeamformGrid{T,S}(lon, lat, x, y, sx, sy, power, nmax, sx_max,
    sy_max, pmax, bazmax, traces, method, n, t1, t2, kind=:normal, corrs=S[],
    corr_t=nothing) where {T,S} =
    BeamformGrid{T,S}(lon, lat, x, y, sx, sy, power, nmax, sx_max,
        sy_max, pmax, bazmax, traces, method, n, t1, t2, king, corrs, corr_t)

"Return `true` if a `BeamformGrid` contains the traces used to construct it"
_hastraces(bf::BeamformGrid) = !isempty(bf.traces)

"Return `true` if a `BeamformGrid` contains the cross-correlation traces calculated
in its construction"
_hascorrs(bf::BeamformGrid) = !isempty(bf.corrs)

"""
    AbstractArrayResponse

Supertype of all array response types.  All subtypes should contain the following
fields:

- `freqs`: Frequencies at which the array response was evaluated (Hz)
- `sx`: Horizontal slowness in s/° in x-direction of grid (first dimension)
- `sy`: Horizontal slowness in s/° in y-direction of grid (second dimension)
- `power`: Power of array response at each `(sx, sy)` point.
"""
abstract type AbstractArrayResponse end

"""
    ArrayResponse{T}

Struct containing the array response function for a set of stations.
Slownesses in an `ArrayResponse` are always in s/km, and coordinates are
in km.

This can be plotted using [`Plots.plot`](@ref).
"""
struct ArrayResponse{T} <: AbstractArrayResponse
    "Frequencies at which the array response was evaluated (Hz)"
    freqs::Vector{T}
    "Cartesian x-coordinates in each station (km)"
    x::Vector{T}
    "Cartesian y-coordinates in each station (km)"
    y::Vector{T}
    "Horizontal slownesses in s/km in x-direction of beamforming grid (first dimension)"
    sx::Vector{T}
    "Horizontal slownesses in s/km in y-direction of beamforming grid (second dimension)"
    sy::Vector{T}
    "Power of array response at each (sx, sy) point"
    power::Array{T,2}
end

"""
    ArrayResponse{T}(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df)

Construct an empty `ArrayResponse` object where the power is set to be zero
at all points.  This is useful when later passing the object to a function
to compute the array response function.
"""
function ArrayResponse{T}(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df) where T
    length(x) == length(y) ||
        throw(DimensionMismatch("x and y coordinate vectors not the same length"))
    freqs = f1:df:f2
    sx = sx1:ds:sx2
    sy = sy1:ds:sy2
    power = zeros(T, length(sx), length(sy))
    ArrayResponse{T}(freqs, x, y, sx, sy, power)
end

"""
    ArrayResponseCrosscorrelation{T}

Struct containing the cross-correlation array response function for a set of
interstation paths.  Slownesses in an `ArrayResponseCrosscorrelation` are
always in s/km, azimuths are in °, and distances are in km.
"""
struct ArrayResponseCrosscorrelation{T} <: AbstractArrayResponse
    "Frequencies at which the array response was evaluated (Hz)"
    freqs::Vector{T}
    "Azimuth between each interstation pair (°)"
    azimuth::Vector{T}
    "Distance between each pair of stations (km)"
    distance::Vector{T}
    "Horizontal slownesses in s/km in x-direction of beamforming grid (first dimension)"
    sx::Vector{T}
    "Horizontal slownesses in s/km in y-direction of beamforming grid (second dimension)"
    sy::Vector{T}
    "Power of cross-correlation array response at each `(sx, sy)` point"
    power::Array{T,2}
end

"""
    ArrayResponseCrosscorrelation{T}(azimuth, distance, sx1, sx2, sy1, sy2, ds, f1, f2, df)

Construct an empty `ArrayResponseCrosscorrelation` object where the power
is set to be zero at all points.  This is useful when later passing the
object to a function to compute the array response function.
"""
function ArrayResponseCrosscorrelation{T}(azimuth, distance,
        sx1, sx2, sy1, sy2, ds, f1, f2, df) where T
    length(azimuth) == length(distance) ||
        throw(DimensionMismatch("azimuth and distance vectors not the same length"))
    sx = sx1:ds:sx2
    sy = sy1:ds:sy2
    freqs = f1:df:f2
    power = zeros(T, length(sx), length(sy))
    ArrayResponseCrosscorrelation{T}(freqs, azimuth, distance, sx, sy, power)
end

"""
    VespaGrid{T, S<:Seis.Trace}

Struct containing the results of vespagram analysis.

This can be plotted using [`Plots.plot`](@ref).
"""
struct VespaGrid{T, S<:Seis.Trace}
    "Reference longitude in ° of array"
    lon::T
    "Reference latitude in ° array"
    lat::T
    "Cartesian x-coordinates of each station (km)"
    x::Vector{T}
    "Cartesian y-coordinates of each station"
    y::Vector{T}
    "Stacking times (s)"
    time::Vector{T}
    "Slownesses (s/°)"
    slow::Vector{T}
    "Beam power at each (time, slow) point"
    power::Array{T,2}
    "Number of maxima identified"
    nmax::Int
    "Times of nmax maxima"
    t_max::Vector{T}
    "Slownesses of nmax maxima"
    s_max::Vector{T}
    "Input traces.  May be empty."
    traces::Vector{S}
    "Stacking method (e.g., :linear, :nthroot, :phaseweight)"
    method::Symbol
    "Order of stack (e.g. in nthroot or phaseweighted stacking)"
    n::Union{Nothing,Int}
    "Stacking wavefront type (e.g., :plane, :circle)"
    wavefront::Symbol
    "Whether vespagram is an envelope stack"
    envelope::Bool
end

"Return `true` if a `BeamformGrid` contains the traces used to construct it"
_hastraces(vesp::VespaGrid) = !isempty(vesp.traces)

# Define ==, isequal and hash for our types
for type in (BeamformGrid, ArrayResponse, ArrayResponseCrosscorrelation, VespaGrid)
    T = nameof(type)
    fields = fieldnames(type)
    @eval begin
        function Base.:(==)(x::$T, y::$T)
            # For parameterised types, since T is a UnionAll
            typeof(x) == typeof(y) || return false
            $([:(x.$f == y.$f || return false) for f in fields]...)
            true
        end

        function Base.isequal(x::$T, y::$T)
            # For parameterised types, since T is a UnionAll
            typeof(x) == typeof(y) || return false
            $([:(isequal(x.$f, y.$f) || return false) for f in fields]...)
            true
        end

        function Base.hash(x::$T, h::UInt)
            $([:(h = hash(x.$f, h)) for f in fields]...)
            h
        end
    end
end
