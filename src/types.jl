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
    "Beam power of nmax maxima"
    pmax::Vector{T}
    "Backazimuth of nmax maxima (°)"
    bazmax::Vector{T}
    "Input traces.  May be empty"
    traces::Vector{S}
    "Stackin method (e.g., :linear, :phaseweight)"
    method::Symbol
    "Order of stack (e.g. in nthroot or phaseweighted stacking)"
    n::Union{Int,Nothing}
    "Start time of stack"
    t1::T
    "End time of stack"
    t2::T
end

"Return `true` if a `BeamformGrid` contains the traces used to construct it"
_hastraces(bf::BeamformGrid) = length(bf.traces) > 0

"""
    ArrayResponse{T}

Struct containing the array response function for a set of stations.
Slownesses in an `ArrayResponse` are always in s/°, and coordinates are
in km.

This can be plotted using [`Plots.plot`](@ref).
"""
struct ArrayResponse{T}
    "Frequencies at which the array responce was evaluated (Hz"
    freqs
    "Cartesian x-coordinates in each station (km)"
    x
    "Cartesian y-coordinates in each station (km)"
    y
    "Horizontal slownesses in s/° in x-direction of beamforming grid (first dimension)"
    sx
    "Horizontal slownesses in s/° in y-direction of beamforming grid (second dimension)"
    sy
    "Power of array response at each (sx, sy) point"
    power::Array{T,2}
end
