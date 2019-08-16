"""
    BeamformGrid{T}

Struct containing the results of beamforming on a regular slowness grid and
array coordinates.
"""
struct BeamformGrid{T}
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
end
