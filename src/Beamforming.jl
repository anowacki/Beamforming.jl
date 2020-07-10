"""
# Beamforming

Compute beam power for an array of stations using the `beamform` function,
and compute array response functions using `array_response_function`.
"""
module Beamforming

using Statistics: mean

import DSP, Geodesy

import Seis, Geodesics

export
    array_response,
    beamform,
    crosscorrelation_array_response,
    crosscorrelation_beamform,
    crosscorrelation_beamform_corrs,
    stack,
    vespagram

# Custom type aliases
const TraceArray{T} = AbstractArray{<:Seis.Trace{T}}

# Constants
const R_EARTH_KM_DEFAULT = 6371.0

"Global verbosity flag.  Set with Beamforming.verbose(::Bool)"
const VERBOSE = Ref(true)

"Set module verbosity"
verbose(tf::Bool) = VERBOSE[] = tf

include("types.jl")
include("util.jl")
include("sphericalgeom.jl")
include("stacking.jl")
include("beamform.jl")
include("correlation.jl")
include("array_response.jl")
include("plots.jl")

# Submodules
include("Synth.jl")
import .Synth

end # module
