"""
# Beamforming

Compute beam power for an array of stations using the `beamform` function,
and compute array response functions using `array_response_function`.
"""
module Beamforming

import DSP, Geodesy

import Seis

export
    array_response_function,
    beamform,
    stack

# Custom type aliases
const TraceArray{T} = AbstractArray{<:Seis.Trace{T}}

# Constants
const R_EARTH_KM_DEFAULT = 6371.0

"Global verbosity flag.  Set with Beamforming.verbose(::Bool)"
const VERBOSE = Ref(true)

"Set module verbosity"
verbose(tf::Bool) = VERBOSE[] = tf

include("types.jl")
include("sphericalgeom.jl")
include("stacking.jl")
include("plots.jl")

# Submodules
include("Synth.jl")
import .Synth

end # module
