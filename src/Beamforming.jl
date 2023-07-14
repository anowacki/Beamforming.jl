"""
# Beamforming

Stack seismic data in various ways, based on [Seis.jl](https://githib.com/anowacki/Seis.jl).

## Exported functions
- `array_response`: Compute the array response function for an array.
- `beamform`: Find the beam power on a grid of slownesses ('beampacking').
- `stack`: Stack up traces using various methods.
- `vespagram`: Form a vespagram (slowness–time).
- The above, but using cross-correlation beamforming:
  - `crosscorrelation_array_response`: Array response.
  - `crosscorrelation_beamform`: Compute cross-correlations between traces, then
    compute beam power on a slowness grid.
  - `crosscorrelation_beamform_corrs`: Use an existing set of cross-correlations
    to compute beam power on a slowness grid.
"""
module Beamforming

import DSP
import Geodesics
import Geodesy
import Seis

using Statistics: mean


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
