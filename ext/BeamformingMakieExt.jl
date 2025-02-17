module BeamformingMakieExt

using Beamforming: Beamforming, AbstractArrayResponse, BeamformGrid, VespaGrid
import Geodesics

@static if isdefined(Base, :get_extension)
    import Makie
else
    import ..Makie
end

include("generic.jl")
include("plot_beamforming.jl")
include("plot_arf.jl")
include("plot_vespagram.jl")

end # module
