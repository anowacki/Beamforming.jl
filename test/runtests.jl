using Beamforming
using Test

@testset "All tests" begin
    include("stack.jl")
    include("beamforming.jl")
    include("vespa.jl")
    include("synth.jl")
    include("arf.jl")
end
