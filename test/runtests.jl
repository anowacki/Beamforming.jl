using Beamforming
using Test

@testset "All tests" begin
    include("fk.jl")
    include("stack.jl")
    include("vespa.jl")
    include("synth.jl")
    include("arf.jl")
end
