using Test

using Seis
using Beamforming, Beamforming.Synth


@testset "Form beam" begin
    let sxi = 2rand()-1, syi = 2rand()-1
        s = Synth.synthetic_arrival(sxi, syi)
        nsta = length(s)
        bf = Beamforming.beamform(s, 40, 60, -2, 2, -2, 2, 0.1)
        @test bf isa Beamforming.BeamformGrid
        @test length(bf.x) == length(bf.y) == nsta
        @test bf.sx == -2:0.1:2
        @test bf.sy == -2:0.1:2
        @test size(bf.power) == (length(bf.sx), length(bf.sy))
        @test bf.nmax == 1
        @test length(bf.bazmax) == length(bf.sx_max) ==
            length(bf.sy_max) == length(bf.pmax)
        @test bf.sx_max[1] ≈ sxi atol=0.05
        @test bf.sy_max[1] ≈ syi atol=0.05
    end
end
