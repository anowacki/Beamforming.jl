using Test

using Seis
using Beamforming, Beamforming.Synth

@testset "Vespagrams" begin
    let slow = 10rand(), baz = 360rand(), maxima = rand(1:5), Δ = 100,
            wavefront = :plane, envelope = false
        sx, sy = slow.*(sind(baz), cosd(baz))
        s = Synth.synthetic_arrival(sx, sy, distance=Δ, noise=0.001)
        nsta = length(s)
        for (method, n) in zip((:linear, :nthroot, :phaseweight), (nothing, 3, 3))
            vespa = vespagram(s, 40, 60, 0, 10, 0.1, maxima=maxima, method=method, n=n)
            @test vespa isa Beamforming.VespaGrid
            @test length(vespa.x) == length(vespa.y) == nsta
            @test vespa.slow ≈ collect(0:0.1:10)
            @test vespa.time ≈ collect(40:s[1].delta:60)
            @test size(vespa.power) == (length(vespa.time), length(vespa.slow))
            @test vespa.nmax == maxima
            @test length(vespa.t_max) == length(vespa.s_max) == vespa.nmax
            @test vespa.method == method
            @test vespa.n == n
            @test vespa.envelope == envelope
            @test vespa.wavefront == wavefront
            @test vespa.s_max[1] ≈ slow atol=0.2
            @test vespa.t_max[1] ≈ 50 atol=15s[1].delta
        end
    end
end
