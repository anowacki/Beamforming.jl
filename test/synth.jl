using Test

using Beamforming
using Beamforming.Synth
import Seis

@testset "Synth" begin
    @testset "synthetic_arrival" begin
        let t = synthetic_arrival(1, 2)
            @test length(t) == 20
            @test (-2.5, -2.5) <= extrema(t.sta.lon) <= (2.5, 2.5)
            @test t[1].delta == 0.01
            @test minimum(minimum, trace.(t)) == 0
            @test all(x -> Seis.nsamples(x) == 10001, t)
            @test all(x -> x === missing, t.evt.lon)
        end
        let t = synthetic_arrival(1, 2, delta=0.02, nsta=10, lon=rand(10),
                lat=rand(10), distance=80, noise=0.1)
            @test length(t) == 10
            @test (0, 0) <= extrema(t.sta.lon) <= (1, 1)
            @test t[1].delta == 0.02
            @test minimum(minimum, trace.(t)) < 0
            @test all(x -> Seis.nsamples(x) == 10001, t)
            @test all(x -> x !== missing, t.evt.lon)
            @test Beamforming.delta(t[1].evt.lon, t[1].evt.lat,
                t[1].sta.lon, t[1].sta.lon, true) ≈ 80 atol=2
            @test all(isapprox.(Seis.backazimuth.(t), atand(1, 2), atol=2))
        end
    end
end