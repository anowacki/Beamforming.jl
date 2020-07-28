using Test

using Beamforming
using Seis

# Test array responses are calculated for a zero-slowness set of plane waves
# at 3, 4, 5, 6 and 7 Hz for three stations of the ARCES array in Norway
# and can be compared to Ruigrok et al., J Seismol, 2017.
# Data are written little-endian
@testset "Array response" begin
    @testset "$T" for T in (Float32, Float64)
        # Geographic stations
        ara1 = Station{T}(net="NO", sta="ARA1", lon=25.5071, lat=69.5363)
        ara2 = Station{T}(net="NO", sta="ARA2", lon=25.5078, lat=69.5338)
        arb2 = Station{T}(net="NO", sta="ARB2", lon=25.5134, lat=69.5357)

        stas = [ara1, ara2, arb2]

        # Geographic traces
        n = length(stas)
        t = [Trace{T}(0, 1, 0) for _ in 1:n]
        t.sta = stas

        s_max = 1
        ds = 0.01
        f1 = 3
        f2 = 7
        df = 1

        x, y, z, lon, lat = Beamforming.array_geometry(t)

        # Cartesian stations
        cstas = [Station{T,Seis.Cartesian{T}}() for _ in 1:n]
        cstas.net = stas.net
        cstas.sta = stas.sta
        cstas.x = x.*1000
        cstas.y = y.*1000
        cstas.z = z.*1000

        # Cartesian traces
        ct = [CartTrace{T}(0, 1, 0) for _ in 1:n]
        ct.sta = cstas

        @testset "Conventional ARF" begin
            true_power = ltoh.(reinterpret(Float64,
                read(joinpath(@__DIR__, "data", "arf_conventional_$(f1)_$(f2)_$(df).bin"))))
            true_power = reshape(true_power, Int(sqrt(length(true_power))), :)
            arf = array_response(t, s_max, ds, f1, f2, df)
            @test arf.power ≈ T.(true_power) rtol=5e-3
            @test arf isa Beamforming.ArrayResponse{T}

            @test array_response(stas, s_max, ds, f1, f2, df) == arf
            
            # Cartesian coordinates
            carf = array_response(ct, s_max, ds, f1, f2, df)
            @testset "$f" for f in fieldnames(typeof(carf))
                @test getfield(carf, f) == getfield(arf, f) ||
                    ≈(getfield(carf, f), getfield(arf, f), rtol=1e-3)
            end
            csarf = array_response(cstas, s_max, ds, f1, f2, df)
            @test carf == csarf

            @test arf.sx == -s_max:ds:s_max
            @test arf.sy == -s_max:ds:s_max
            @test arf.freqs == f1:df:f2
            @test arf.x == x
            @test arf.y == y
        end

        @testset "Cross-correlation ARF" begin
            true_power = ltoh.(reinterpret(Float64,
                read(joinpath(@__DIR__, "data", "arf_crosscorrelation_$(f1)_$(f2)_$(df).bin"))))
            true_power = reshape(true_power, Int(sqrt(length(true_power))), :)
            arf = crosscorrelation_array_response(t, s_max, ds, f1, f2, df)
            @test arf.power ≈ T.(true_power) rtol=5e-3
            @test arf isa Beamforming.ArrayResponseCrosscorrelation{T}

            @test crosscorrelation_array_response(stas, s_max, ds, f1, f2, df) == arf
            @test crosscorrelation_array_response(ct, s_max, ds, f1, f2, df) == arf
            @test crosscorrelation_array_response(cstas, s_max, ds, f1, f2, df) == arf

            α, Δ = Beamforming._azimuth_distance_ccbf(x, y)
            @test all(arf.sx .== T.(-s_max:ds:s_max))
            @test all(arf.sy .== T.(-s_max:ds:s_max))
            @test all(arf.freqs .== f1:df:f2)
            @test arf.azimuth == α
            @test arf.distance == Δ
        end
    end

    @testset "_azimuth_distance_ccbf" begin
        x = [0, 1]
        y = [1, 2]
        α = [45, 45-180]
        Δ = [√2, √2]
        @test Beamforming._azimuth_distance_ccbf(x, y) == (α, Δ)
    end
end
