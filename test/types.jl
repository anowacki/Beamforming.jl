using Test

using Beamforming
using Seis

@testset "Types" begin
    n = 5
    sx = -2:0.1:3
    sy = -1:0.1:2
    x = 1:n
    y = 0:(n - 1)
    lon = rand()
    lat = rand()
    freqs = 1:0.2:2
    @testset "$T" for T in (Float32, Float64)
        power = zeros(T, length(sx), length(sy))
        @testset "==, isequal, hash" begin
            @testset "BeamformGrid" begin
                fk = Beamforming.BeamformGrid{T, Trace{T, Vector{T}, Seis.Geographic{T}}}(
                    lon, lat, x, y, sx, sy, power, 1, [1], [2], [1], [120], [],
                    :linear, nothing, -10, 20, :normal, [], nothing
                    )
                fk′ = deepcopy(fk)
                @test fk == fk′
                @test fk !== fk′
                @test isequal(fk, fk′)
                @test hash(fk) == hash(fk′)
                @test unique([fk, fk′]) == [fk]
                # Modify
                fk′.power[1] = 1
                @test fk != fk′
                @test !isequal(fk, fk′)
                @test unique([fk, fk′]) == [fk, fk′]
            end

            @testset "$ARF" for ARF in (Beamforming.ArrayResponse, Beamforming.ArrayResponseCrosscorrelation)
                arf = ARF{T}(freqs, x, y, sx, sy, power)
                arf′ = deepcopy(arf)
                @test arf == arf′
                @test arf !== arf′
                @test isequal(arf, arf′)
                @test hash(arf) == hash(arf′)
                @test unique([arf, arf′]) == [arf]
                # Modify
                arf′.power[1] = 1
                @test arf != arf′
                @test !isequal(arf, arf′)
                @test unique([arf, arf′]) == [arf, arf′]
            end

            @testset "VespaGrid" begin
                times = freqs
                slows = sx
                grid = zeros(T, length(times), length(slows))
                vg = Beamforming.VespaGrid{T, Trace{T, Vector{T}, Seis.Geographic{T}}}(
                    lon, lat, x, y, times, slows,
                    grid, 1, [1], [1], [], :nthroot, 3, :plane, false)
                vg′ = deepcopy(vg)
                @test vg == vg′
                @test vg !== vg′
                @test isequal(vg, vg′)
                @test hash(vg) == hash(vg′)
                @test unique([vg, vg′]) == [vg]
                # Modify
                vg.power[1] = 1
                @test vg != vg′
                @test !isequal(vg, vg′)
                @test unique([vg, vg′]) == [vg, vg′]
            end
        end
    end
end