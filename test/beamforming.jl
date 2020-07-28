using Test

using Seis
using Beamforming, Beamforming.Synth

@testset "Beamforming" begin
    @testset "Conventional" begin
        let sxi = 2rand()-1, syi = 2rand()-1
            s = Synth.synthetic_arrival(sxi, syi)
            nsta = length(s)
            s_max = 2
            ds = 0.1
            @testset "$method" for (method, n) in zip(
                    (:linear, :nthroot, :phaseweight), (nothing, 3, 3))
                bf = Beamforming.beamform(s, 40, 60, -s_max, s_max, -s_max, s_max, ds)
                @test bf isa Beamforming.BeamformGrid
                @test length(bf.x) == length(bf.y) == nsta
                @test bf.sx == -s_max:ds:s_max
                @test bf.sy == -s_max:ds:s_max
                @test size(bf.power) == (length(bf.sx), length(bf.sy))
                @test bf.nmax == 1
                @test length(bf.bazmax) == length(bf.sx_max) ==
                    length(bf.sy_max) == length(bf.pmax)
                @test bf.sx_max[1] ≈ sxi atol=ds
                @test bf.sy_max[1] ≈ syi atol=ds
            end
        end
    end

    @testset "Cross-correlation BF" begin
        let sx = 20*(rand() - 0.5), sy = 20*(rand() - 0.5), nsta = 8
            s = Synth.synthetic_arrival(sx, sy, nsta=nsta,
                lon=2rand(nsta), lat=2rand(nsta), distance=90)
            t1, t2 = maximum(starttime.(s)), minimum(endtime.(s))
            s_max = 11
            ds = 0.5
            t = 25
            bf = crosscorrelation_beamform(s, t1, t2,
                -s_max, s_max, -s_max, s_max, ds, t=t)
            @test bf isa Beamforming.BeamformGrid
            @test length(bf.x) == length(bf.y) == nsta
            @test bf.sx == -s_max:ds:s_max
            @test bf.sy == -s_max:ds:s_max
            @test size(bf.power) == (length(bf.sx), length(bf.sy))
            @test bf.nmax == 1
            @test length(bf.bazmax) == length(bf.sx_max) ==
                length(bf.sy_max) == length(bf.pmax)
            @test bf.sx_max[1] ≈ sx atol=ds
            @test bf.sy_max[1] ≈ sy atol=ds
            # Check the max slowness version is the same
            bf′ = crosscorrelation_beamform(s, t1, t2, s_max, ds, t=t)
            @testset "$f" for f in fieldnames(Beamforming.BeamformGrid)
                @test getfield(bf, f) == getfield(bf′, f)
            end
        end
    end

    @testset "Cross-correlation BF with correlations" begin
        let sx = 5*(rand() - 0.5), sy = 5*(rand() - 0.5), nsta = 8
            s = Synth.synthetic_arrival(sx, sy, nsta=nsta,
                lon=rand(nsta), lat=rand(nsta), distance=30)
            t1, t2 = 40, 60
            corrs = Beamforming.compute_crosscorrelations(s, t1, t2)
            win = 1
            smax = 3
            ds = 0.1
            bf = crosscorrelation_beamform(s, t1, t2, smax, ds, t=win)
            bf′ = crosscorrelation_beamform_corrs(corrs, win, smax, ds)
            # Skip some fields in comparison
            @testset "$f" for f in fieldnames(Beamforming.BeamformGrid)
                f in (:traces, :t1, :t2) && continue
                if f in (:x, :y)
                    # Order of unique stations will be different,
                    # but values should be the same.
                    @test sort(getfield(bf, f)) ≈ sort(getfield(bf′, f))
                else
                    @test getfield(bf, f) == getfield(bf′, f) ||
                        # Catch floating-point roundoff errors in different order
                        getfield(bf, f) ≈ getfield(bf′, f)
                end
            end
        end
    end
end
