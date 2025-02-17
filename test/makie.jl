using Test

using Beamforming
using Seis

@testset "Makie plots" begin
    t = filter(is_vertical, sample_data(:regional))

    bf = beamform(t, 100, 150, 1, 0.2)
    arf = array_response(t, 1, 0.2, 0.1, 1, 0.1, true)
    vesp = vespagram(t, 100, 150, 1, 0.1)

    @testset "Makie not loaded" begin
        @testset "New plot" begin
            @testset "plot_array_response[!]" begin
                @test_throws (
                    "Note: If you have not already, you must first load a Makie " *
                    "backend such as GLMakie, CairoMakie or WGLMakie before" *
                    "you can use plot_array_response! to create plots."
                ) plot_array_response(arf)
                @test_throws (
                    "Note: If you have not already, you must first load a Makie " *
                    "backend such as GLMakie, CairoMakie or WGLMakie before" *
                    "you can use plot_array_response! to create plots."
                ) plot_array_response!(arf)
            end

            @testset "plot_beamforming[!]" begin
                @test_throws (
                    "Note: If you have not already, you must first load a Makie " *
                    "backend such as GLMakie, CairoMakie or WGLMakie before" *
                    "you can use plot_beamforming! to create plots."
                ) plot_beamforming(arf)
                @test_throws (
                    "Note: If you have not already, you must first load a Makie " *
                    "backend such as GLMakie, CairoMakie or WGLMakie before" *
                    "you can use plot_beamforming! to create plots."
                ) plot_beamforming!(arf)
            end

            @testset "plot_vespagram[!]" begin
                @test_throws (
                    "Note: If you have not already, you must first load a Makie " *
                    "backend such as GLMakie, CairoMakie or WGLMakie before" *
                    "you can use plot_vespagram! to create plots."
                ) plot_vespagram(arf)
                @test_throws (
                    "Note: If you have not already, you must first load a Makie " *
                    "backend such as GLMakie, CairoMakie or WGLMakie before" *
                    "you can use plot_vespagram! to create plots."
                ) plot_vespagram!(arf)
            end
        end
    end
end