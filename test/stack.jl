using Test
using DelimitedFiles: readdlm

import Seis
using Beamforming

@testset "Stack" begin
    # These test data are generated with the script test/SAC/sac_stack.jl
    let t = Seis.sample_data(:array)
        # Time limits of stack, relative to A marker
        t1, t2 = -5, 10

        # Stacking weights
        weights = readdlm(joinpath(@__DIR__, "SAC", "stack_weights.txt"))

        for stack_type in ("linear", "nthroot", "phaseweight"), n in 1:2:5
            stack_file = joinpath(@__DIR__, "SAC", "$(stack_type)")
            if stack_type == "linear"
                n > 1 && continue
                nstring = ""
            else
                stack_file *= "-$n"
                nstring = "$n"
            end
            stack_file *= ".SHZ"

            sacstack = Seis.read_sac(stack_file)
            s = Beamforming.stack(t, (t1, t2), :A,
                method=Symbol(stack_type), n=n)

            @test Seis.trace(Seis.normalise(s)) ≈
                Seis.trace(Seis.normalise(sacstack)) rtol=0.2
        end
    end
end
