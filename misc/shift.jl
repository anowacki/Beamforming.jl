import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise

using Beamforming
using Plots
import Seis, Seis.Plot
using BenchmarkTools
import FFTW

"""
    synthetic_arrival(sx, sy) -> ::Vector{Seis.Trace}

Calculate synthetic traces with an arrival at a set slowness (`sx`, `sy`)
s/°
"""
function synthetic_arrival(sx, sy)
    # Synthetic arrival
    data = exp.(-(-50:0.01:50).^2/0.01)
    # Grid of stations
    delta = 0.01
    nsta = 20
    s = [Seis.Trace(0, delta, data) for _ in 1:nsta]
    s.sta.lon = lon = 5rand(nsta) - 2.5
    s.sta.lat = lat = 5rand(nsta) - 2.5
    mlon, mlat = mean(lon), mean(lat)
    # Shift arrivals to correct time
    for (i, ss) in enumerate(s)
        dist = Beamforming.delta(ss.sta.lon, ss.sta.lat, mlon, mlat, true)
        az = Beamforming.azimuth(mlon, mlat, ss.sat.lon, ss.sta.lat, true)
        dx = dist*sind(az) # °
        dy = dist*cosd(az)
        dt = sx*dx + sy*dy # s
        ss.t = circshift(ss.t, -round(Int, dt/delta))
    end
    s
end

function shift!(s::Seis.Trace, Δ)
    s.t .= circshift(s.t, round(Int, Δ/s.delta))
    s
end

shift!(s::AbstractArray{<:Seis.Trace}, Δ) = shift!.(s, Δ)

"Cyclically shift a SAC trace backward by Δ s in time"
function fdshift!(s::Seis.Trace{T}, Δ) where T
    S = FFTW.fft(Seis.trace(s))
    δ = s.delta
    nf = Int32(length(S))
    coeff = 2π*im*Δ/(δ * 2*(nf - 1))
    @inbounds for i in 1:nf
        S[i] = S[i]*exp(coeff*(i-1))
    end
    s.t .= real.(FFTW.ifft(S))
    s
end

function fdshift!(s::AbstractArray{<:Seis.Trace}, Δ::AbstractArray)
    @assert length(s) == length(Δ)
    S = FFTW.fft(hcat(Seis.trace.(s)...), 1)
    δ = s[1].delta
    nf, nsta = size(S)
    coeff = Complex{Float32}(2π*im/(δ*2*(nf-1)))
    @inbounds for j in 1:nsta, i in 1:nf
        S[i,j] = S[i,j]*exp(-coeff*Float32(Δ[j])*(i - 1))
    end
    FFTW.ifft!(S, 1)
    @inbounds for j in 1:nsta
        s[j].t .= real.(S[:,j])
    end
    s
end

# using ArrayFire
# function fdshift_af!(s::SACtr, Δ)
#     δ = s[:delta]
#     S = AFArray(Complex{Float32}.(s.t))
#     fft!(S)
#     nf = length(S)
#     coeffs = AFArray(Complex{Float32}.(2π.*im.*Δ./δ.*collect(0:(nf-1))./(2*(nf-1))))
#     S .= S.*exp(coeffs)
#     ifft!(S)
#     s[:t] = real(S)
#     s
# end
#
# function fdshift_af!(s::AbstractArray{SACtr}, Δ::AbstractArray)
#     @assert length(s) == length(Δ)
#     @assert all(s[:npts] .== s[1][:npts])
#     @assert all(isapprox.(s[:b], s[1][:b], atol=s[1][:delta]/4))
#     δ = s[1][:delta]
#     nsta = length(s)
#     npts = s[1][:npts]
#     np = nextpow2(s[1][:npts])
#     nf = np
#     S = zeros(Complex{Float32}, np, nsta)
#     S[1:npts,:] .= hcat(s[:t]...)
#     S = AFArray(S)
#     fft!(S)
#     coeffs = AFArray(Complex{Float32}.(2π.*im.*Δ'./δ.*collect(0:(nf-1))./(2*(nf-1))))
#     S = S.*exp(coeffs)
#     S = Array(real(ifft!(S, 1)))
#     @inbounds for i in 1:nsta
#         s[i].t[:] .= S[1:npts,i]
#     end
#     SAC.update_headers!(s)
# end

function beamformfd(s::AbstractArray{<:Seis.Trace}, t1, t2, sx1, sx2, sy1, sy2, ds)
    nsta = length(s)
    s = Seis.cut.(s, t1, t2)
    npts = Seis.nsamples(s[1])
    @assert all(Seis.nsamples.(s) .== npts)
    δ = s[1].delta
    sxs, sys = sx1:ds:sx2, sy1:ds:sy2
    sx = Beamforming.s_per_km.(sxs)
    sy = Beamforming.s_per_km.(sys)
    slow = hcat(vec(sxs), vec(sys))
    nslow = size(slow, 1)
    pow = zeros(Float32, nslow)
    Δ = Beamforming.delay_times(s, sx, sy)
    spec = FFTW.fft(hcat(Seis.trace.(s)...), 1)
    nf = size(spec, 1)
    # Steering vectors
    steers = 2π*im*collect(0:(nf-1))./(δ*2*(nf - 1))
    steer = Array{Complex{Float32}}(undef, nf, nsta)
    shift_spec = similar(spec)
    for k in eachindex(sys), j in eachindex(sxs)
        for i in 1:nf
            shift_spec[i,:] .= spec[i,:].*exp.(steers[i].*Δ[:,j,k])
        end
        for ista in 1:nsta
            temp = zero(eltype(shift_spec))
            for ifreq in 1:nf
                temp += shift_spec[ifreq,ista]
            end
            pow[j,k] += abs(temp^2)
        end
    end
    sxs, sys, pow
end

# function beamform_af(s::AbstractArray{SACtr}, Δ::AbstractArray, slow)
#     @assert length(s) == length(Δ)
#     @assert all(s[:npts] .== s[1][:npts])
#     @assert all(isapprox.(s[:b], s[1][:b], atol=s[1][:delta]/4))
#     δ = s[1][:delta]
#     nsta = length(s)
#     npts = s[1][:npts]
#     np = nextpow2(s[1][:npts])
#     nf = np
#     S = zeros(Complex{Float32}, np, nsta)
#     S[1:npts,:] .= hcat(s[:t]...)
#     S = AFArray(S)
# end

s = Seis.sample_data(:array)
Seis.cut!.(s, maximum(s.b), minimum(Seis.endtime.(s)))
Seis.cut!.(s, s.b, s.b.+s[1].delta.*(Seis.nsamples(s[1])-1))
# @btime fdshift!($s, $(-s[:a]))
# @btime fdshift_af!($s, $(-s[:a]))
@time s′ = shift!(deepcopy(s), -s.picks.A.time)
@time s″ = fdshift!(deepcopy(s), -2s.picks.A.time)

@time beamformfd(s, 1050, 1200, -3, 3, -3, 3, 0.5)
@time Beamforming.beamform(s, 1050, 1200, -3, 3, -3, 3, 0.5)

# Plots.plot([Seis.Plot.section(ss, align=ss.b) for ss in (s, s′, s″)]...)
