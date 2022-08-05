"""
    array_response(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df, s_deg=false) -> arf::ArrayResponse
    array_response(x, y, smax, ds, f1, f2, df, s_deg=false) -> arf::ArrayResponse
    array_response(traces::Array{<:Seis.Trace}, args...) -> arf::ArrayResponse
    array_response(stations::Array{<:Seis.Station}, args...) -> arf

Return the array response function `arf.power` (an [`ArrayResponse`](@ref)) at
slownesses with east and north components `arf.sx` and `arf.sy` respectively,
in s/km.  The ARF is computed between frequencies `f1` and `f2` Hz, spaced by `df` Hz.

Provide either explicit limits on slowness `(sx1, sx2)` and `(sy1, sy2)` in the
east and north components respectively, or the maximum slowness `smax`,
plus the slowness spacing `ds`.

`x` and `y` are the station coordinates in km; or an array of `Seis.Trace` `traces`
can be given.

Default input units of slowness are s/km.  Give `s_deg` as `true`
to use s/°.  Note that values in `arf` are always given in s/km.
"""
function array_response(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df, s_deg=false)
    if s_deg
        sx1, sx2, sy1, sy2, ds = s_per_km.((sx1, sx2, sy1, sy2, ds))
    end
    arf = ArrayResponse{eltype(x)}(x, y, sx1, sx2, sy1, sy2, ds, f1, f2, df)
    _compute_arf!(arf)
end

array_response(x, y, smax, ds, f1, f2, df, s_deg::Bool=false) =
    array_response(x, y, -smax, smax, -smax, smax, ds, f1, f2, df, s_deg)
array_response(s::AbstractArray{<:Seis.Station}, args...) =
    array_response(array_geometry(s)[1:2]..., args...)
array_response(s::TraceArray, args...) =
    array_response(array_geometry(s)[1:2]..., args...)

"""
    _compute_arf!(arf)

Do the actual ARF calculation for an `ArrayResponse`.
"""
function _compute_arf!(arf::ArrayResponse{T}) where T
    @inbounds for (j, slowy) in enumerate(arf.sy), (i, slowx) in enumerate(arf.sx)
        for f in arf.freqs
            ω = T(2)*π*f
            pow = zero(T)*im
            @simd for istat in eachindex(arf.x, arf.y)
                pow += cis(ω*(slowx*arf.x[istat] + slowy*arf.y[istat]))
            end
            arf.power[i,j] += abs2(pow)
        end
    end
    arf
end
