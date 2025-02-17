# Common functionality for the Makie extension

"Default color map used for heatmap-style plots where the colour scale is one-sided"
const _DEFAULT_LINEAR_COLORMAP = :turbo

"Default color map used for heatmap-style plots where the colour scale is centred on 0"
const _DEFAULT_DIVERGING_COLORMAP = :RdBu

"Default scatterplot settings"
const _DEFAULT_SCATTER_KWARGS = (color=:orange, strokecolor=:white, strokewidth=1.5)

"""
    _scale_power(power, normalize, powscale) -> power′

Create a new matrix `power′`, which is a copy of `power` where the values
are normalised to maximum 1 if `normalize` is `true`, the scaled accoring
to `powscale`.
"""
function _scale_power(power, normalize, powscale)
    # Normalisation
    scale = normalize ? maximum(power) : one(eltype(power))

    power = if powscale === :linear
        power./scale
    elseif powscale === :log10
        log10.(power./scale)
    elseif powscale === :dB
        10 .* log10.(power./scale)
    elseif powscale isa Symbol
        throw(ArgumentError("unrecognised powscale option '$powscale'"))
    else
        try
            powscale(power)
        catch err
            throw(ArgumentError("cannot apply provided powscale function as an argument:\n$err"))
        end
    end

    power
end

"""
    _power_scale_unit(powscale) -> name::String

Choose a name for the power scaling function, e.g. with which to annotate
a colour scale bar.
"""
function _power_scale_unit(powscale)
    if powscale === :linear
        ""
    elseif powscale === :log10
        " (log10)"
    elseif powscale === :dB
        " (dB)"
    elseif powscale isa Symbol
        @warn("Unhandled power scale '$powscale'")
        ""
    else
        ""
    end
end

"""
    _plot_polar_grid!(ax, sx, sy, dslowness, dazimuth; lines=(color=:white, linewidth=1.5))

Add lines to the existing `Makie.Axis` `ax` of constant slowness and constant
azimuth, spaced respectively `dslowness` apart (units same as `sx` and `sy`)
and `dazimuth`° apart.

`sx` and `sy` are the grid of slowness values being plotted.
"""
function _plot_polar_grid!(ax, sx, sy, dslowness, dazimuth;
    lines=(color=:white, linewidth=1.5)
)
    max_abs_s = √2*maximum(
        x -> hypot(x...),
        Iterators.product(extrema(sx), extrema(sy))
    )
    slowness_ring_values = if isnothing(dslowness)
        Makie.get_tickvalues(Makie.WilkinsonTicks(5; k_min=4), 0, max_abs_s)
    else
        0:dslowness:(cld(max_abs_s, dslowness)*dslowness)
    end
    dθ = 1 # degree
    θs = 0:dθ:360
    for slowness in slowness_ring_values
        Makie.lines!(
            ax, slowness.*sind.(θs), slowness.*cosd.(θs); lines...
        )
    end

    for azi in 0:dazimuth:(360 - dazimuth)
        Makie.lines!(
            ax, [0, max_abs_s].*sind(azi), [0, max_abs_s].*cosd(azi);
            lines...
        )
    end

end
