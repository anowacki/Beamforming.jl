# Default plot overloads
Makie.plot(arf::AbstractArrayResponse; kwargs...) = Beamforming.plot_array_response(arf; kwargs...)
Makie.plot(gp::Makie.GridPosition, arf::AbstractArrayResponse; kwargs...) = Beamforming.plot_array_response(gp, arf; kwargs...)
Makie.plot!(ax::Makie.Axis, arf::AbstractArrayResponse; kwargs...) = Beamforming.plot_array_response!(ax, arf; kwargs...)

function Beamforming.plot_array_response!(ax::Makie.Axis, arf::AbstractArrayResponse;
    colormap=_DEFAULT_LINEAR_COLORMAP,
    colorrange=Makie.Automatic(),
    normalize=true,
    powscale=:linear,
    s_km=false,
    dazimuth=30,
    dslowness=nothing,
    heatmap=(),
    resize=true,
)
    # Function to transform to correct units
    slowness_trans = if s_km
        identity
    else
        Beamforming.s_per_degree
    end

    # Convert power to appropriate scale
    power = _scale_power(arf.power, normalize, powscale)

    sx, sy = slowness_trans.(arf.sx), slowness_trans.(arf.sy)

    # Beam power
    hm = Makie.heatmap!(
        ax, sx, sy, power;
        colormap,
        colorrange,
        heatmap...
    )

    # Slowness and azimuth lines
    _plot_polar_grid!(ax, sx, sy, dslowness, dazimuth)

    if resize
        Makie.xlims!(ax, extrema(sx)...)
        Makie.ylims!(ax, extrema(sy)...)
    end

    hm
end

function Beamforming.plot_array_response(gp::Makie.GridPosition, arf::AbstractArrayResponse;
    axis=(),
    s_km=false,
    kwargs...
)
    slowness_unit = s_km ? "s/km" : "s/°"
    ax = Makie.Axis(gp;
        aspect=Makie.DataAspect(),
        xlabel="East slowness ($slowness_unit)",
        ylabel="North slowness ($slowness_unit)",
        axis...)
    hm = Beamforming.plot_array_response!(ax, arf; s_km, kwargs...)
    Makie.AxisPlot(ax, hm)
end

function Beamforming.plot_array_response(arf::AbstractArrayResponse;
    axis=(),
    colormap=_DEFAULT_LINEAR_COLORMAP,
    colorrange=Makie.Automatic(),
    figure=(),
    normalize=true,
    powscale=:linear,
    resize=true,
    s_km=false,
    kwargs...
)
    slowness_unit = s_km ? "s/km" : "s/°"

    fig = Makie.Figure(; figure...)
    ax = Makie.Axis(fig[1,1];
        aspect=Makie.DataAspect(),
        xlabel="East slowness ($slowness_unit)",
        ylabel="North slowness ($slowness_unit)",
        axis...
    )
    hm = Beamforming.plot_array_response!(ax, arf;
        colormap, colorrange, normalize, powscale, resize, s_km, kwargs...
    )

    power_name = normalize ? "Normalized beam power" : "Beam power"
    power_units = _power_scale_unit(powscale)

    Makie.Colorbar(fig[1,2], hm; label=power_name*power_units)

    if resize
        Makie.resize_to_layout!(fig)
    end

    Makie.FigureAxisPlot(fig, ax, hm)
end
