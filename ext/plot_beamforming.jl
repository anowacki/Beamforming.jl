
# Default plot overloads
Makie.plot(bf::BeamformGrid; kwargs...) = Beamforming.plot_beamforming(bf; kwargs...)
Makie.plot(gp::Makie.GridPosition, bf::BeamformGrid; kwargs...) = Beamforming.plot_beamforming(gp, bf; kwargs...)
Makie.plot!(ax::Makie.Axis, bf::BeamformGrid; kwargs...) = Beamforming.plot_beamforming!(ax, bf; kwargs...)

function Beamforming.plot_beamforming!(ax::Makie.Axis, bf::BeamformGrid;
    colormap=_DEFAULT_LINEAR_COLORMAP,
    colorrange=Makie.Automatic(),
    phases=[],
    model="iasp91",
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
        Beamforming.s_per_km
    else
        identity
    end

    # Convert power to appropriate scale
    power = _scale_power(bf.power, normalize, powscale)

    sx, sy = slowness_trans.(bf.sx), slowness_trans.(bf.sy)

    # Beam power
    hm = Makie.heatmap!(
        ax, sx, sy, power;
        colormap,
        colorrange,
        heatmap...
    )

    # Slowness and azimuth lines
    _plot_polar_grid!(ax, sx, sy, dslowness, dazimuth)

    # Predicted 1D phase arrivals, if event information is present
    phases = phases isa Union{Tuple,AbstractArray} ? phases : (phases,)
    phases_x, phases_y = Float64[], Float64[]
    phase_annotations = String[]

    if !isempty(phases)
        if ! Beamforming._hastraces(bf)
            throw(ArgumentError(
                "beamforming object does not contain traces, " *
                "which are needed to add phase arrivals"
            ))
        end

        evt = bf.traces[1].evt
        β = Geodesics.azimuth(bf.lon, bf.lat, evt.lon, evt.lat, true)
        for phase in phases
            # FIXME: Find a better way to conditionally calculate this depening
            # on whether SeisTau is loaded
            try
                arr = Main.SeisTau.travel_time(first(bf.traces), phase, model=model)
                isempty(arr) && continue
                dtdd = first(arr).dtdd
                s_km && (dtdd = Beamforming.s_per_km(dtdd))
                push!(phases_x, sind(β)*dtdd)
                push!(phases_y, cosd(β)*dtdd)
                push!(phase_annotations, phase)
            catch err
                if err isa UndefVarError && err.var == :SeisTau
                    @warn("SeisTau is not loaded; no phase arrivals plotted.  " *
                          "To add phase arrivals, do `using SeisTau`, then " *
                          "try this plotting command again.")
                else
                    rethrow(err)
                end
            end
        end

        Makie.scatter!(ax, phases_x, phases_y;
            _DEFAULT_SCATTER_KWARGS...
        )

        Makie.text!(ax, phases_x, phases_y; text=phase_annotations, color=:white)
    end

    if resize
        Makie.xlims!(ax, extrema(sx)...)
        Makie.ylims!(ax, extrema(sy)...)
    end

    hm
end

function Beamforming.plot_beamforming(gp::Makie.GridPosition, bf::BeamformGrid;
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
    hm = Beamforming.plot_beamforming!(ax, bf; kwargs...)
    Makie.AxisPlot(ax, hm)
end

function Beamforming.plot_beamforming(bf::BeamformGrid;
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
    hm = Beamforming.plot_beamforming!(ax, bf;
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
