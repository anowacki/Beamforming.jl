# Default plot overloads
Makie.plot(vesp::VespaGrid; kwargs...) = Beamforming.plot_vespagram(vesp; kwargs...)
Makie.plot(gp::Makie.GridPosition, vesp::VespaGrid; kwargs...) = Beamforming.plot_vespagram(gp, vesp; kwargs...)
Makie.plot!(ax::Makie.Axis, vesp::VespaGrid; kwargs...) = Beamforming.plot_vespagram!(ax, vesp; kwargs...)

function Beamforming.plot_vespagram!(ax::Makie.Axis, vesp::VespaGrid;
    ampscale=:linear,
    colormap=_vespagram_colormap(vesp),
    heatmap=(),
    model="iasp91",
    normalize=true,
    phases=[],
    resize=true,
    s_km=false,
    colorrange=_vespagram_colorrange(vesp, normalize),
)
    # Function to transform to correct units
    slowness_trans = if s_km
        Beamforming.s_per_km
    else
        identity
    end

    # Convert amplitude to appropriate scale
    amplitude = _scale_power(vesp.power, normalize, ampscale)

    slow = slowness_trans.(vesp.slow)

    # Vespagram amplitude
    hm = Makie.heatmap!(
        ax, vesp.time, slow, amplitude;
        colormap,
        colorrange,
        heatmap...
    )

    # Predicted 1D phase arrivals, if event information is present
    phases = phases isa Union{Tuple,AbstractArray} ? phases : (phases,)
    phases_time, phases_slow = Float64[], Float64[]
    phase_annotations = String[]

    if !isempty(phases)
        if ! Beamforming._hastraces(vesp)
            throw(ArgumentError(
                "vespagram object does not contain traces, " *
                "which are needed to add phase arrivals"
            ))
        end

        evt = vesp.traces[1].evt
        β = Geodesics.azimuth(vesp.lon, vesp.lat, evt.lon, evt.lat, true)
        for phase in phases
            # FIXME: Find a better way to conditionally calculate this depening
            # on whether SeisTau is loaded
            try
                arr = Main.SeisTau.travel_time(first(vesp.traces), phase, model=model)
                isempty(arr) && continue
                dtdd = first(arr).dtdd
                time = first(arr).time
                s_km && (dtdd = Beamforming.s_per_km(dtdd))
                push!(phases_time, time)
                push!(phases_slow, dtdd)
                push!(phase_annotations, phase)
            catch err
                if err isa UndefVarError && err.var == :SeisTau
                    @warn("SeisTau is not loaded; no phase arrivals plotted.  " *
                          "To add phase arrivals, do `using SeisTau`, then " *
                          "try this plotting command again.")
                    break
                else
                    rethrow(err)
                end
            end
        end

        Makie.scatter!(ax, phases_time, phases_slow;
            _DEFAULT_SCATTER_KWARGS...
        )

        Makie.text!(ax, phases_time, phases_slow; text=phase_annotations, color=:black)
    end

    if resize
        Makie.xlims!(ax, extrema(vesp.time)...)
        Makie.ylims!(ax, extrema(slow)...)
    end

    hm
end

function Beamforming.plot_vespagram(gp::Makie.GridPosition, vesp::VespaGrid;
    axis=(),
    s_km=false,
    kwargs...
)
    slowness_unit = s_km ? "s/km" : "s/°"
    ax = Makie.Axis(gp;
        xlabel="Time (s)",
        ylabel="Slowness ($slowness_unit)",
        axis...)
    hm = Beamforming.plot_vespagram!(ax, vesp; s_km, kwargs...)
    Makie.AxisPlot(ax, hm)
end

function Beamforming.plot_vespagram(vesp::VespaGrid;
    ampscale=:linear,
    axis=(),
    colormap=_vespagram_colormap(vesp),
    figure=(),
    model="iasp91",
    normalize=true,
    phases=[],
    resize=true,
    s_km=false,
    colorrange=_vespagram_colorrange(vesp, normalize),
    kwargs...
)
    slowness_unit = s_km ? "s/km" : "s/°"

    fig = Makie.Figure(; figure...)
    ax = Makie.Axis(fig[1,1];
        xlabel="Time (s)",
        ylabel="Slowness ($slowness_unit)",
        axis...
    )
    hm = Beamforming.plot_vespagram!(ax, vesp;
        colormap, colorrange, model, normalize, phases, ampscale, resize, s_km,
        kwargs...
    )

    amplitude_name = vesp.envelope ? "Envelope" : "Amplitude"
    if normalize
        amplitude_name = "Normalized " * lowercase(amplitude_name)
    end
    amplitude_units = _power_scale_unit(ampscale)

    Makie.Colorbar(fig[1,2], hm; label=amplitude_name*amplitude_units)

    if resize
        Makie.resize_to_layout!(fig)
    end

    Makie.FigureAxisPlot(fig, ax, hm)
end

function _vespagram_colorrange(vesp, normalize)
    maxval = normalize ? one(eltype(vesp.power)) : maximum(abs, vesp.power)
    vesp.envelope ? (zero(eltype(vesp.power)), maxval) : (-maxval, maxval)
end

_vespagram_colormap(vesp) = vesp.envelope ?
    _DEFAULT_LINEAR_COLORMAP :
    _DEFAULT_DIVERGING_COLORMAP
