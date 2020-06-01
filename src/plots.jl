# Recipes for beamforming grids, vespagrams and ARFs

using RecipesBase: RecipesBase, @recipe, @series

@recipe function f(bf::BeamformGrid; phases=[], model="iasp91", powscale=:linear, dazimuth=30, dslowness=2)
    xlims --> extrema(bf.sx)
    ylims --> extrema(bf.sy)
    aspect_ratio --> :equal
    framestyle --> :box
    title --> "Normalised beam power"

    # Power heatmap
    @series begin
        seriescolor --> :lightrainbow
        seriestype := :heatmap
        power = if powscale == :linear
            bf.power'
        elseif powscale == :log10
            log10.(bf.power')
        elseif powscale == :dB
            DSP.pow2db.(bf.power')
        else
            @error("unrecognised plot powscale '$powscale'")
        end
        bf.sx, bf.sy, power
    end

    # Backazimuth and slowness lines
    @series begin
        linecolor --> :white
        linewidth --> 1
        linestyle --> :dash
        label := ""
        x = Vector{Float64}[]
        y = similar(x)
        pmax = √2*maximum([bf.sx; bf.sy])
        azimuths = range(0, stop=360, step=dazimuth)
        slownesses = range(0, stop=pmax, step=dslowness)
        # Constant azimuth lines
        for azimuth in azimuths
            push!(x, [0, sind(azimuth)*pmax])
            push!(y, [0, cosd(azimuth)*pmax])
        end
        # Constant slowness circles
        for slowness in slownesses
            azis = 0:360
            push!(x, sind.(azis)*slowness)
            push!(y, cosd.(azis)*slowness)
        end
        x, y
    end

    # 1D phase annotations, if event information is present
    phases = phases isa Union{Tuple,AbstractArray} ? phases : (phases,)
    phases_x, phases_y = Float64[], Float64[]
    phase_annotations = String[]
    if !isempty(phases)
        ! _hastraces(bf) && error("beamforming object does not contain traces, " *
                                 "which are needed to add phase arrivals")
        evt = bf.traces[1].evt
        β = Geodesics.azimuth(bf.lon, bf.lat, evt.lon, evt.lat, true)
        for phase in phases
            # FIXME: Find a better way to conditionally calculate this depening
            # on whether SeisTau is loaded
            try
                arr = Main.SeisTau.travel_time(first(bf.traces), phase, model=model)
                isempty(arr) && continue
                dtdd = first(arr).dtdd
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
    end

    # Phase location
    @series begin
        seriestype := :scatter
        markercolor --> :white
        markersize --> 4
        markershape --> :square
        markerstrokecolor --> :black
        markerstrokewidth --> 2
        label := ""
        phases_x, phases_y
    end

    # Phase label
    @series begin
        seriestype := :scatter
        primary := false
        markersize := 0
        label := ""
        # FIXME: Find a better way to set font properties
        series_annotations := Main.Plots.series_annotations(phase_annotations,
            Main.Plots.font("Helvetica", 9, :left, :bottom))
        phases_x, phases_y
    end
end

@recipe function f(arf::ArrayResponse)
    xlims --> extrema(arf.sx)
    ylims --> extrema(arf.sy)
    aspect_ratio --> :equal
    framestyle --> :box
    @series begin
        seriescolor --> :lightrainbow
        seriestype := :heatmap
        arf.sx, arf.sy, arf.power'
    end
end

@recipe function f(vespa::VespaGrid; phases=[], model="iasp91", normalise=true)
    xlims --> extrema(vespa.time)
    ylims --> extrema(vespa.slow)
    framestyle --> :box
    maxamp = maximum(abs, vespa.power)
    clims --> (normalise ? (-1, 1) : (-maxamp, maxamp))
    grid = normalise ? vespa.power' ./ maxamp : vespa.power'
    xlabel --> "Time / s"
    ylabel --> "Slowness / s/°"
    plot_title = (normalise ? "Normalised v" : "V") * "espagram " *
        (vespa. envelope ? "envelope" : "amplitude")
    title --> (vespa.envelope ? "Vespagram envelope" : "Vespagram amplitude")
    @series begin
        seriestype := :heatmap
        seriescolor --> (vespa.envelope ? :lightrainbow : :RdBu)
        vespa.time, vespa.slow, grid
    end
    # 1D phase annotations, if event information is present
    phases = phases isa Union{Tuple,AbstractArray} ? phases : (phases,)
    phases_x, phases_y = Float64[], Float64[]
    phase_annotations = String[]
    if !isempty(phases)
        ! _hastraces(vespa) && error("beamforming object does not contain traces, " *
                                     "which are needed to add phase arrivals")
        evt = vespa.traces[1].evt
        for phase in phases
            # FIXME: Find a better way to conditionally calculate this depening
            # on whether SeisTau is loaded
            try
                arr = Main.SeisTau.travel_time(first(vespa.traces), phase, model=model)
                isempty(arr) && continue
                dtdd = first(arr).dtdd
                tt = first(arr).time
                push!(phases_x, tt)
                push!(phases_y, dtdd)
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
    end

    # Phase location
    @series begin
        seriestype := :scatter
        markercolor --> :white
        markersize --> 4
        markershape --> :square
        markerstrokecolor --> :black
        markerstrokewidth --> 2
        label := ""
        phases_x, phases_y
    end

    # Phase label
    @series begin
        seriestype := :scatter
        primary := false
        markersize := 0
        label := ""
        # FIXME: Find a better way to set font properties
        series_annotations := Main.Plots.series_annotations(phase_annotations,
            Main.Plots.font("Helvetica", 9, :left, :bottom))
        phases_x, phases_y
    end
end

"""
    plot(beamform_grid; kwargs...) -> ::Plot

Create a plot of a [`BeamformGrid`](@ref), as created by [`beamform`](@ref).

# Keyword arguments:
- `phases = []`: Vector of phases to add to plot.  Before using this option,
  you must have imported `SeisTau`.
- `model = "iasp91"`: Model to use in finding predicted arrival slownesses.
- `powscale = :linear`: Scaling of beam power.  Choices are: `:linear`,
  `:log10` and `:dB`.
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness` = 2: Spacing of constant-slowness lines on plot in s/°.

---

    plot(vespagram_grid; phases=[], model="iasp91") -> ::Plot

Create a plot of a vespagram computed by [`vespagram`](@ref).

Optionally provide a `Vector` of phases to plot.  Before using this option,
you must have imported `SeisTau`.  Specfiy the model for predicted
slownesses and travel times with `model`.

---

    plot(arf; kwargs...) -> ::Plot

Create a plot of an array response function as computed by [`array_response`](@ref).
""" plot
