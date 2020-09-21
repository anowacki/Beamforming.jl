# Recipes for beamforming grids, vespagrams and ARFs

using RecipesBase: RecipesBase, @recipe, @series

#=
    Beamforming grids
=#
@recipe function f(bf::BeamformGrid;
        phases=[], model="iasp91", powscale=:linear,
        s_km=false,
        dazimuth=30,
        dslowness=_default_dslowness(bf, s_km),
    )
    # Slowness grid
    sx, sy = if s_km
        s_per_km.(bf.sx), s_per_km.(bf.sy)
    else
        bf.sx, bf.sy
    end

    # Plotting defaults
    xlims --> extrema(sx)
    ylims --> extrema(sy)
    aspect_ratio --> :equal
    framestyle --> :box
    title --> "Normalised beam power ($(powscale))"
    fontfamily --> "Helvetica"

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
        sx, sy, power
    end

    # Backazimuth and slowness lines
    @series begin
        linecolor --> :white
        linewidth --> 1
        linestyle --> :dash
        label := ""
        x, y = _guidelines(sx, sy, dslowness, dazimuth)
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
                s_km && (dtdd = s_per_km(dtdd))
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

#=
    ARFs
=#
@recipe function f(arf::AbstractArrayResponse;
        powscale=:linear, s_km=false,
        dslowness=_default_dslowness(arf, s_km), # s_km ? 0.2 : 2,
        dazimuth=30
    )
    sx, sy = s_km ? (arf.sx, arf.sy) : s_per_degree.((arf.sx, arf.sy))
    xlabel, ylabel = _axis_labels(s_km)
    xlims --> extrema(sx)
    ylims --> extrema(sy)
    aspect_ratio --> :equal
    framestyle --> :box
    fontfamily --> "Helvetica"
    xguide --> xlabel
    yguide --> ylabel
    f1, f2 = round.(extrema(arf.freqs), sigdigits=2)
    freqlabel = length(arf.freqs) == 1 ? "$f1" : "$f1 - $f2"
    title --> "Array response ($freqlabel Hz)"
    colorbar_title --> "Beam power ($powscale)"
    # Power
    @series begin
        seriescolor --> :lightrainbow
        seriestype := :heatmap
        sx, sy, _scaled_power(arf, powscale, norm=true)
    end
    # Azimuth and slowness contours
    @series begin
        primary := false
        linecolor --> :white
        linewidth --> 1
        linestyle --> :dash
        x, y = _guidelines(sx, sy, dslowness, dazimuth)
        x, y
    end
end

#=
    Vespagrams
=#
@recipe function f(vespa::VespaGrid; phases=[], model="iasp91", normalise=true)
    xlims --> extrema(vespa.time)
    ylims --> extrema(vespa.slow)
    framestyle --> :box
    fontfamily --> "Helvetica"
    maxamp = maximum(abs, vespa.power)
    clims --> (normalise ? (-1, 1) : (-maxamp, maxamp))
    grid = normalise ? vespa.power' ./ maxamp : vespa.power'
    xguide --> "Time / s"
    yguide --> "Slowness / s/°"
    plot_title = (normalise ? "Normalised v" : "V") * "espagram " *
        (vespa.envelope ? "envelope" : "amplitude")
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
Return a scaled beam power correctly transposed for plotting.

`powscale` can be:
- `:linear`: No transformation applied.
- `:log10`: Base-10 log of power is taken.
- `:dB`: Power converted to dB = 10log₁₀(power).

If `norm` is `true`, power is scaled such that the maximum absolute power is
unity.
"""
function _scaled_power(bf::Union{BeamformGrid,AbstractArrayResponse}, powscale; norm=true)
    scale = norm ? maximum(abs, bf.power) : one(eltype(bf.power))
    if powscale == :linear
        bf.power'./scale
    elseif powscale == :log10
        log10.(bf.power'./scale)
    elseif powscale == :dB
        DSP.pow2db.(bf.power'./scale)
    else
        throw(ArgumentError("unrecognised plot powscale '$powscale'"))
    end
end

"""
Create paths for lines of constant slowness and constant azimuth on a grid of
`sx` and `sy` slowness points.

Constant-azimuth lines are placed `dazimuth`° apart, whilst constant-slowness
lines are placed `dslowness` apart; units must match those of `sx` and `sy`.
"""
function _guidelines(sx, sy, dslowness, dazimuth)
    x = Vector{Float64}[]
    y = similar(x)
    pmax = √2*maximum(abs, [sx; sy])
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

"Return a pair of axis labels with units s/km if `s_km` is `true`, and s/° otherwise."
function _axis_labels(s_km)
    unit = s_km ? "s/km" : "s/°"
    "East slowness ($unit)", "North slowness ($unit)"
end

"""
Pick a nice round number to give about 5 isoslowness contours for a beam power
plot.
"""
function _default_dslowness(grid::Union{AbstractArrayResponse,BeamformGrid}, s_km)
    maxslow = maximum(abs, [grid.sx; grid.sy])
    !s_km && (maxslow = s_per_degree(maxslow))
    dslow = maxslow/5
    n = floor(Int, log10(dslow))
    first_digit = floor(Int, dslow*10.0^-n)
    if first_digit == 0
        1/10.0^-n
    elseif first_digit in 1:2
        first_digit/10.0^-n
    elseif first_digit in 3:9
        5/10.0^-n
    end
end

"""
    plot(beamform_grid; kwargs...) -> ::Plot

Create a plot of a [`BeamformGrid`](@ref), as created by [`beamform`](@ref)
or [`crosscorrelation_beamform`](@ref).

# Keyword arguments:
- `phases = []`: Vector of phases to add to plot.  Before using this option,
  you must have imported `SeisTau`.
- `model = "iasp91"`: Model to use in finding predicted arrival slownesses.
- `powscale = :linear`: Scaling of beam power.  Choices are: `:linear`,
  `:log10` and `:dB`.
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness = 2`: Spacing of constant-slowness lines on plot,
  in s/° if `s_km` is `false` (the default), or s/km if `s_km` is `true`.
- `s_km = false`: Plot grid in s/km rather than s/° if `true`.

---

    plot(vespagram_grid; phases=[], model="iasp91") -> ::Plot

Create a plot of a vespagram computed by [`vespagram`](@ref).

Optionally provide a `Vector` of phases to plot.  Before using this option,
you must have imported `SeisTau`.  Specify the model for predicted
slownesses and travel times with `model`.

---

    plot(arf; kwargs...) -> ::Plot

Create a plot of an array response function as computed by [`array_response`](@ref)
or [`crosscorrelation_array_response`](@ref).

# Keyword arguments
- `powscale = :linear`: Scaling of array response power.  Choices are:
  `:linear`, `:log10` and `:dB`.
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness = 2`: Spacing of constant-slowness lines on plot,
  in s/° if `s_km` is `false` (the default), or s/km if `s_km` is `true`.
- `s_km = false`: Plot grid in s/km rather than s/° if `true`.

""" plot
