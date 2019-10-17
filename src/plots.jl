# Recipes for beamforming grids and ARFs

using RecipesBase: RecipesBase, @recipe, @series

@recipe function f(bf::BeamformGrid; phases=[], model="iasp91", powscale=:linear, dazimuth=30, dslowness=2)
    xlim --> extrema(bf.sx)
    ylim --> extrema(bf.sy)
    aspect_ratio --> :equal
    framestyle --> :box
    title --> "Beam power"

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
    @series begin
        seriestype := :scatter
        markercolor --> :white
        markersize --> 4
        markershape --> :square
        markerstrokecolor --> :black
        markerstrokewidth --> 2
        label := ""
        x = Float64[]
        y = Float64[]
        # Deal with single phase case
        phases = phases isa Union{Tuple,AbstractArray} ? phases : (phases,)
        if _hastraces(bf) !isempty(phases)
            evt = bf.traces[1].evt
            β = Geodesics.azimuth(bf.lon, bf.lat, evt.lon, evt.lat, true)
            for phase in phases
                arr = Main.SeisTau.travel_time(bf.traces[1], phase, model=model)
                isempty(arr) && continue
                dtdd = first(arr).dtdd
                push!(x, sind(β)*dtdd)
                push!(y, cosd(β)*dtdd)
            end
        end
        x, y
    end
end

@recipe function f(arf::ArrayResponse)
    xlim --> extrema(arf.sx)
    ylim --> extrema(arf.sy)
    aspect_ratio --> :equal
    framestyle --> :box
    @series begin
        seriescolor --> :lightrainbow
        seriestype := :heatmap
        arf.sx, arf.sy, arf.power'
    end
end

"""
    plot(beamform_grid; kwargs...) -> ::Plot

Create a plot of a [`BeamformGrid`](@ref), as created by [`beamform`](@ref).

    plot(arf; kwargs...) -> ::Plot

Create a plot of an 
""" plot