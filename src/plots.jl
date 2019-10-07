# Recipes for beamforming grids and ARFs

using RecipesBase: RecipesBase, @recipe, @series

@recipe function f(bf::BeamformGrid)
	xlim --> extrema(bf.sx)
	ylim --> extrema(bf.sy)
	aspect_ratio --> :equal
	framestyle --> :box
	title --> "Beam power"
	@series begin
		seriescolor --> :lightrainbow
		seriestype := :heatmap
		bf.sx, bf.sy, bf.power'
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