# Makie plots, activated by using a backend

"""
    plot_array_response!(ax::Makie.Axis, beamform_grid; kwargs...) -> ::Makie.Heatmap
    Makie.plot!(ax, arf; kwargs...) -> ::Makie.Heatmap

Add a plot of a [`ArrayResponse`](@ref) or [`ArrayResponseCrosscorrelation`],
as created by [`array_response`](@ref)
or [`crosscorrelation_array_response`](@ref) to an existing `Makie.Axis` object.

# Keyword arguments
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness`: Spacing of constant-slowness lines on plot,
  in s/° if `s_km` is `false` (the default), or s/km if `s_km` is `true`.
- `heatmap = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.heatmap` which creates the power heatmap.
- `normalize = true`: Perform further normalization (beyond that done in
  `beamform`) to plot normalised beam power.
- `powscale = :linear`: Scaling of array response power
  - If passed a `Symbol`, choices are: `:linear`, `:log10` and `:dB` and
    respectively scale the power linearly, by base-10 logarithm or as dB.
  - Otherwise, `powscale` is assumed to be a function which is applied to
    the *whole response power matrix* before plotting.  In other words, the
    function should take a single argument of type `Matrix`.
- `resize = true`: Whether or not to resize the axis to fit the array_response
  grid.
- `s_km = false`: Plot grid in s/km rather than s/° if `true`.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> smax, ds, fmin, fmax, df, s_per_deg = 4, 0.1, 0.1, 1, 0.1, true;

julia> arf = array_response(t, smax, ds, fmin, fmax, df, s_per_deg);

julia> fig = Figure(); ax = Axis(fig[1,1], aspect=DataAspect());

julia> plot_array_response!(ax, arf); fig
```
"""
function plot_array_response! end

"""
    plot_array_response(arf; kwargs...) -> (fig, ax, hm)::Makie.FigureAxisPlot
    Makie.plot(arf; kwargs...) -> (fig, ax, hm)::Makie.FigureAxisPlot

Create a plot of a [`ArrayResponse`](@ref) or [`ArrayResponseCrosscorrelation`](@ref),
as created by [`array_response`](@ref) or [`crosscorrelation_array_response`](@ref),
returning the figure `fig`, axis `ax` and heamap object `hm` as a `Makie.FigureAxisPlot`.

# Keyword arguments
- `axis = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.Axis` which creates the beam power axis.
- `colormap = :turbo`: Colormap for beam power heatmap.
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness`: Spacing of constant-slowness lines on plot,
  in s/° if `s_km` is `false` (the default), or s/km if `s_km` is `true`.
- `figure = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.Figure` which creates the figure.
- `heatmap = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.heatmap` which creates the power heatmap.
- `normalize = true`: Perform further normalization (beyond that done in
  `beamform`) to plot normalised beam power.
- `powscale = :linear`: Scaling of beam power
  - If passed a `Symbol`, choices are: `:linear`, `:log10` and `:dB` and
    respectively scale the power linearly, by base-10 logarithm or as dB.
  - Otherwise, `powscale` is assumed to be a function which is applied to
    the *whole beam power matrix* before plotting.  In other words, the
    function should take a single argument of type `Matrix`.
- `resize = true`: Whether or not to resize the figure to fit the plot.
- `s_km = false`: Plot grid in s/km rather than s/° if `true`.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> smax, ds, fmin, fmax, df, s_per_deg = 4, 0.1, 0.1, 1, 0.1, true;

julia> arf = array_response(t, smax, ds, fmin, fmax, df, s_per_deg);

julia> plot_array_response(arf)
```

---

    plot_array_response(::Makie.GridPosition, arf; kwargs...) -> (ax, hm)::Makie.AxisPlot
    Makie.plot(::Makie.GridPosition, arf; kwargs...) -> (ax, hm):Makie.AxisPlot

Add an array response power plot to an existing `Makie.Figure` by specifing
the figure grid position to use, returning a `Makie.AxisPlot` which can be
iterated to give the new `Makie.Axis` and `Makie.Heatmap` objects.

The `hm` object returned can be used with `Makie.Colorbar` to plot a colour
scale bar for the array response power.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> smax, ds, fmin, fmax, df, s_per_deg = 4, 0.1, 0.1, 1, 0.1, true;

julia> arf = crosscorrelation_array_response(t, smax, ds, fmin, fmax, df, s_per_deg);

julia> fig = Makie.Figure();

julia> ax, hm = plot(fig[1,1], arf); fig
```
"""
function plot_array_response end


"""
    plot_beamforming(beamform_grid; kwargs...) -> (fig, ax, hm)::Makie.FigureAxisPlot
    Makie.plot(beamform_grid::BeamformGrid; kwargs...) -> (fig, ax, hm)::Makie.FigureAxisPlot

Create a plot of a [`BeamformGrid`](@ref), as created by [`beamform`](@ref)
or [`crosscorrelation_beamform`](@ref), returning the figure `fig`,
axis `ax` and heatmap object `hm` as a `Makie.FigureAxisPlot`.

# Keyword arguments
- `axis = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.Axis` which creates the beam power axis.
- `colormap = :turbo`: Colormap for beam power heatmap.
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness`: Spacing of constant-slowness lines on plot,
  in s/° if `s_km` is `false` (the default), or s/km if `s_km` is `true`.
- `figure = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.Figure` which creates the figure.
- `heatmap = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.heatmap` which creates the power heatmap.
- `model = "iasp91"`: Model to use in finding predicted arrival slownesses.
- `normalize = true`: Perform further normalization (beyond that done in
  `beamform`) to plot normalised beam power.
- `phases = []`: Vector of phases to add to plot.  Before using this option,
  you must have imported `SeisTau`.
- `powscale = :linear`: Scaling of beam power
  - If passed a `Symbol`, choices are: `:linear`, `:log10` and `:dB` and
    respectively scale the power linearly, by base-10 logarithm or as dB.
  - Otherwise, `powscale` is assumed to be a function which is applied to
    the *whole beam power matrix* before plotting.  In other words, the
    function should take a single argument of type `Matrix`.
- `resize = true`: Whether or not to resize the figure to fit the plot.
- `s_km = false`: Plot grid in s/km rather than s/° if `true`.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> bf = beamform(t, 1100, 1150, 5, 0.1);

julia> plot_beamforming(bf)
```

---

    plot_beamforming(::Makie.GridPosition, beamform_grid; kwargs) -> (ax, hm)::Makie.AxisPlot
    Makie.plot(::Makie.GridPosition, beamform_grid::BeamformGrid; kwargs) -> (ax, hm)::Makie.AxisPlot

Add a beamforming power plot to an existing `Makie.Figure` by specifing
the figure grid position to use, returning a `Makie.AxisPlot` which can be
iterated to give the new `Makie.Axis` and `Makie.Heatmap` objects.

The `hm` object returned can be used with `Makie.Colorbar` to plot a colour
scale bar for the array response power.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> bf = beamform(t, 1100, 1150, 5, 0.1);

julia> fig = Makie.Figure();

julia> ax, hm = plot(fig[1,1], bf)
```
"""
function plot_beamforming end

"""
    plot_beamforming!(ax::Makie.Axis, beamform_grid; kwargs...) -> hm::Makie.Heatmap
    Makie.plot!(ax::Makie.Axis, beamform_grid::BeamformGrid; kwargs...) -> hm::Makie.Heatmap

Add a plot of a [`BeamformGrid`](@ref), as created by [`beamform`](@ref)
or [`crosscorrelation_beamform`](@ref) to an existing `Makie.Axis` object.

# Keyword arguments
- `dazimuth = 30`: Spacing of constant-azimuth lines on plot in °.
- `dslowness`: Spacing of constant-slowness lines on plot,
  in s/° if `s_km` is `false` (the default), or s/km if `s_km` is `true`.
- `heatmap = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.heatmap` which creates the power heatmap.
- `model = "iasp91"`: Model to use in finding predicted arrival slownesses.
- `normalize = true`: Perform further normalization (beyond that done in
  `beamform`) to plot normalised beam power.
- `phases = []`: Vector of phases to add to plot.  Before using this option,
  you must have imported `SeisTau`.
- `powscale = :linear`: Scaling of beam power
  - If passed a `Symbol`, choices are: `:linear`, `:log10` and `:dB` and
    respectively scale the power linearly, by base-10 logarithm or as dB.
  - Otherwise, `powscale` is assumed to be a function which is applied to
    the *whole beam power matrix* before plotting.  In other words, the
    function should take a single argument of type `Matrix`.
- `resize = true`: Whether or not to resize the axis to fit the beamforming
  grid.
- `s_km = false`: Plot grid in s/km rather than s/° if `true`.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> bf = beamform(t, 1100, 1150, 5, 0.1);

julia> fig = Figure(); ax = Axis(fig[1,1], aspect=DataAspect());

julia> plot_beamforming!(ax, bf)
"""
function plot_beamforming! end

"""
    plot_vespagram(vesp::VespaGrid; kwargs...) -> (fig, ax, hm)::Makie.FigureAxisPlot
    Makie.plot(vesp::VespaGrid; kwargs...) -> (fig, ax, hm)::Makie.FigureAxisPlot

Create a plot of a [`VespaGrid`](@ref) object, as created by [`vespagram`](@ref),
returning the figure `fig`, axis `ax` and heatmap object `hm` as a
`Makie.FigureAxisPlot`.

# Keyword arguments
- `heatmap = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.heatmap` which creates the power heatmap.
- `model = "iasp91"`: Model to use in finding predicted arrival slownesses.
- `normalize = true`: Perform further normalization (beyond that done in
  `vespagram`) to plot normalised vespagram amplitude.
- `phases = []`: Vector of phases to add to plot.  Before using this option,
  you must have imported `SeisTau`.
- `ampscale = :linear`: Scaling of beam amplitude
  - If passed a `Symbol`, choices are: `:linear`
  - Otherwise, `ampscale` is assumed to be a function which is applied to
    the *whole vespagram matrix* before plotting.  In other words, the
    function should take a single argument of type `Matrix`.
- `resize = true`: Whether or not to resize the axis to fit the vespagram
  grid.
- `s_km = false`: Plot slowness in s/km rather than s/° if `true`.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> t1, t2, s1, s2, ds = 1100, 1200, 1, 4, 0.025;

julia> vesp = vespagram(t, t1, t2, s1, s2, ds);

julia> fig, ax, hm = plot_vespagram(vesp)
```

---

    plot_vespagram(::Makie.GridPosition, vesp; kwargs) -> (ax, hm)::Makie.AxisPlot
    Makie.plot(::Makie.GridPosition, vesp; kwargs) -> (ax, hm)::Makie.AxisPlot

Add a vespagram plot to an existing `Makie.Figure` by specifing
the figure grid position to use, returning a `Makie.AxisPlot` which can be
iterated to give the new `Makie.Axis` and `Makie.Heatmap` objects.

The `hm` object returned can be used with `Makie.Colorbar` to plot a colour
scale bar for the vespagram amplitude.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> vesp = vespagram(t, 1100, 1200, 1, 4, 0.025);

julia> fig = Makie.Figure()

julia> ax, hm = plot(fig[1,1], bf)
```

"""
function plot_vespagram end

"""
    plot_vespagram!(ax::Makie.Axis, vesp; kwargs...) -> hm::Makie.Heatmap
    Makie.plot!(ax::Makie.Axis, vesp; kwargs...) -> hm::Makie.Heatmap

Add a plot of a [`VespaGrid`](@ref), as created by [`vespagram`](@ref)
to an existing `Makie.Axis` object.

The `hm` object returned can be used with `Makie.Colorbar` to plot a colour
scale bar for the vespagram amplitude.

# Keyword arguments
- `heatmap = ()`: Keyword arguments (as a named tuple or dictionary) which
  are passed to the call to `Makie.heatmap` which creates the power heatmap.
- `model = "iasp91"`: Model to use in finding predicted arrival slownesses.
- `normalize = true`: Perform further normalization (beyond that done in
  `beamform`) to plot normalised beam power.
- `phases = []`: Vector of phases to add to plot.  Before using this option,
  you must have imported `SeisTau`.
- `ampscale = :linear`: Scaling of beam power
  - If passed a `Symbol`, choices are: `:linear`.
  - Otherwise, `powscale` is assumed to be a function which is applied to
    the *whole vespagram amplitude matrix* before plotting.  In other words, the
    function should take a single argument of type `Matrix`.
- `resize = true`: Whether or not to resize the axis to fit the beamforming
  grid.
- `s_km = false`: Plot slowness in s/km rather than s/° if `true`.

# Example
```
julia> using Seis, GLMakie, Beamforming

julia> t = sample_data(:array);

julia> vesp = vespagram(t, 1100, 1200, 1, 4, 0.025);

julia> fig = Figure()

julia> ax = Axis(fig[1,1]);

julia> hm = plot_vespagram!(ax, vesp)
"""
function plot_vespagram! end
