# Script to test the use of cross-correlation beamforming

using Plots

using Beamforming
using Seis, Seis.Plot
using SeisTau

"Return a random selection of `n` items from `set`"
function random_selection(set, n)
    n > length(set) &&
        throw(ArgumentError("cannot select more than all the items in set"))
    out = Vector{eltype(set)}(undef, n)
    indices = collect(eachindex(set))
    for i in 1:n
        index = rand(1:length(indices))
        out[i] = set[indices[index]]
        deleteat!(indices, index)
    end
    out
end

# Plotting defaults
default(fontfamily="Helvetica")

# Choose case to look at
case = :synth
case = :real_array
case = :real_local

# Synthetic case
if case == :synth
    # Synthetic slowness of arrival
    sx, sy = 1, 2
    t = Beamforming.Synth.synthetic_arrival(sx, sy, nsta=8, lon=2rand(8), lat=2rand(8),
        distance=90)

    # Beamforming parameters
    t1, t2 = 40, 60
    window_time = 1
    s_max = 3
    ds = 0.05

    # Plotting parameters
    phases = []
# Real data example
elseif case == :real_array
    # Find a random selection of 6 traces from the UK network
    tt = sample_data(:array)
    t = random_selection(tt, 6)
    remove_trend!.(t)
    taper!.(t)
    bandpass!.(t, 0.05, 2, twopass=true)

    # Beamforming parameters
    t1, t2 = 1110, 1150
    window_time = 0.5
    s_max = 4
    ds = 0.05
    sx = sy = NaN

    phases = ["PKiKP", "PKIKP"]
elseif case == :real_local
    # Vertical stations only
    t = filter(x -> endswith(x.sta.sta, r"z"i), sample_data(:local))
    for tt in t
        remove_trend!(tt)
        taper!(tt)
        bandpass!(tt, 0.05, 1, twopass=true)
    end
    t1, t2 = 1, 10
    window_time = 2
    s_max = 20
    ds = 0.2
    sx = sy = NaN

    phases = ["p", "s"]
else
    error("choose a case to look at")
end

# Generic BF options for both CCBF and BF
method = :linear
n = 3

# Cross-correlation beamforming
bf = Beamforming.crosscorrelation_beamform(t, t1, t2,
    -s_max, s_max, -s_max, s_max, ds, t=window_time, method=method, n=n)

# Conventional beamforming for comparison
bf′ = Beamforming.beamform(t, t1, t2, s_max, ds, method=method, n=n)

# Plot parameters
colors = cgrad([:white, :blue, :green, :yellow, :orange, :red])
bf_plot_kwargs = (powscale=:dB, c=colors, clim=(-10,0),
    xlabel="X-slowness (s/°)", ylabel="Y-slowness (s/°)",
    phases=phases
    )

# Assemble plots
begin
    # Input trace plot
    p_in = section(bf.traces, title="Input traces")
    vspan!(p_in, [t1, t2], label="", fillcolor=:blue, fillalpha=0.1)

    # Cross-correlations plots
    p_cc_d = section(bf.corrs, distance_km, ylabel="Inter-station distance / km",
        title="Cross-correlations")
    vspan!(p_cc_d, [-window_time, window_time], label="", fillcolor=:blue, fillalpha=0.1)

    p_cc_az = section(bf.corrs, azimuth, ylabel="Inter-station azimuth / °",
        title="Cross-correlations")
    vspan!(p_cc_az, [-window_time, window_time], label="", fillcolor=:blue, fillalpha=0.1)

    # Map of stations with pair lines plotted
    p_map = plot([[x.evt.lon, x.sta.lon] for x in bf.corrs],
        [[x.evt.lat, x.sta.lat] for x in bf.corrs], linecolor=:black, primary=false,
        framestyle=:box, grid=false)
    scatter!(p_map, bf.traces.sta.lon, bf.traces.sta.lat, markershape=:dtriangle,
        markercolor=:blue, aspect_ratio=:equal, xlabel="Longitude", ylabel="Latitude",
        title="Array map", label="")

    # CC beamforming plot
    p_bf = plot(bf; title="CC beam power (dB)", bf_plot_kwargs...)
    if case == :synth
        scatter!(p_bf, [bf.sx_max[1] sx], [bf.sy_max[1] sy],
            label=["Recovered" "True"], markercolor=[:black :white], legend=:bottomleft)
    end

    # BF plot
    p_bf′ = plot(bf′; title="Beam power (dB)", bf_plot_kwargs...)
    if case == :synth
        scatter!(p_bf′, [sx], [sy], label="", markercolor=:black)
    end

    # Put all plots together (and display)
    plot(p_in, p_cc_d, p_cc_az, p_map, p_bf, p_bf′, layout=(2, 3), size=(1200,800))
end
