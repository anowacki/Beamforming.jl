# Reproduce Figure 2 of Ruigrok et al., J Seismol, 2017

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Beamforming, Seis, Plots

"Compute all unique pairs of stations"
station_pairs(stas) = [(stas[i], stas[j]) for i in 1:length(stas) for j in i+1:length(stas)]

# Plotting defaults
default(fontfamily="Helvetica")

# Three ARCES stations.
# Coordinates from https://www.fdsn.org/media/meetings/2015/NORSAR_Network_Report.pdf
ara1 = Station(sta="ARA1", lon=25.5071, lat=69.5363)
ara2 = Station(sta="ARA2", lon=25.5078, lat=69.5338)
arb2 = Station(sta="ARB2", lon=25.5134, lat=69.5357)

stas = [ara1, ara2, arb2]
n = length(stas)
t = [Trace(0, 1, 0) for _ in 1:n]
t.sta = stas

# Station cartesian coordinates and mean longitude and latitude
x, y, z, lon, lat = Beamforming.array_geometry(t)

# Cartesian stations
cstas = [CartStation(x=1000xx, y=1000yy, z=1000zz, sta=sta)
         for (xx, yy, zz, sta) in zip(x, y, z, stas.sta)]


# Beamforming parameters for all
sx = sy = -1:0.01:1
smax = max(maximum.((sx, sy))...)
ds = step(sx)
ps = 0
θs = 0

# Plotting arguments
cpt = cgrad([:white, :blue, :cyan, :yellow, :orange, :red, :darkred])
kwargs = (c=cpt, clim=(-10,0), s_km=true, powscale=:dB)

# Single 5 Hz wave, absolute power
f = 5
df = 1
P_5 = Beamforming.crosscorrelation_array_response(stas, smax, ds, f, f, df, powf=abs)
arf_5 = plot(P_5; kwargs...)

# 3 – 7 Hz, absolute power
f1 = 3
f2 = 7
df = 0.1
P_3_7 = Beamforming.crosscorrelation_array_response(stas, smax, ds, f1, f2, df, powf=abs)
arf_3_7 = plot(P_3_7; kwargs...)

# 3 – 7 Hz, squared power
P_3_7² = Beamforming.crosscorrelation_array_response(stas, smax, ds, f1, f2, df, powf=abs2)
arf_3_7² = plot(P_3_7²; kwargs...)

# Map of array
pairs = station_pairs(cstas)
map = plot([first.(pairs).x last.(pairs).x]', [first.(pairs).y last.(pairs).y]',
    lc=:black, grid=false, framestyle=:box, aspect_ratio=:equal, primary=false,
    xlabel="Easting (m)", ylabel="Northing (m)")
scatter!(map, cstas.x, cstas.y, m=(:dtriangle, :green), label="")
annotate!(map, cstas.x .+ 30, cstas.y, text.(cstas.sta, 8))

# Plot of azimuths and distances
X, Y = [], []
for (θ, h2) in zip(Beamforming._azimuth_distance_ccbf(x, y)...)
    h = h2/2
    push!(X, [0, h*cosd(θ), NaN, 0, -h*cosd(θ)])
    push!(Y, [0, h*sind(θ), NaN, 0, -h*sind(θ)])
end
azmap = plot(X, Y, primary=false, aspect_ratio=:equal, framestyle=:box, grid=false, arrow=:arrow,
    xlabel="East distance (km)", ylabel="North distance (km)")

# Colorbar
cbar = heatmap([0, 2], -10:0.01:0, [collect(-10:0.01:0) collect(-10:0.01:0)], c=cpt, xticks=nothing, cbar=false,
    aspect_ratio=8, xlim=(0,2), ylabel="Beam power (dB)")

# All plots assembled together
allplots = plot(map, azmap, cbar, arf_5, arf_3_7, arf_3_7²,
    layout=(2,3), size=(900,500), label="", cbar=false, titlefontsize=12)

# Save plot if run as a script, otherwise return (and show)
if !isinteractive() && abspath(@__FILE__) == abspath(PROGRAM_FILE)
    savefig(allplots, joinpath(@__DIR__, "Ruigrok_et_al_2017_Fig2.pdf"))
else
    allplots
end
