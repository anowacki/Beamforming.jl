# Compare stacking using circular and plane wave approximations
using Beamforming, Seis, TauPy, Plots
import Printf

function synthetic_prem_arrivals(β, Δ, phase;
		aperture=20, nsta=20, depth=0, pad=100, freq=0.1, delta=1/20freq)
	lon = aperture.*(rand(nsta) .- 0.5)
	lat = aperture.*(rand(nsta) .- 0.5)
	x, y, z, mean_lon, mean_lat = Beamforming.array_geometry(lon, lat, zeros(nsta))
	# Event
	evt = Seis.Event()
	evt.lon, evt.lat = Beamforming.geodesic_endpoint(mean_lon, mean_lat, β, Δ)
	evt.dep = depth
	# Distances
	Δs = Beamforming.delta.(lon, lat, evt.lon, evt.lat)
	# Get arrival times
	arr = try
		first.(TauPy.travel_time.(evt.dep, Δs, phase))
	catch err
		error("cannot find arrivals for $phase")
	end
	times = [a.time for a in arr]
	t1, t2 = extrema(times) .+ (-pad, pad)
	npts = round(Int, (t2 - t1)/delta) + 1
	# Traces
	s = [Seis.Synth.gaussian1(t1, delta, npts, freq, time) for time in times]
	s.sta.lon = lon
	s.sta.lat = lat
	s.sta.elev = 0
	s.sta.net = "AN"
	s.sta.sta = [Printf.@sprintf("A%01d", i) for i in 1:length(s)]
	s.sta.cha = "BXR"
	s.evt = evt
	s
end

s = synthetic_prem_arrivals(45, 120, "SKS", pad=500, freq=0.2, aperture=10)
arf = array_response(s, 1.5, 0.05, 0.01, 0.2, 0.01, true)
bfs = [beamform(s, 1500, 1600, 1, 4, 1, 4, 0.05, method=:linear, n=3,
	            wavefront=wave) for wave in (:plane, :circle)]
colormap = cgrad([:red,:orange,:yellow,:green,:blue,:white] |> reverse)
pl = plot(
	scatter(s.sta.lon, s.sta.lat, aspect_ratio=:equal, label="",
		framestyle=:box, grid=false, m=(:utriangle,:blue)),
	plot(arf, c=colormap),
	plot.(bfs, phases="SKS", c=colormap)...,
	title=["Array map" "ARF" "Plane" "Circle"],
	size=(700,700))