# Beamforming

Beamforming is a Julia module for the array analysis of seismic
data, using [Seis.jl](https://github.com/anowacki/Seis.jl).

## Build status
[![Build Status](https://github.com/anowacki/Beamforming.jl/workflows/CI/badge.svg)](https://github.com/anowacki/Beamforming.jl/actions)
[![codecov](https://codecov.io/gh/anowacki/Beamforming.jl/branch/master/graph/badge.svg?token=d0ePcA1m54)](https://codecov.io/gh/anowacki/Beamforming.jl)



## Installation

```julia
julia> ] # press ']' to enter pkg mode

(v1.7) pkg> add https://github.com/anowacki/Geodesics.jl https://github.com/anowacki/Seis.jl https://github.com/anowacki/Beamforming.jl
```

Beamforming.jl requires Julia v1.6 or above.

## Usage overview

The main functions exported are:

- For traditional beamforming:
  - `array_response`: Compute the array response function for a set
    of stations
  - `beamform`: Compute the beam power across a grid of slowness points
    (as in f–k analysis)
  - `vespagram`: Compute a slowness vespagram
- For cross-correlation beamforming:
  - `crosscorrelation_array_response`
  - `crosscorrelation_beamform`

Docstrings are exhaustive, and can be consulted for usage of these
functions.  To bring up docstrings in the REPL, type `?` followed
by the name of the function and press return.

### Array response
The array response for a set of stations can be computed like so:
```julia
using Beamforming, Seis
t = sample_data(:array); # 60 UK stations
# Array response parameters
arf = array_response(
    t,   # Set of traces
    5,   # Maximum slowness in s/°
    0.1, # Slowness grid spacing
    0.1, # Minimum frequency in Hz
    1,   # Maximum frequency
    0.1, # Frequency spacing
    true # `true` for s/° slowness; `false` for s/km
)
using CairoMakie
plot(arf, powscale=:dB)   # If a Makie backend is installed
```
<img src="doc/images/array_response.png" width="600">

### Beamforming
Beamforming can be performed like so:
```julia
bf = beamform(
    t,
    1100, # Start time in s
    1160, # End time in s
    -1,   # Minimum east slowness in s/°
    0.5,  # Maximum east slowness
    0,    # Minimum north slowness
    4,    # Maximum north slowness
    0.01  # Slowness spacing
)
using SeisTau       # See note below
plot(bf, phases=["PKiKP", "PKIKP"])
```


<img src="doc/images/beamforming.png" width="300">

### Vespagrams
Compute a vespagram like so:
```julia
vesp = vespagram(
    t,
    1110, # Start time in s
    1150, # End time
    1,    # Minimum slowness in s/°
    4,    # Maximum slowness
    0.01, # Slowness spacing
)
plot(vesp, phases=["PKiKP", "PKIKP"])
```
<img src="doc/images/vespagram.png" width="600">

## Optional extras

### Plotting
#### Makie.jl
Plots can be created using the [Makie.jl](https://docs.makie.org/stable/)
ecosystem.  To enable this, load a [Makie backend](https://docs.makie.org/stable/explanations/backends/backends) module by doing e.g. `using GLMakie` in
your REPL or module.  You will likely need to add this backend as a package
to your environment, e.g. by doing `import Pkg; Pkg.add("GLMakie").`

The following plotting functions are provided through Makie:
- `plot_array_response`
- `plot_array_response!`
- `plot_beamfomring`
- `plot_beamfomring!`
- `plot_vespagram`
- `plot_vespagram!`

`Makie.plot` and `Makie.plot!` are overloaded such that
`plot(arf_or_bf_or_vesp)` and `plot!(axis, arf_or_bf_or_vesp)` will produce
the right plot, whether the result of `array_response`, `beamform` or
`vespagram` is passed in, as seen above.

#### Plots.jl
As an alternative to Makie plotting, Beamforming.jl supports using
[Plots.jl](https://docs.juliaplots.org/stable/) as a legacy implementation.
If you have installed Plots to your environment, then you can easily
create plots as shown above.

For more help on the Plots.jl plotting commands provided in this package, call up
the online help for `Beamforming.plot`:

```julia
help?> Beamforming.plot
```

### Travel times
To plot the predicted location in slowness and time of arrivals,
you can install [SeisTau.jl](https://github.com/anowacki/SeisTau.jl).
(See the [installation notes](https://github.com/anowacki/SeisTau.jl#installation) for more.)

Simply do `using SeisTau`, then the `phases` option to `plot` will
enable easy predicted phase arrival plotting for `BeamformGrid`s
and `VespaGrid`s, produced respectively by `beamform` and
`vespagram`.
