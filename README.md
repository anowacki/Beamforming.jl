# Beamforming

Beamforming is a Julia module for the array analysis of seismic
data, using [Seis.jl](https://github.com/anowacki/Seis.jl).

## Build status
[![Build Status](https://travis-ci.org/anowacki/Beamforming.jl.svg?branch=master)](https://travis-ci.org/anowacki/Beamforming.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/voywfa1yu2f72fuk/branch/master?svg=true)](https://ci.appveyor.com/project/AndyNowacki/beamforming-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/anowacki/Beamforming.jl/badge.svg?branch=master)](https://coveralls.io/github/anowacki/Beamforming.jl?branch=master)


## Installation

```julia
julia> ] # press ']' to enter pkg mode

(v1.3) pkg> add https://github.com/anowacki/Geodesics.jl https://github.com/anowacki/Seis.jl https://github.com/anowacki/Beamforming.jl
```

## Using

The main functions exported are:

- `array_response`: Compute the array response function for a set
  of stations
- `beamform`: Compute the beam power across a grid of slowness points
  (as in f–k analysis)
- `vespagram`: Compute a slowness vespagram

If you install [Plots.jl](https://github.com/JuliaPlots/Plots.jl),
then you can visualise the output of each of these functions
easily using the `plot` function.

Docstrings are exhaustive, and can be consulted for usage of these
functions.

## Optional extras

To plot the predicted location in slowness and time of arrivals,
you can install [SeisTau.jl](https://github.com/anowacki/SeisTau.jl).
(See the [installation notes](https://github.com/anowacki/SeisTau.jl#installation) for more.)

Simply do `using SeisTau`, then the `phases` option to `plot` will
enable easy predicted phase arrival plotting for `BeamformGrid`s
and `VespaGrid`s, produced respectively by `beamform` and
`vespagram`.
