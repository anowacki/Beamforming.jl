# Spherical and cartesian coordinate conversions

"""
    azimuth(lon1, lat1, lon2, lat2, degrees::Bool=true) -> az

Compute the azimuth `az` from point (`lon1`, `lat1`) to (`lon2`, `lat2`) on the sphere.
Points and azimuth are read and returned in degrees by default; use `degrees=false`
for radians.
"""
function azimuth(lon1, lat1, lon2, lat2, degrees::Bool=true)
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    azimuth = atan(sin(lon2-lon1)*cos(lat2),
                   cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon2-lon1))
    degrees ? rad2deg(azimuth) : azimuth
end

"""
    cart2geog(x, y, z, degrees::Bool=true) -> lon, lat, r

Compute the longitude, latitude and radius given the cartesian coordinates `x`,
`y` and `z`, where `x` is at (lon,lat) = (0,0), `y` is at (90°,0) and `z` is
through lat = 90°.
"""
function cart2geog(x, y, z, degrees::Bool=true)
    r = sqrt(x^2 + y^2 + z^2)
    r == 0. && return zero(x), zero(x), zero(x)
    lon = atan(y, x)
    lat = asin(z/r)
    degrees ? (rad2deg(lon), rad2deg(lat), r) : (lon, lat, r)
end

function cart2geog(x::AbstractArray, y::AbstractArray, z::AbstractArray, degrees::Bool=true)
    dims = size(x)
    dims == size(y) == size(z) || throw(ArgumentError("All arrays must have same length"))
    T = promote_type(float.(eltype.((x, y, z)))...)
    lon, lat, r = Array{T}(dims), Array{T}(dims), Array{T}(dims)
    for i in eachindex(lon)
        lon[i], lat[i], r[i] = cart2geog(x[i], y[i], z[i], degrees)
    end
    lon, lat, r
end

"""
    delta(lon1, lat1, lon2, lat2, degrees::Bool=true) -> d

Compute the angular distance `d` on the sphere between two points, (`lon1`, `lat1`)
and (`lon2`, `lat2`).  Points and distance are read and returned in degrees
by default; use `degrees=false` for radians.
"""
function delta(lon1, lat1, lon2, lat2, degrees::Bool=true)
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    d = atan(sqrt(
               (cos(lat2)*sin(lon2-lon1))^2 + (cos(lat1)*sin(lat2) -
                sin(lat1)*cos(lat2)*cos(lon2-lon1))^2),
               sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1)
              )
    degrees ? rad2deg(d) : d
end

"""
    geog2cart(lon, lat, r, degrees::Bool=true) -> x, y, z
    geog2cart(lon, lat, degrees::Bool=true) -> x, y, z

Return the cartesian coordinates given the geographic longitude, latitude and
radius `lon`, `lat` and `r`.

If `r` is not given, points are returned on the unit sphere.
"""
function geog2cart(lon, lat, r, degrees::Bool=true)
    points_valid(lon, lat, degrees) || error("geog2cart: Points are not on the sphere")
    degrees && begin lon, lat = deg2rad(lon), deg2rad(lat) end
    x = r*cos(lon)*cos(lat)
    y = r*sin(lon)*cos(lat)
    z = r*sin(lat)
    x, y, z
end
function geog2cart(lon::AbstractArray, lat::AbstractArray, r::AbstractArray, degrees::Bool=true)
    size(lon) == size(lat) == size(r) ||
        throw(ArgumentError("Sizes of lon, lat and r must be the same"))
    T = promote_type(float.(eltype.((lon, lat, r)))...)
    x = Array{T}(undef, size(lon))
    y, z = similar(x), similar(x)
    for i in eachindex(lon)
        x[i], y[i], z[i] = geog2cart(lon[i], lat[i], r[i], degrees)
    end
    x, y, z
end
geog2cart(lon::AbstractArray, lat::AbstractArray, degrees::Bool=true) =
    geog2cart(lon, lat, fill(one(eltype(lon)), size(lon)), degrees)

"""
    points_valid(lon, lat, degrees::Bool=true) -> ::Bool

Return `true` if all points in arrays `lon` and `lat` are on the sphere.
`lon` is ignored, but `lat` is checked to see if points are in the range
-90°–90° (-π–π).  Points are read and returned in degrees by default; use
`degrees=false` for radians.
"""
points_valid(lon, lat, degrees::Bool=true) =
    degrees ? !any(abs.(lat) .> 90.) : !any(abs.(lat) .> pi/2.)
