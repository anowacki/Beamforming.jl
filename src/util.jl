"""
    s_per_degree(s_per_km, r_earth=$R_EARTH_KM_DEFAULT)

Convert horizontal slowness in s/km at the surface of the Earth to s/°.

`r_earth` is the radius of the spherical body in km.
"""
s_per_degree(s_km, r_earth=R_EARTH_KM_DEFAULT) = s_km*2π*r_earth/360

"""
    s_per_km(s_per_degree, r_earth=$R_EARTH_KM_DEFAULT)

Convert horizontal slowness in s/° at the surface of the Earth to s/km.

`r_earth` is the radius of the spherical body in km.
"""
s_per_km(s_degree, r_earth=R_EARTH_KM_DEFAULT) = s_degree*360/(2π*r_earth)
