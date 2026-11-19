# cases/registry/_planar_common.jl — helpers shared by the planar cases
# (planar_slope.jl, planar_embankment.jl, planar_wall.jl), which generate their terrain
# with test/planar_embankment/generate_planar_embankment_dem.py. (Registry
# files are included in name order; this one comes first.)

"""
    _planar_point(; domain_km = (4.0, 2.0), centre = (-0.017, 51.0), frac = (0.1, 0.5))
        -> (lat, lon)

The injection point `generate_planar_embankment_dem.py` writes to its params
file (`injection_lat`, `injection_lon`), computed as the generator does it
(equirectangular, 111 320 m per degree): x = frac[1]·W from the west edge,
y = frac[2]·H from the south edge. The job arguments need it before the
generator has run; the cases check it against the params file.
"""
function _planar_point(; domain_km = (4.0, 2.0), centre = (-0.017, 51.0), frac = (0.1, 0.5))
    m_lon = 111320.0 * cos(deg2rad(centre[2]))
    lon0 = centre[1] - (domain_km[1] * 1000.0 / 2.0) / m_lon
    lat0 = centre[2] - (domain_km[2] * 1000.0 / 2.0) / 111320.0
    return (lat0 + frac[2] * domain_km[2] * 1000.0 / 111320.0,
            lon0 + frac[1] * domain_km[1] * 1000.0 / m_lon)
end

"The generator's parameters (`<dir>/planar_embankment_params.json`) as a Dict."
_planar_params(dir) = FA.FO.JSON3.read(read(joinpath(dir, "planar_embankment_params.json"), String),
                                       Dict{String,Any})

"Cells' y extent (m from the DEM's south edge, generator convention) from their polygons."
function _planar_yrange(mesh, p)
    lat0 = Float64(p["lat0"])
    ymin = [minimum(view(P, 2, :)) for P in mesh.polygons]
    ymax = [maximum(view(P, 2, :)) for P in mesh.polygons]
    return (ymin .- lat0) .* 111320.0, (ymax .- lat0) .* 111320.0
end
