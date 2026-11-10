# analysis/shape.jl — point-spread and planar-slope shape and direction metrics.
#
# Pure functions of per-cell arrays (x, y in metres about the source, depth h,
# volume V, velocity u, v), plus Run-level conveniences. They replace the
# copies in test/wp1_metrics.jl, test/test_point_spread.jl,
# test/test_octagonal_anisotropy.jl, test/test_planar_slope.jl and
# test/test_planar_symmetry.jl, which are now thin wrappers over these.
#
# Angle conventions (two are in use; each function says which):
#   math     φ = atan(y, x): 0 = east, counter-clockwise (+90° = north).
#            harmonic, sector_ratio, plume_deflection, flow directions.
#   bearing  β = atand(x, y): 0 = north, clockwise (+90° = east).
#            front_radius_by_sector, sample_at_bearings.
#
# Coordinate conventions (`shape_xy`) differ between the scripts; each passes
# its own choices so that its printed numbers are unchanged:
#   wp1_metrics.jl             R = 6_371_008.8 m (IUGG mean radius), origin
#                              snapped to the source cell's centre.
#   test_point_spread.jl,      R = 6_371_000.0 m (R_EARTH = A5Grid._EARTH_R),
#   test_octagonal_anisotropy  origin at the nominal source point.
#   test_planar_symmetry.jl    y = Δlat · 111_320 m/° (M_PER_DEG_LAT, i.e.
#                              R ≈ 6_378_137 m), origin at the nominal latitude.
# The Run-level conveniences default to R = R_EARTH and snap = true (snapping
# removes the up-to-half-a-cell offset between the nominal point and the cell
# FloodA5 actually injects into, which otherwise reads as a deflection while
# the plume is small).

export M_PER_DEG_LAT, shape_xy, harmonic, sector_ratio, front_radius_by_sector,
       radius_cv, cardinal_mask, cardinal_diagonal_means, sample_at_bearings,
       convex_hull, polygon_area_perimeter, polsby_popper, hull_circularity,
       ns_split, plume_deflection, circular_mean_std, flow_directions,
       direction_stats, spread_metrics, planar_metrics, point_spread_shape,
       flow_direction_metrics

"Metres per degree of latitude used by test/test_planar_symmetry.jl (≈ 2π·6_378_137/360)."
const M_PER_DEG_LAT = 111_320.0

# ───────────────────────────── coordinates ─────────────────────────────────

"""
    shape_xy(mesh, lon0, lat0; R = R_EARTH, snap = true, m_per_deg = nothing) -> (x, y)

Cell-centre coordinates (m, x east, y north) about the source.

- `snap = true`: the origin is the centre of the source cell (`source_cell`,
  the cell nearest (lon0, lat0), where FloodA5 injects); `false`: the
  nominal point (lon0, lat0).
- Equirectangular with Earth radius `R` (`local_xy`), unless `m_per_deg` is
  given: then y = Δlat·m_per_deg and x = Δlon·m_per_deg·cos(lat0) (the
  test_planar_symmetry.jl convention, `m_per_deg = M_PER_DEG_LAT`).
"""
function shape_xy(m::MeshGeom, lon0::Real, lat0::Real; R::Real = R_EARTH,
                  snap::Bool = true, m_per_deg = nothing)
    if snap
        i = source_cell(m, lon0, lat0)
        lon0, lat0 = m.lon[i], m.lat[i]
    end
    m_per_deg === nothing && return local_xy(m.lon, m.lat, lon0, lat0; R)
    return (m.lon .- lon0) .* m_per_deg .* cosd(lat0), (m.lat .- lat0) .* m_per_deg
end

# ─────────────────────── point spread: moments, sectors ────────────────────

"""
    harmonic(h, x, y, m; sinus = false) -> Float64

Angular harmonic of a depth field about the origin (math angle φ = atan(y, x),
r = hypot(x, y)):

    A_m = Σ h r² cos(mφ) / Σ h r²      (sinus = true: sin(mφ))

NaN if Σ h r² ≤ 0. On a flat bed with a central source A4 < 0 means
diagonal elongation, A4 > 0 cardinal; A6 is the hexagonal analogue. On the
planar slope (downslope = +x), A2 with `sinus = true` is the signed lateral
skew of the plume (> 0: along the NE–SW diagonal). Zero non-finite h first.
"""
function harmonic(h::AbstractVector, x::AbstractVector, y::AbstractVector, m::Integer;
                  sinus::Bool = false)
    r = hypot.(x, y); φ = atan.(y, x)
    w = h .* r .^ 2
    s = sum(w)
    s > 0 || return NaN
    return sum(w .* (sinus ? sin.(m .* φ) : cos.(m .* φ))) / s
end

"""
    sector_ratio(h, x, y, threshold; min_wet = 20) -> Float64

Cardinal/diagonal ratio of the front radius: cells with h > threshold are
binned into 8 math-angle sectors of 45° centred on E, NE, N, … (sector k
covers φ ∈ [45k − 22.5°, 45k + 22.5°)); a sector's radius is its largest r.
Returns mean(E, N, W, S radii) / mean(NE, NW, SW, SE radii): 1 isotropic,
> 1 the front is further out along the axes. NaN if fewer than `min_wet`
cells are wet or a sector is empty. (test/wp1_metrics.jl `R(h>x)`.)
"""
function sector_ratio(h::AbstractVector, x::AbstractVector, y::AbstractVector,
                      threshold::Real; min_wet::Integer = 20)
    r = hypot.(x, y); φ = atan.(y, x)
    wet = findall(>(threshold), h)
    length(wet) < min_wet && return NaN
    sec = mod.(floor.(Int, (mod.(rad2deg.(φ[wet]), 360.0) .+ 22.5) ./ 45.0), 8)
    rad = [any(==(k), sec) ? maximum(r[wet][sec .== k]) : NaN for k in 0:7]
    return mean(rad[1:2:end]) / mean(rad[2:2:end])
end

"""
    front_radius_by_sector(x, y; n_dirs = 8, centred = true) -> (bearings, radii, counts)

Directional front radius of the points (x, y) (normally the wet cells), by
compass bearing β = atand(x, y) (0 = north, clockwise) in `n_dirs` sectors
of width w = 360°/n_dirs:

- `centred = true`: sector k is centred on bearing (k − 1)·w (N, NE, E, …),
  so cardinal and diagonal are literal sectors (test_octagonal_anisotropy.jl);
- `centred = false`: sector k covers [(k − 1)·w, k·w), centred on (k − ½)·w,
  i.e. each straddles a cardinal/diagonal boundary (test_point_spread.jl).

`radii[k]` is the largest distance √(x² + y²) in sector k (0 if empty),
`counts[k]` the number of points in it; points with r < 1e-6 m are skipped.
`bearings` are the sector centres (°).
"""
function front_radius_by_sector(xs::AbstractVector, ys::AbstractVector;
                                n_dirs::Int = 8, centred::Bool = true)
    radii  = zeros(n_dirs)
    counts = zeros(Int, n_dirs)
    sector_width = 360.0 / n_dirs
    half = sector_width / 2
    for (x, y) in zip(xs, ys)
        r = sqrt(x^2 + y^2)
        r < 1e-6 && continue
        bearing = mod(atand(x, y), 360.0)               # 0° = N, clockwise
        shifted = centred ? mod(bearing + half, 360.0) : bearing
        sector  = clamp(Int(floor(shifted / sector_width)) + 1, 1, n_dirs)
        counts[sector] += 1
        radii[sector] = max(radii[sector], r)
    end
    bearings = centred ? [(k - 1) * sector_width for k in 1:n_dirs] :
                         [(k - 0.5) * sector_width for k in 1:n_dirs]
    return bearings, radii, counts
end

"""
    radius_cv(radii) -> (cv, n_empty)

Coefficient of variation std/mean (sample std, n − 1) of the non-zero sector
radii — 0 for an isotropic front; NaN if fewer than two sectors are
non-empty. `n_empty` = number of zero (empty) sectors, which are excluded.
"""
function radius_cv(radii::AbstractVector)
    nonzero = radii .> 0.0
    cv = count(nonzero) >= 2 ? std(radii[nonzero]) / mean(radii[nonzero]) : NaN
    return cv, count(!, nonzero)
end

"`cardinal_mask(bearings)`: true where a bearing (°) is N, E, S or W (a multiple of 90°)."
cardinal_mask(bearings::AbstractVector) = [b % 90 == 0 for b in bearings]

"""
    cardinal_diagonal_means(bearings, radii)
        -> (cardinal, n_cardinal, diagonal, n_diagonal, ratio)

Mean of the non-zero radii in the cardinal sectors (bearing a multiple of
90°) and in the others, and their ratio (> 1: the front is further out
along the axes, as expected on a 4-connected grid). With 16 sectors
"diagonal" is everything not N/E/S/W. A mean is NaN when no sector qualifies.
"""
function cardinal_diagonal_means(bearings::AbstractVector, radii::AbstractVector)
    card = cardinal_mask(bearings)
    c = radii[card];   c = c[c .> 0]
    d = radii[.!card]; d = d[d .> 0]
    mc = isempty(c) ? NaN : mean(c)
    md = isempty(d) ? NaN : mean(d)
    return (cardinal = mc, n_cardinal = length(c), diagonal = md, n_diagonal = length(d),
            ratio = mc / md)
end

"""
    sample_at_bearings(x, y, radius; n_dirs = 8) -> Vector of (bearing, cell, dist)

For each compass bearing (k − 1)·360°/n_dirs (0 = north, clockwise), the
cell whose centre (x, y) is nearest the point at distance `radius` on that
bearing from the origin, over all cells given; `dist` is how far that
centre is from the target point (should be well under a cell size).
"""
function sample_at_bearings(xs::AbstractVector, ys::AbstractVector, radius::Real;
                            n_dirs::Int = 8)
    out = NamedTuple{(:bearing, :cell, :dist),Tuple{Float64,Int,Float64}}[]
    for b in [(k - 1) * (360.0 / n_dirs) for k in 1:n_dirs]
        tx = radius * sind(b)
        ty = radius * cosd(b)
        best_i = 0
        best_d2 = Inf
        for i in eachindex(xs)
            d2 = (xs[i] - tx)^2 + (ys[i] - ty)^2
            if d2 < best_d2
                best_d2 = d2
                best_i = i
            end
        end
        push!(out, (bearing = b, cell = best_i, dist = sqrt(best_d2)))
    end
    return out
end

# ───────────────────── point spread: hull, circularity ─────────────────────

"""
    convex_hull(x, y) -> Vector{Tuple{Float64,Float64}}
    convex_hull(points::Vector{Tuple{Float64,Float64}})

Convex hull (Andrew's monotone chain) of the points, counter-clockwise,
open (first vertex not repeated); collinear points are dropped. Fewer than
three distinct points are returned as they are.
"""
function convex_hull(pts::Vector{Tuple{Float64,Float64}})
    pts = sort(unique(pts))
    n = length(pts)
    n <= 2 && return pts
    turn(o, a, b) = (a[1]-o[1])*(b[2]-o[2]) - (a[2]-o[2])*(b[1]-o[1])
    lower = Tuple{Float64,Float64}[]
    for p in pts
        while length(lower) >= 2 && turn(lower[end-1], lower[end], p) <= 0
            pop!(lower)
        end
        push!(lower, p)
    end
    upper = Tuple{Float64,Float64}[]
    for p in reverse(pts)
        while length(upper) >= 2 && turn(upper[end-1], upper[end], p) <= 0
            pop!(upper)
        end
        push!(upper, p)
    end
    return vcat(lower[1:end-1], upper[1:end-1])
end
convex_hull(x::AbstractVector, y::AbstractVector) =
    convex_hull(collect(zip(Float64.(x), Float64.(y))))

"""
    polygon_area_perimeter(poly) -> (area, perimeter)

Shoelace area (absolute) and perimeter of a closed polygon given as an open
list of (x, y) vertices; (0, 0) for fewer than three vertices.
"""
function polygon_area_perimeter(hull::AbstractVector{<:Tuple})
    n = length(hull)
    n < 3 && return (0.0, 0.0)
    area = 0.0
    perim = 0.0
    for k in 1:n
        x1, y1 = hull[k]
        x2, y2 = hull[mod1(k+1, n)]
        area  += x1 * y2 - x2 * y1
        perim += sqrt((x2-x1)^2 + (y2-y1)^2)
    end
    return (abs(area) / 2.0, perim)
end

"Polsby–Popper compactness 4πA/P² (1 = circle, π/4 = square); 0 if P ≤ 0."
polsby_popper(area::Real, perimeter::Real) =
    perimeter <= 0.0 ? 0.0 : (4 * π * area) / perimeter^2

"""
    hull_circularity(x, y) -> (circularity, area, perimeter, hull)

Polsby–Popper circularity of the convex hull of the points (normally the
wet-cell centres): the circularity assessment of the FOSS4G 2026 paper.
"""
function hull_circularity(x::AbstractVector, y::AbstractVector)
    hull = convex_hull(x, y)
    area, perim = polygon_area_perimeter(hull)
    return (circularity = polsby_popper(area, perim), area = area, perimeter = perim, hull = hull)
end

# ───────────────────────── planar slope: plume shape ───────────────────────

"""
    ns_split(V, y; band = 50.0) -> NamedTuple

North/south volume split about the source latitude: cells with y > band
are north, y < −band south, the rest the excluded band (cells straddling
the symmetry line). Returns `V_north`, `V_south`, `V_band`,
`V_ns` = V_north + V_south, `n_north`, `n_south`, `n_band` and

    asym = (V_north − V_south) / V_ns     (0 if V_ns ≤ 1e-9)

On a pure east–west slope any asymmetry is numerical. Non-finite V are
summed as they are (zero them first if needed).
"""
function ns_split(V::AbstractVector, y::AbstractVector; band::Real = 50.0)
    north = y .> band
    south = y .< -band
    mid   = .!north .& .!south
    Vn = sum(V[north]); Vs = sum(V[south])
    Vns = Vn + Vs
    return (V_north = Vn, V_south = Vs, V_band = sum(V[mid]), V_ns = Vns,
            n_north = count(north), n_south = count(south), n_band = count(mid),
            asym = Vns > 1e-9 ? (Vn - Vs) / Vns : 0.0)
end

"""
    plume_deflection(V, x, y) -> degrees

Direction of the volume centroid from the origin, math convention:
atand(ΣVy, ΣVx), degrees north of east (the planar test's downslope axis).
Less sensitive than the N/S split to how coarse cells straddle the source
latitude.
"""
plume_deflection(V::AbstractVector, x::AbstractVector, y::AbstractVector) =
    atand(sum(V .* y), sum(V .* x))

# ──────────────────────── planar slope: flow direction ─────────────────────

"""
    circular_mean_std(angles_deg; weights = nothing) -> (mean_deg, std_deg, R)

Circular (Mardia & Jupp) mean direction, standard deviation and mean
resultant length of angles in degrees, optionally weighted (e.g. by speed):
with C = Σw cos θ / Σw, S = Σw sin θ / Σw, R = √(C² + S²) ∈ [0, 1],
mean = atand(S, C), std = √(−2 ln R) in degrees (Inf if R < 1e-9). A small
R means spread-out or bimodal directions (the mean is then not meaningful).
(NaN, NaN, NaN) for no angles or Σw < 1e-12.
"""
function circular_mean_std(angles_deg::AbstractVector;
                           weights::Union{AbstractVector,Nothing} = nothing)
    isempty(angles_deg) && return (NaN, NaN, NaN)
    xs = cosd.(angles_deg)
    ys = sind.(angles_deg)
    if weights === nothing
        mx, my = mean(xs), mean(ys)
    else
        wsum = sum(weights)
        wsum < 1e-12 && return (NaN, NaN, NaN)
        mx = sum(weights .* xs) / wsum
        my = sum(weights .* ys) / wsum
    end
    R = sqrt(mx^2 + my^2)
    mean_deg = atand(my, mx)
    std_deg  = R < 1e-9 ? Inf : rad2deg(sqrt(max(0.0, -2.0 * log(R))))
    return (mean_deg, std_deg, R)
end

"""
    flow_directions(u, v, h; depth_threshold = 1e-3, speed_threshold = 1e-4)
        -> (angles, speeds, n_wet, n_moving)

Velocity direction atand(v, u) (math convention: 0 = east, +90° = north)
and speed √(u² + v²) of every "moving" cell: wet (h > depth_threshold) with
speed > speed_threshold (a near-stationary cell's direction is ill-defined).
`n_wet`, `n_moving`: the cell counts.
"""
function flow_directions(u::AbstractVector, v::AbstractVector, h::AbstractVector;
                         depth_threshold::Real = 1e-3, speed_threshold::Real = 1e-4)
    wet = h .> depth_threshold
    speed = sqrt.(u.^2 .+ v.^2)
    moving = wet .& (speed .> speed_threshold)
    angles = [atand(v[i], u[i]) for i in eachindex(u) if moving[i]]
    speeds = [speed[i]          for i in eachindex(u) if moving[i]]
    return angles, speeds, count(wet), count(moving)
end

"""
    direction_stats(angles, speeds) -> NamedTuple

Flow-direction statistics relative to east (the planar test's downslope):
unweighted circular `mean_dir`, `circ_std`, `R`; speed-weighted
`mean_dir_w`, `circ_std_w`, `R_w` (cells with real momentum dominate; the
thin front contributes mostly noise); `mean_abs_dev` = mean |angle|
(≈ 90° for random directions); `frac_within_30` = share of angles within
±30° of east (a low share with a near-zero mean flags a bimodal split).
"""
function direction_stats(angles::AbstractVector, speeds::AbstractVector)
    mean_dir, circ_std, R = circular_mean_std(angles)
    mean_dir_w, circ_std_w, R_w = circular_mean_std(angles; weights = speeds)
    return (mean_dir = mean_dir, circ_std = circ_std, R = R,
            mean_dir_w = mean_dir_w, circ_std_w = circ_std_w, R_w = R_w,
            mean_abs_dev = mean(abs.(angles)),
            frac_within_30 = count(a -> abs(a) <= 30.0, angles) / length(angles))
end

# ───────────────────────────── Run level ───────────────────────────────────

_finite0(a::AbstractVector) = (a = copy(a); a[.!isfinite.(a)] .= 0.0; a)

"""
    spread_metrics(r, k; lon0, lat0, R = R_EARTH, snap = true,
                   thresholds = (1e-4, 1e-2, 5e-2), wet_threshold = 1e-4) -> NamedTuple

Point-spread shape at frame k about the source (lon0, lat0) (`shape_xy`):
`t`, `n_wet` (h > wet_threshold), `A4`, `A6` (`harmonic`) and `ratio`, the
`sector_ratio` at each depth threshold. Non-finite depths count as dry.
test/wp1_metrics.jl `spread` uses R = 6_371_008.8, snap = true.
"""
function spread_metrics(r::Run, k::Integer; lon0::Real, lat0::Real, R::Real = R_EARTH,
                        snap::Bool = true, thresholds = (1e-4, 1e-2, 5e-2),
                        wet_threshold::Real = 1e-4)
    h = _finite0(depth(r, k))
    x, y = shape_xy(r.mesh, lon0, lat0; R, snap)
    return (t = r.t[k], n_wet = count(>(wet_threshold), h),
            A4 = harmonic(h, x, y, 4), A6 = harmonic(h, x, y, 6),
            ratio = map(th -> sector_ratio(h, x, y, th), Tuple(thresholds)))
end

"""
    planar_metrics(r, k; lon0, lat0, R = R_EARTH, snap = true, band = 50.0,
                   m_per_deg = nothing) -> NamedTuple

Planar-slope plume shape at frame k about the source (`shape_xy`): `t`,
the `ns_split` fields (`asym`, `V_north`, `V_south`, …),
`A2` = `harmonic(h, x, y, 2; sinus = true)` and `deflection`
(`plume_deflection`, °). Non-finite volumes and depths count as zero.
test/wp1_metrics.jl `planar` uses R = 6_371_008.8, snap = true;
test_planar_symmetry.jl uses `ns_split` with snap = false,
m_per_deg = M_PER_DEG_LAT and non-finite volumes left as they are.
"""
function planar_metrics(r::Run, k::Integer; lon0::Real, lat0::Real, R::Real = R_EARTH,
                        snap::Bool = true, band::Real = 50.0, m_per_deg = nothing)
    x, y = shape_xy(r.mesh, lon0, lat0; R, snap, m_per_deg)
    V = _finite0(volume(r, k))
    h = _finite0(depth(r, k))
    return (; t = r.t[k], ns_split(V, y; band)...,
            A2 = harmonic(h, x, y, 2; sinus = true), deflection = plume_deflection(V, x, y))
end

"""
    point_spread_shape(r, k; lon0, lat0, R = R_EARTH, snap = true,
                       saturation_threshold = 1.0, n_dirs = 8, centred = true) -> NamedTuple

Shape of the wetted region (cells with saturation ≥ saturation_threshold)
at frame k about the source (`shape_xy`): `t`, `n_wet`, the hull
`circularity` (Polsby–Popper), `area`, `perimeter`; the directional front
radius `bearings`, `radii`, `counts` (`front_radius_by_sector`); its
`cv_radius` and `n_empty_sectors` (`radius_cv`).
test_point_spread.jl uses R = R_EARTH, snap = false, centred = false;
test_octagonal_anisotropy.jl R = R_EARTH, snap = false, centred = true.
"""
function point_spread_shape(r::Run, k::Integer; lon0::Real, lat0::Real, R::Real = R_EARTH,
                            snap::Bool = true, saturation_threshold::Real = 1.0,
                            n_dirs::Int = 8, centred::Bool = true)
    wet = field(r, "saturation", k) .>= saturation_threshold
    x, y = shape_xy(r.mesh, lon0, lat0; R, snap)
    xs, ys = x[wet], y[wet]
    hc = hull_circularity(xs, ys)
    bearings, radii, counts = front_radius_by_sector(xs, ys; n_dirs, centred)
    cv, n_empty = radius_cv(radii)
    return (t = r.t[k], n_wet = count(wet), circularity = hc.circularity, area = hc.area,
            perimeter = hc.perimeter, cv_radius = cv, n_empty_sectors = n_empty,
            bearings = bearings, radii = radii, counts = counts)
end

"""
    flow_direction_metrics(r, k; depth_threshold = 1e-3, speed_threshold = 1e-4) -> NamedTuple

Velocity-direction statistics at frame k (needs `vel_u`, `vel_v`): `t`,
`n_wet`, `n_moving` (`flow_directions`), the `direction_stats` fields, and
the `angles` and `speeds` of the moving cells.
"""
function flow_direction_metrics(r::Run, k::Integer; depth_threshold::Real = 1e-3,
                                speed_threshold::Real = 1e-4)
    has_field(r, "vel_u") || error("$(r.label): no vel_u/vel_v in the output")
    angles, speeds, n_wet, n_moving =
        flow_directions(field(r, "vel_u", k), field(r, "vel_v", k), depth(r, k);
                        depth_threshold, speed_threshold)
    return (; t = r.t[k], n_wet, n_moving, direction_stats(angles, speeds)...,
            angles, speeds)
end
