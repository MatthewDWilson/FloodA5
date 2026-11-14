# cases/scripts/sgs_strip_geometry.jl — the SGS strip's terrain and its
# analytic hydraulics (case `sgs_strip`; docs/TEST_CASES.md §2c). Pure Julia,
# no packages: included by cases/scripts/sgs_strip_dem.jl (the DEM),
# cases/scripts/sgs_strip_tables.jl (SGS tables against the terrain) and
# cases/registry/sgs_strip.jl (the case's analytic answers), so all three use
# one definition. Names carry the registry's `_ss_` prefix (registry files
# share the FloodA5Cases namespace).
#
# Parameters `p` (a NamedTuple; `_ss_params(dict)` from sgs_strip_params.json):
# L, B (flume length, width), W, D (channel width, depth), S (bed slope),
# P (DEM node spacing), H (wall height), YC (channel centre line y), SX
# (floodplain cross slope), Z0 (floodplain at x = 0), XB, S2 (below x = XB the
# bed slope is S + S2: a steep, supercritical outlet reach).

"Parameters as a NamedTuple from a Dict with sgs_strip_params.json's keys."
_ss_params(d) = (L = Float64(d["length_m"]), B = Float64(d["width_m"]), W = Float64(d["channel_width_m"]),
                 D = Float64(d["channel_depth_m"]), S = Float64(d["slope"]), P = Float64(d["pixel_m"]),
                 H = Float64(d["wall_m"]), YC = Float64(d["channel_offset_m"]), SX = Float64(d["cross_slope"]),
                 Z0 = Float64(d["z0_m"]), XB = Float64(get(d, "break_x_m", 1e30)),
                 S2 = Float64(get(d, "break_slope", 0.0)))

"Floodplain elevation at x (m) (no cross slope)."
_ss_zf(p, x) = p.Z0 - p.S * x - p.S2 * max(x - p.XB, 0.0)

"The bed at a DEM node (x, y): the generator's piecewise definition (walls, channel, floodplain)."
function _ss_bed(p, x, y)
    zf = _ss_zf(p, x)
    (x < -p.P || abs(y) > p.B / 2 + p.P) && return zf + p.H
    a = abs(y - p.YC)
    a < p.W / 2 && return zf - p.D
    return zf + p.SX * (a - p.W / 2)
end

"""
    _ss_surface(p, x, y)

The terrain the model samples: bilinear interpolation of `_ss_bed` between the
DEM nodes at (k + ½)·P — what `A5Grid._bilinear` returns from the DEM, up to
Float32 rounding. Each vertical wall becomes a ramp of width P centred on it.
"""
function _ss_surface(p, x, y)
    fx = x / p.P - 0.5; fy = y / p.P - 0.5
    i = floor(fx); j = floor(fy); a = fx - i; b = fy - j
    x0 = (i + 0.5) * p.P; x1 = x0 + p.P; y0 = (j + 0.5) * p.P; y1 = y0 + p.P
    return (1 - a) * (1 - b) * _ss_bed(p, x0, y0) + a * (1 - b) * _ss_bed(p, x1, y0) +
           (1 - a) * b * _ss_bed(p, x0, y1) + a * b * _ss_bed(p, x1, y1)
end

"Bed relative to the floodplain across a section, z(x, y) − z_f(x) (independent of x for x > P)."
_ss_zrel(p, y) = _ss_surface(p, 1000.0 * p.P, y) - _ss_zf(p, 1000.0 * p.P)

"""
    _ss_section(p; ylim = p.B / 2, dy = 0.005) -> (y, zrel)

Cross-section sample points (midpoint rule) for |y| ≤ ylim: `ylim = B/2` is a
Cartesian mesh of width B (its edges are the closed boundary); `B/2 + 2P`
includes the side-wall ramps (the flume an A5 mesh sees).
"""
function _ss_section(p; ylim = p.B / 2, dy = 0.005)
    y = collect((-ylim + dy / 2):dy:ylim)
    return y, [_ss_zrel(p, v) for v in y], dy
end

"""
    _ss_conveyance(p, δ, mode; ylim, rows) -> K (m³/s at unit slope, n = 1)

Conveyance of the section at water level δ above the floodplain:
- `:twod`: Σ over lateral strips of h^(5/3) dy — 2-D depth-averaged flow with
  no lateral exchange (what a fine 2-D model converges to at uniform flow)
- `:cells`: rows of width `rows` (m) from y = −ylim, each with its mean bed:
  Σ D·(δ − z̄)^(5/3) — a standard (cell-mean) grid of those rows
- `:lumped`: rows with their true section: Σ A^(5/3)/P_w^(2/3), A the flow area
  and P_w the wetted width — the SGS R–A edge formula with exact edge tables
"""
function _ss_conveyance(p, δ, mode::Symbol; ylim = p.B / 2, rows = 0.0)
    y, z, dy = _ss_section(p; ylim)
    mode === :twod && return sum(max(δ - zi, 0.0)^(5 / 3) for zi in z) * dy
    nr = round(Int, 2ylim / rows)
    K = 0.0
    for r in 1:nr
        sel = [k for k in eachindex(y) if -ylim + (r - 1) * rows <= y[k] < -ylim + r * rows]
        if mode === :cells
            zm = sum(z[sel]) / length(sel)
            K += rows * max(δ - zm, 0.0)^(5 / 3)
        else
            A = sum(max(δ - z[k], 0.0) for k in sel) * dy
            Pw = count(k -> δ > z[k], sel) * dy
            Pw > 0 && (K += A^(5 / 3) / Pw^(2 / 3))
        end
    end
    return K
end

"Storage per metre (m²) at level δ: true water volume (`:water`) or above the row-mean beds (`:cells`)."
function _ss_storage(p, δ, mode::Symbol = :water; ylim = p.B / 2, rows = 0.0)
    y, z, dy = _ss_section(p; ylim)
    mode === :water && return sum(max(δ - zi, 0.0) for zi in z) * dy
    nr = round(Int, 2ylim / rows)
    return sum(rows * max(δ - sum(z[k] for k in eachindex(y) if -ylim + (r - 1) * rows <= y[k] < -ylim + r * rows) /
                          count(k -> -ylim + (r - 1) * rows <= y[k] < -ylim + r * rows, eachindex(y)), 0.0)
               for r in 1:nr)
end

"""
    _ss_normal_level(p, Q, n, mode; kw...) -> δ

Uniform-flow water level above the floodplain (m) at discharge Q, Manning n,
bed slope p.S: Q = K(δ)·√S/n solved by bisection on δ ∈ (−D, 20].
"""
function _ss_normal_level(p, Q, n, mode::Symbol; kw...)
    lo, hi = -p.D + 1e-9, 20.0
    for _ in 1:60
        m = (lo + hi) / 2
        _ss_conveyance(p, m, mode; kw...) * sqrt(p.S) / n < Q ? (lo = m) : (hi = m)
    end
    return (lo + hi) / 2
end
