# cases/scripts/sgs_strip_tables.jl — an SGS mesh's tables against the
# strip's analytic terrain, without a simulation (case `sgs_strip`,
# docs/TEST_CASES.md §2c).
#
#   FLOODA5_NO_VIS=1 julia --project=. cases/scripts/sgs_strip_tables.jl MESH.parquet PARAMS.json OUT.json
#
# Uses the model's own table lookups (A5Grid.sgs_table, face_apk and
# the volume curve as wse_from_volume inverts it) and compares, at water
# levels δ above the local floodplain (z_f at the cell centre, or at the face
# midpoint for an edge):
#   cells  stored volume V(η) and wetted area against the analytic surface
#          (cases/scripts/sgs_strip_geometry.jl), integrated over the cell
#          polygon on a 0.1 m grid
#   edges  the face table (flow area A, wetted width P, conveyance
#          K = Σ h^(5/3); A5Grid.face_apk, one table per face since 0.7.0)
#          against the face's analytic A, P and K, sampled every 1 cm
# Edges are grouped as "channel" (the face crosses the channel) and
# "floodplain". For a Cartesian strip (rows along x) also the uniform-flow
# discharge its tables give: Q(δ) = Σ over the faces of one cross-section of
# K √S / n (2-D form; `Q_lumped`: A (A/P)^(2/3) √S / n) at the upstream
# cell's level, so the normal level the tables imply.
# Writes OUT.json: per level the error statistics, curves for one cell in the
# middle of the strip and its downstream and lateral faces, and Q(δ).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
include(joinpath(@__DIR__, "sgs_strip_geometry.jl"))
using Logging, Printf, Statistics

const MESHF, PARAMSF, OUTF = ARGS[1], ARGS[2], ARGS[3]
const NMAN = 0.03
const LEVELS = [-0.75, -0.5, -0.25, -0.1, 0.1, 0.25, 0.4, 0.6, 1.0]   # δ (m) above the floodplain
const pd = JSON3.read(read(PARAMSF, String), Dict{String,Any})
const p = _ss_params(pd)
const LON0, LAT0, RE = Float64(pd["lon0"]), Float64(pd["lat0"]), Float64(pd["earth_radius_m"])
xy(lon, lat) = (deg2rad(lon - LON0) * RE * cosd(LAT0), deg2rad(lat - LAT0) * RE)

const mesh = with_logger(() -> load_mesh_geoparquet(MESHF), SimpleLogger(stderr, Logging.Warn))
const n = length(mesh.cells)
const ids = [c.id for c in mesh.cells]
const idx = Dict(id => i for (i, id) in enumerate(ids))
const tabs = [sgs_table(mesh, i) for i in 1:n]
const poly = [[xy(v[1], v[2]) for v in c.boundary] for c in mesh.cells]
const ctr = [xy(c.center_lon, c.center_lat) for c in mesh.cells]

"Stored volume at WSE η, as wse_from_volume inverts it (linear between knots; full cell area above z_max)."
function vol_at(t, η)
    eb, vc = t.elev_bins, t.vol_curve
    η <= t.z_min && return 0.0
    η >= eb[end] && return vc[end] + (η - eb[end]) * t.cell_area
    k = searchsortedlast(eb, η)
    return vc[k] + (η - eb[k]) / (eb[k+1] - eb[k] + 1e-15) * (vc[k+1] - vc[k])
end

function inpoly(x, y, P)
    c = false; m = length(P)
    for i in 1:m
        j = i == 1 ? m : i - 1
        (xi, yi), (xj, yj) = P[i], P[j]
        ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) && (c = !c)
    end
    return c
end

"Analytic surface points of a cell (0.1 m grid over its polygon): (z values, area per point)."
function cell_points(i; h = 0.1)
    P = poly[i]
    x0, x1 = extrema(first.(P)); y0, y1 = extrema(last.(P))
    z = Float64[]
    for x in (x0 + h / 2):h:x1, y in (y0 + h / 2):h:y1
        inpoly(x, y, P) && push!(z, _ss_surface(p, x, y))
    end
    return z, h^2
end

"Analytic face profile (1 cm spacing): (z values, spacing, midpoint x, y)."
function face_profile(i, j)
    e = A5Grid._shared_edge(mesh.cells[i].boundary, mesh.cells[j].boundary)
    e === nothing && return nothing
    (xa, ya), (xb, yb) = xy(e[1], e[2]), xy(e[3], e[4])
    len = hypot(xb - xa, yb - ya); m = max(2, round(Int, len / 0.01))
    s = [(k - 0.5) / m for k in 1:m]
    z = [_ss_surface(p, xa + f * (xb - xa), ya + f * (yb - ya)) for f in s]
    return z, len / m, (xa + xb) / 2, (ya + yb) / 2, len
end

# ── cells: V(η) and wetted area ───────────────────────────────────────────
reach(x) = 0.1p.L <= x <= 0.9p.L
cell_err = Dict(δ => Float64[] for δ in LEVELS); area_err = Dict(δ => Float64[] for δ in LEVELS)
cell_pts = Dict{Int,Any}()
for i in 1:n
    reach(ctr[i][1]) || continue
    z, dA = cell_points(i); cell_pts[i] = (z, dA)
    zf = _ss_zf(p, ctr[i][1]); At = tabs[i].cell_area
    for δ in LEVELS
        η = zf + δ
        Va = sum(max(η - zz, 0.0) for zz in z) * dA
        Aw = count(<(η), z) * dA
        push!(cell_err[δ], (vol_at(tabs[i], η) - Va) / max(Va, 0.05 * At))     # floor: 5 cm over the cell
        push!(area_err[δ], (A5Grid.wetted_area_from_wse(tabs[i], η) - Aw) / At)
    end
end

# ── edges: A_edge, R_edge ─────────────────────────────────────────────────
slot(i, j) = findfirst(==(ids[j]), mesh.adjacency[ids[i]])
edges = Tuple{Int,Int}[]
for i in 1:n, id in mesh.adjacency[ids[i]]
    j = get(idx, id, 0); j > i && push!(edges, (i, j))
end
edge_ratio = Dict((c, δ) => Float64[] for c in ("channel", "floodplain"), δ in LEVELS)
edge_rratio = Dict((c, δ) => Float64[] for c in ("channel", "floodplain"), δ in LEVELS)
edge_kratio = Dict((c, δ) => Float64[] for c in ("channel", "floodplain"), δ in LEVELS)
"The face table of the face i|j as the kernel reads it (the lower-index cell's slot): (A, P, K) at η."
face_tab(i, j, η) = i < j ? A5Grid.face_apk(tabs[i], slot(i, j), η) : A5Grid.face_apk(tabs[j], slot(j, i), η)
for (i, j) in edges
    (reach(ctr[i][1]) && reach(ctr[j][1])) || continue
    fp = face_profile(i, j); fp === nothing && continue
    z, ds, xm, _, _ = fp
    cls = minimum(z) < _ss_zf(p, xm) - p.D / 2 ? "channel" : "floodplain"
    si, sj = slot(i, j), slot(j, i)
    for δ in LEVELS
        η = _ss_zf(p, xm) + δ
        Aa = sum(max(η - zz, 0.0) for zz in z) * ds; Pa = count(<(η), z) * ds
        Ka = sum(max(η - zz, 0.0)^(5 / 3) for zz in z) * ds
        Aa < 1e-3 && continue
        At, Pt, Kt = face_tab(i, j, η)
        push!(edge_ratio[(cls, δ)], At / Aa)
        push!(edge_rratio[(cls, δ)], (Pt > 0 ? At / Pt : 0.0) / (Aa / Pa))
        push!(edge_kratio[(cls, δ)], Kt / Ka)
    end
end

# ── one cell in the middle of the strip: curves ───────────────────────────
ic = argmin([hypot(c[1] - p.L / 2, c[2] - p.YC) for c in ctr])
t = tabs[ic]; zf_c = _ss_zf(p, ctr[ic][1])
zc, dAc = haskey(cell_pts, ic) ? cell_pts[ic] : cell_points(ic)
grid = collect(range(t.z_min - 0.05, zf_c + 1.5; length = 155))
curve(f) = [f(η) for η in grid]
"Downstream (largest x of the face midpoint) and one lateral neighbour of the cell."
function pick_faces(i)
    best = nothing; lat = nothing
    for id in mesh.adjacency[ids[i]]
        j = get(idx, id, 0); j == 0 && continue
        fp = face_profile(i, j); fp === nothing && continue
        dx, dy = fp[3] - ctr[i][1], fp[4] - ctr[i][2]
        (best === nothing || dx > best[2]) && (best = (j, dx))
        (abs(dy) > abs(dx) && (lat === nothing || abs(dy) > lat[2])) && (lat = (j, abs(dy)))
    end
    return best === nothing ? 0 : best[1], lat === nothing ? 0 : lat[1]
end
jd, jl = pick_faces(ic)
function face_curves(i, j)
    j == 0 && return nothing
    z, ds, xm, ym, len = face_profile(i, j)
    At = curve(η -> face_tab(i, j, η)[1]); Pt = curve(η -> face_tab(i, j, η)[2]); Kt = curve(η -> face_tab(i, j, η)[3])
    Aa = curve(η -> sum(max(η - zz, 0.0) for zz in z) * ds)
    Pa = curve(η -> count(<(η), z) * ds)
    Ka = curve(η -> sum(max(η - zz, 0.0)^(5 / 3) for zz in z) * ds)
    return Dict("neighbour" => j, "x_mid" => xm, "y_mid" => ym, "length_m" => len,
                "sill" => minimum(z), "sill_table" => tabs[i].face_z0[slot(i, j)],
                "top_table" => tabs[i].face_z0[slot(i, j)] + tabs[i].face_zr[slot(i, j)],
                "A_edge" => At, "A_analytic" => Aa, "P_edge" => Pt, "P_analytic" => Pa,
                "K_edge" => Kt, "K_analytic" => Ka,
                "z_max_i" => tabs[i].z_max, "z_max_j" => tabs[j].z_max, "z_min_j" => tabs[j].z_min)
end
cell_curves = Dict("index" => ic, "x" => ctr[ic][1], "y" => ctr[ic][2], "z_floodplain" => zf_c,
    "z_min" => t.z_min, "z_max" => t.z_max, "cell_area" => t.cell_area, "elevation" => mesh.static_vars["elevation"][ic],
    "wse" => grid, "V_table" => curve(η -> vol_at(t, η)),
    "V_analytic" => curve(η -> sum(max(η - zz, 0.0) for zz in zc) * dAc),
    "Awet_table" => curve(η -> A5Grid.wetted_area_from_wse(t, η)),
    "Awet_analytic" => curve(η -> count(<(η), zc) * dAc),
    "elev_bins" => t.elev_bins, "vol_curve" => t.vol_curve, "area_curve" => t.area_curve,
    "face_downstream" => face_curves(ic, jd), "face_lateral" => face_curves(ic, jl))

# ── Cartesian strip: uniform-flow discharge the tables give ───────────────
# Faces between cells of one row (|Δy| < 1 m) whose midpoint is at x_f, the
# section; each evaluated at the upstream cell's level: η_i = z_f(x_i) + δ.
qcurve = nothing
rows_along = [(i, j) for (i, j) in edges if abs(ctr[i][2] - ctr[j][2]) < 1.0]
if haskey(mesh.static_vars, "cell_size_m") && !isempty(rows_along)
    xs = [(ctr[i][1] + ctr[j][1]) / 2 for (i, j) in rows_along]
    xf = xs[argmin(abs.(xs .- p.L / 2))]
    sec = [(i, j) for (k, (i, j)) in enumerate(rows_along) if abs(xs[k] - xf) < 0.5]
    δs = collect(range(-p.D + 0.02, 3.0; length = 300))
    # 2-D form (the default): Q = K·√S/n; lumped: A·(A/P)^(2/3)·√S/n
    Qof(lumped) = map(δs) do δ
        sum(sec) do (i, j)
            up = ctr[i][1] < ctr[j][1] ? i : j
            A, P, K = face_tab(i, j, _ss_zf(p, ctr[up][1]) + δ)
            (lumped ? (P > 0 ? A * (A / P)^(2 / 3) : 0.0) : K) * sqrt(p.S) / NMAN
        end
    end
    Qd = Qof(false)
    global qcurve = Dict("x_face" => xf, "n_faces" => length(sec), "delta" => δs, "Q" => Qd, "Q_lumped" => Qof(true),
                         "Q_max" => maximum(Qd), "delta_at_Q_max" => δs[argmax(Qd)])
end

stat(v) = isempty(v) ? Dict("n" => 0) :
    Dict("n" => length(v), "median" => median(v), "min" => minimum(v), "max" => maximum(v),
         "max_abs" => maximum(abs, v))
out = Dict("mesh" => MESHF, "n_cells" => n, "levels" => LEVELS,
    "volume_rel_err" => Dict(string(δ) => stat(cell_err[δ]) for δ in LEVELS),
    "wet_area_err" => Dict(string(δ) => stat(area_err[δ]) for δ in LEVELS),
    "edge_area_ratio" => Dict(c => Dict(string(δ) => stat(edge_ratio[(c, δ)]) for δ in LEVELS) for c in ("channel", "floodplain")),
    "edge_radius_ratio" => Dict(c => Dict(string(δ) => stat(edge_rratio[(c, δ)]) for δ in LEVELS) for c in ("channel", "floodplain")),
    "edge_conveyance_ratio" => Dict(c => Dict(string(δ) => stat(edge_kratio[(c, δ)]) for δ in LEVELS) for c in ("channel", "floodplain")),
    "cell" => cell_curves, "uniform_flow" => qcurve)
open(OUTF, "w") do f; JSON3.write(f, out); end
@printf("%s: %d cells; V error at δ = +0.4: max %.3f; channel edge A ratio at δ = +0.4: median %.3f; floodplain %.3f\n",
        OUTF, n, stat(cell_err[0.4])["max_abs"], get(stat(edge_ratio[("channel", 0.4)]), "median", NaN),
        get(stat(edge_ratio[("floodplain", 0.4)]), "median", NaN))
qcurve === nothing || @printf("  uniform flow at the bed slope: at most %.2f m³/s (δ = %.2f m)\n", qcurve["Q_max"], qcurve["delta_at_Q_max"])
