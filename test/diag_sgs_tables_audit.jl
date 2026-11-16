# test/diag_sgs_tables_audit.jl — SGS Phase 0 audit of a mesh's cell tables
# (FloodA5_SGS_Conveyance_Fix_Plan_2026-11.md, Phase 0 items 2–4). Measures;
# changes nothing. Results: docs/notes/SGS_Phase0_Audit_2026-11.md.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/diag_sgs_tables_audit.jl MESH.parquet DEM.tif OUT.json
#       [--run SGS_RUN.h5] [--compare-run STD_RUN.h5]     wet cells and Courant numbers at the last frame
#       [--strip-params sgs_strip_params.json]            also against the strip's analytic terrain
#       [--embankment-params planar_embankment_params.json]   list the cells on the embankment
#       [--subset 500] [--halton-hi 4096] [--seed 1]      cells rebuilt with more Halton candidates
#
# A. Hypsometry. Reference: the DEM's bilinear surface (A5Grid._bilinear, the
#    surface the model samples) on a regular grid over the cell polygon, at
#    h = min(pixel/8, √A/64); each point stands for A/n (pixel/2 is not
#    enough where the DEM has ramps one pixel wide: the strip's channel walls). Compared with the
#    stored table as the model inverts it (FA.sgs_table_volume / _wse), with a
#    100-level table built from the dense samples (level interpolation alone)
#    and, for a subset, with tables rebuilt from 512 and 4096 Halton candidates.
# B. Sub-basins of each cell on the same grid (FA.basin_merge_tree).
# C. Courant numbers at the run's mean step for three depth measures.
#
# The table loader, DEM reader, Halton points and point-in-polygon test are
# A5Grid's own; everything else is in analysis/sgs_audit.jl.

ENV["FLOODA5_NO_VIS"] = "1"
ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
include(joinpath(@__DIR__, "..", "cases", "scripts", "sgs_strip_geometry.jl"))
using Logging, Printf, Statistics, Random
const AG = A5Grid.ArchGDAL

# ── arguments ──────────────────────────────────────────────────────────────
function parse_cli(args)
    pos = String[]; opt = Dict{String,String}()
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            opt[args[i][3:end]] = args[i+1]; i += 2
        else
            push!(pos, args[i]); i += 1
        end
    end
    length(pos) == 3 || error("usage: diag_sgs_tables_audit.jl MESH.parquet DEM.tif OUT.json [options]")
    return pos, opt
end
const (POS, OPT) = parse_cli(ARGS)
const MESHF, DEMF, OUTF = POS
const NSUB = parse(Int, get(OPT, "subset", "500"))
const NHI = parse(Int, get(OPT, "halton-hi", "4096"))
const SEED = parse(Int, get(OPT, "seed", "1"))

const LEVELS = [0.05, 0.1, 0.25, 0.5, 1.0]            # m above the table's z_min; then z_max
const LEVEL_NAMES = ["z_min + 0.05 m", "z_min + 0.10 m", "z_min + 0.25 m", "z_min + 0.50 m", "z_min + 1.00 m", "z_max"]
const MEAN_DEPTHS = [0.01, 0.05, 0.20]                # m of water over the cell area
const AREA_T = [0.01, 0.05, 0.20]                     # basin thresholds: share of the cell's plan area
const DEPTH_T = [0.10, 0.25, 0.50]                    # … and depth below the merge level (m)
const VFLOOR = 0.01                                   # relative volume errors: floor of 1 cm over the cell

# ── mesh, tables, DEM ──────────────────────────────────────────────────────
const mesh = with_logger(() -> load_mesh_geoparquet(MESHF), SimpleLogger(stderr, Logging.Warn))
const n = length(mesh.cells)
const EB = mesh.array_vars["sgs_elev_bins"]; const VC = mesh.array_vars["sgs_vol_curve"]
const AC = mesh.array_vars["sgs_area_curve"]; const SILLS = mesh.array_vars["sgs_edge_sills"]
const AREA = mesh.static_vars["sgs_cell_area"]
const ZMIN = mesh.static_vars["sgs_z_min"]; const ZMAX = mesh.static_vars["sgs_z_max"]
const NBINS = size(EB, 1)
const raster = A5Grid._load_dem_raster(DEMF)
const GEOGRAPHIC = abs(raster.gt[2]) < 1.0                 # the heuristic of build_sgs_tables!
const PX_MODEL = GEOGRAPHIC ? abs(raster.gt[2]) * A5Grid._EARTH_R * π / 180.0 : abs(raster.gt[2])   # its px_size_m (no cos φ)
const PX_M = GEOGRAPHIC ? abs(raster.gt[6]) * A5Grid._EARTH_R * π / 180.0 : abs(raster.gt[2])       # the pixel, from its north–south size
const RE = A5Grid._EARTH_R
dense_h(A) = min(PX_M / 8, sqrt(A) / 64)

"DEM surface at lon/lat points: the model's transform and bilinear lookup."
function sample_dem(lons::Vector{Float64}, lats::Vector{Float64}, xform)
    xs, ys = lons, lats
    if xform !== nothing
        xs = copy(lons); ys = copy(lats); zs = zeros(length(xs))
        AG.transform!(xs, ys, zs, xform)
    end
    z = Vector{Float64}(undef, length(xs))
    for j in eachindex(xs)
        c, r = A5Grid._affine_to_pixel(xs[j], ys[j], raster.gt)
        z[j] = A5Grid._bilinear(raster.band, c, r, raster.nx, raster.ny, raster.nodata_val)
    end
    return z
end

"Regular grid of spacing ≈ h (m) over the cell's bounding box: (lons, lats, linear indices) of the points inside the polygon, and the grid size."
function cell_grid(bnd, h)
    lo0, lo1 = extrema(Float64(v[1]) for v in bnd); la0, la1 = extrema(Float64(v[2]) for v in bnd)
    lam = (la0 + la1) / 2
    nx = max(2, ceil(Int, deg2rad(lo1 - lo0) * RE * cosd(lam) / h))
    ny = max(2, ceil(Int, deg2rad(la1 - la0) * RE / h))
    lons = Float64[]; lats = Float64[]; idx = Int[]
    for j in 1:ny, i in 1:nx
        lo = lo0 + (i - 0.5) / nx * (lo1 - lo0); la = la0 + (j - 0.5) / ny * (la1 - la0)
        if A5Grid._point_in_ring(lo, la, bnd)
            push!(lons, lo); push!(lats, la); push!(idx, (j - 1) * nx + i)
        end
    end
    return lons, lats, idx, nx, ny
end

"The model's Halton candidates for a cell (build_sgs_tables! step 1): sorted finite sample elevations and the number inside the polygon."
function halton_samples(bnd, h2, h3, xform)
    lo0, lo1 = extrema(Float64(v[1]) for v in bnd); la0, la1 = extrema(Float64(v[2]) for v in bnd)
    cl = lo0 .+ h2 .* (lo1 - lo0); ca = la0 .+ h3 .* (la1 - la0)
    inside = [A5Grid._point_in_ring(cl[j], ca[j], bnd) for j in eachindex(cl)]
    z = sample_dem(cl[inside], ca[inside], xform)
    return sort!(filter(isfinite, z)), count(inside)
end

"Signed errors of a table against a reference curve at the audit's levels and volumes (datum: z0, z1)."
function curve_errors(eb, vc, A, ref::FA.HypsCurve, z0, z1)
    lv = vcat(z0 .+ LEVELS, z1)
    rv = map(lv) do η
        Vr = FA.hyps_volume(ref, η)
        (FA.sgs_table_volume(eb, vc, A, η) - Vr) / max(Vr, VFLOOR * A)
    end
    dη = [FA.sgs_table_wse(eb, vc, A, d * A) - FA.hyps_wse(ref, d * A) for d in MEAN_DEPTHS]
    return rv, dη
end
"The same for one reference curve against another."
function curve_errors(c::FA.HypsCurve, ref::FA.HypsCurve, A, z0, z1)
    lv = vcat(z0 .+ LEVELS, z1)
    rv = [(FA.hyps_volume(c, η) - FA.hyps_volume(ref, η)) / max(FA.hyps_volume(ref, η), VFLOOR * A) for η in lv]
    dη = [FA.hyps_wse(c, d * A) - FA.hyps_wse(ref, d * A) for d in MEAN_DEPTHS]
    return rv, dη
end

# ── per-cell columns ───────────────────────────────────────────────────────
const col = Dict{String,Vector{Float64}}()
C(name) = get!(() -> fill(NaN, n), col, name)
function store!(prefix, i, rv, dη)
    for k in eachindex(rv); C("$(prefix)_relV$k")[i] = rv[k]; end
    for k in eachindex(dη); C("$(prefix)_deta$k")[i] = dη[k]; end
end

const strip = haskey(OPT, "strip-params") ? JSON3.read(read(OPT["strip-params"], String), Dict{String,Any}) : nothing
const sp = strip === nothing ? nothing : _ss_params(strip)
strip_xy(lon, lat) = (deg2rad(lon - Float64(strip["lon0"])) * Float64(strip["earth_radius_m"]) * cosd(Float64(strip["lat0"])),
                      deg2rad(lat - Float64(strip["lat0"])) * Float64(strip["earth_radius_m"]))
"Analytic strip terrain on a 0.05 m grid over the cell polygon (as cases/scripts/sgs_strip_tables.jl, finer)."
function strip_points(bnd; h = 0.05)
    P = [strip_xy(Float64(v[1]), Float64(v[2])) for v in bnd]
    x0, x1 = extrema(first.(P)); y0, y1 = extrema(last.(P))
    ring = [[q[1], q[2]] for q in P]
    z = Float64[]
    for x in (x0 + h / 2):h:x1, y in (y0 + h / 2):h:y1
        A5Grid._point_in_ring(x, y, ring) && push!(z, _ss_surface(sp, x, y))
    end
    return z
end

const subset = Set(n <= NSUB ? collect(1:n) : shuffle(MersenneTwister(SEED), collect(1:n))[1:NSUB])
const conv_subset = Set(collect(subset)[1:min(200, length(subset))])
const trees = Vector{Any}(undef, n)                   # per cell: basin summary
const face_dz = Float64[]                             # per face side: sill − the cell's table z_min

function audit_cells(xform)
    h2 = A5Grid._halton(2, 512, 0); h3 = A5Grid._halton(3, 512, 0)
    g2 = A5Grid._halton(2, NHI, 0); g3 = A5Grid._halton(3, NHI, 0)
    t0 = time()
    for i in 1:n
        trees[i] = nothing
        bnd = mesh.cells[i].boundary
        A = AREA[i]
        isfinite(ZMIN[i]) || continue
        eb = view(EB, :, i); vc = view(VC, :, i)
        z0, z1 = ZMIN[i], ZMAX[i]
        # dense reference
        h = dense_h(A)
        lons, lats, idx, nx, ny = cell_grid(bnd, h)
        zd = sample_dem(lons, lats, xform)
        nfin = count(isfinite, zd)
        C("dense_n")[i] = nfin; C("dense_nodata_frac")[i] = 1 - nfin / length(zd)
        nfin == length(zd) || continue                # cells the DEM covers only in part are left out
        ref = FA.HypsCurve(zd, A)
        C("relief")[i] = z1 - z0; C("spacing")[i] = (z1 - z0) / (NBINS - 1)
        C("relief_dense")[i] = ref.z[end] - ref.z[1]
        C("dzmin")[i] = z0 - ref.z[1]; C("dzmax")[i] = z1 - ref.z[end]
        store!("tab", i, curve_errors(eb, vc, A, ref, z0, z1)...)
        # wetted area at the same levels (share of the cell)
        for (k, η) in enumerate(vcat(z0 .+ LEVELS, z1))
            C("tab_dAwet$k")[i] = (FA.sgs_table_wet_area(eb, view(AC, :, i), A, η) - FA.hyps_wet_area(ref, η)) / A
        end
        # level interpolation alone: 100 uniform levels from the dense samples
        ue, uv, _ = FA.uniform_level_table(ref.z, A; n_bins = NBINS)
        store!("uni", i, curve_errors(ue, uv, A, ref, z0, z1)...)
        C("uni_mid_max")[i] = maximum(1:NBINS-1) do k
            η = (ue[k] + ue[k+1]) / 2; Vr = FA.hyps_volume(ref, η)
            abs(FA.sgs_table_volume(ue, uv, A, η) - Vr) / max(Vr, VFLOOR * A)
        end
        C("tab_mid_max")[i] = maximum(1:NBINS-1) do k
            η = (eb[k] + eb[k+1]) / 2; Vr = FA.hyps_volume(ref, η)
            abs(FA.sgs_table_volume(eb, vc, A, η) - Vr) / max(Vr, VFLOOR * A)
        end
        # the model's 512 candidates: points inside, and a rebuild as a check of this script
        zs, nin = halton_samples(bnd, h2, h3, xform)
        C("halton_in")[i] = nin
        if !isempty(zs)
            re, rv, _ = FA.uniform_level_table(zs, A; n_bins = NBINS)
            C("rebuild_dz")[i] = maximum(abs, re .- eb)
            C("rebuild_dV")[i] = maximum(abs, rv .- vc) / max(vc[end], 1e-12)
        end
        # faces: lowest sampled face point against the table's z_min
        s = filter(isfinite, view(SILLS, :, i))
        if !isempty(s)
            C("sill_min_minus_zmin")[i] = minimum(s) - z0
            C("sill_min_minus_dense_zmin")[i] = minimum(s) - ref.z[1]
            append!(face_dz, s .- z0)
        end
        # sampling convergence: 4096 candidates
        if i in subset
            zs4, nin4 = halton_samples(bnd, g2, g3, xform)
            e4, v4, _ = FA.uniform_level_table(zs4, A; n_bins = NBINS)
            store!("hi", i, curve_errors(e4, v4, A, ref, z0, z1)...)
            C("hi_dzmin")[i] = zs4[1] - ref.z[1]; C("hi_dzmax")[i] = zs4[end] - ref.z[end]
            C("hi_in")[i] = nin4
        end
        # is the dense reference itself converged? half the spacing
        if i in conv_subset
            l2, a2, _, _, _ = cell_grid(bnd, h / 2)
            ref2 = FA.HypsCurve(sample_dem(l2, a2, xform), A)
            store!("conv", i, curve_errors(ref, ref2, A, z0, z1)...)
            C("conv_dzmin")[i] = ref.z[1] - ref2.z[1]; C("conv_dzmax")[i] = ref.z[end] - ref2.z[end]
        end
        # the strip's analytic terrain
        if sp !== nothing
            an = FA.HypsCurve(strip_points(bnd), A)
            store!("tab_an", i, curve_errors(eb, vc, A, an, z0, z1)...)
            store!("dense_an", i, curve_errors(ref, an, A, z0, z1)...)
            C("an_dzmin")[i] = z0 - an.z[1]; C("an_dzmax")[i] = z1 - an.z[end]
        end
        # sub-basins
        zg = fill(NaN, nx, ny); zg[idx] .= zd
        t = FA.basin_merge_tree(zg)
        lb = FA.lowest_point_basin(t)
        trees[i] = (n = t.n, nx = nx, ny = ny,
            one = [FA.significant_basins(t; area = a, depth = d) for a in AREA_T, d in DEPTH_T],
            two = [FA.significant_basins(t; area = a, depth = d, two_sided = true) for a in AREA_T, d in DEPTH_T],
            lowest_in_main = lb === nothing,
            low_area = lb === nothing ? NaN : lb.area_frac, low_depth = lb === nothing ? NaN : lb.depth,
            low_main_depth = lb === nothing ? NaN : lb.main_depth)
        i % 5000 == 0 && @printf(stderr, "  %d / %d cells (%.0f s)\n", i, n, time() - t0)
    end
end

if isempty(raster.proj)
    audit_cells(nothing)
else
    AG.createcoordtrans(A5Grid._crs_gis(4326), A5Grid._crs_gis(raster.proj)) do xform
        audit_cells(xform)
    end
end
const ok = [trees[i] !== nothing for i in 1:n]

# ── summaries ──────────────────────────────────────────────────────────────
absq(v) = (w = abs.(filter(isfinite, v)); isempty(w) ? (NaN, NaN, NaN) : (median(w), quantile(w, 0.95), maximum(w)))
sgnmed(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(w))
const out = Dict{String,Any}("mesh" => MESHF, "dem" => DEMF, "n_cells" => n, "n_audited" => count(ok),
    "pixel_m" => PX_M, "n_bins" => NBINS, "levels_m" => LEVELS, "mean_depths_m" => MEAN_DEPTHS,
    "area_thresholds" => AREA_T, "depth_thresholds_m" => DEPTH_T, "argv" => ARGS)

println("\n## ", basename(MESHF), " on ", basename(DEMF))
@printf("%d cells (%d audited), DEM pixel %.2f m, cell area median %.0f m² (√A %.1f m), dense grid h = %.3f m (median %d points per cell; nodata share max %.3f)\n",
        n, count(ok), PX_M, median(AREA), sqrt(median(AREA)), dense_h(median(AREA)),
        round(Int, sgnmed(C("dense_n"))), maximum(filter(isfinite, C("dense_nodata_frac")); init = 0.0))
@printf("build_sgs_tables! takes the pixel as %.2f m, so samples faces every %.2f m
", PX_MODEL, PX_MODEL / 2)
@printf("check of this script: tables rebuilt from the model's 512 Halton candidates differ from the stored ones by at most %.2e m (levels), %.2e (volume, relative to the top)\n",
        maximum(filter(isfinite, C("rebuild_dz")); init = 0.0), maximum(filter(isfinite, C("rebuild_dV")); init = 0.0))

function err_table(title, prefix; sel = trues(n), key = prefix)
    println("\n", title)
    println("| quantity | median abs | p95 abs | max abs | median signed |")
    println("|---|---|---|---|---|")
    d = Dict{String,Any}()
    for k in 1:6
        v = C("$(prefix)_relV$k")[sel]; a = absq(v)
        @printf("| volume at %s (%%) | %.3f | %.3f | %.2f | %+.3f |\n", LEVEL_NAMES[k], 100a[1], 100a[2], 100a[3], 100sgnmed(v))
        d["relV_" * LEVEL_NAMES[k]] = Dict("median_abs" => a[1], "p95_abs" => a[2], "max_abs" => a[3], "median" => sgnmed(v))
    end
    for k in 1:3
        v = C("$(prefix)_deta$k")[sel]; a = absq(v)
        @printf("| level at mean depth %g cm (mm) | %.2f | %.2f | %.1f | %+.2f |\n", 100MEAN_DEPTHS[k], 1e3a[1], 1e3a[2], 1e3a[3], 1e3sgnmed(v))
        d["level_mean_depth_$(MEAN_DEPTHS[k])"] = Dict("median_abs" => a[1], "p95_abs" => a[2], "max_abs" => a[3], "median" => sgnmed(v))
    end
    out[key] = d
end
function dist_table(title, rows)
    println("\n", title)
    println("| quantity | median | p5 | p95 | min | max |")
    println("|---|---|---|---|---|---|")
    for (name, v, scale, key) in rows
        w = filter(isfinite, v) .* scale
        isempty(w) && continue
        @printf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f |\n", name, median(w), quantile(w, 0.05), quantile(w, 0.95), minimum(w), maximum(w))
        out[key] = FA.dist_summary(filter(isfinite, v))
    end
end

err_table("A(i, ii) stored table against the dense reference, all cells", "tab")
println("\nA(i') wetted area, table − dense, share of the cell area (%)")
println("| level | median abs | p95 abs | max abs | median signed |"); println("|---|---|---|---|---|")
for k in 1:6
    v = C("tab_dAwet$k"); a = absq(v)
    @printf("| %s | %.2f | %.2f | %.1f | %+.2f |\n", LEVEL_NAMES[k], 100a[1], 100a[2], 100a[3], 100sgnmed(v))
end
dist_table("A(iii, iv, vi) extremes, sample counts, relief", [
    ("table z_min − dense z_min (m)", C("dzmin"), 1.0, "dzmin"),
    ("table z_max − dense z_max (m)", C("dzmax"), 1.0, "dzmax"),
    ("Halton points inside the polygon (of 512)", C("halton_in"), 1.0, "halton_in"),
    ("relief z_max − z_min, table (m)", C("relief"), 1.0, "relief"),
    ("relief, dense (m)", C("relief_dense"), 1.0, "relief_dense"),
    ("level spacing (mm)", C("spacing"), 1e3, "spacing"),
    ("lowest face point − table z_min, per cell (m)", C("sill_min_minus_zmin"), 1.0, "sill_min_minus_zmin"),
    ("lowest face point − dense z_min, per cell (m)", C("sill_min_minus_dense_zmin"), 1.0, "sill_min_minus_dense_zmin")])
let v = filter(isfinite, C("dzmin")), w = filter(isfinite, C("dzmax")), r = filter(isfinite, C("relief"))
    @printf("\ntable z_min above the dense z_min by more than 1 / 5 / 10 / 25 cm: %.1f / %.1f / %.1f / %.1f %% of cells\n",
            (100 * count(>(t), v) / length(v) for t in (0.01, 0.05, 0.10, 0.25))...)
    @printf("table z_max below the dense z_max by more than 1 / 5 / 10 / 25 cm: %.1f / %.1f / %.1f / %.1f %% of cells\n",
            (100 * count(<(-t), w) / length(w) for t in (0.01, 0.05, 0.10, 0.25))...)
    @printf("relief above 1 / 2 / 5 / 10 m: %.1f / %.1f / %.1f / %.1f %% of cells; level spacing above 1 / 2 / 5 cm: %.1f / %.1f / %.1f %%\n",
            (100 * count(>(t), r) / length(r) for t in (1.0, 2.0, 5.0, 10.0))...,
            (100 * count(>(t * (NBINS - 1)), r) / length(r) for t in (0.01, 0.02, 0.05))...)
    out["dzmin_share_above"] = Dict(string(t) => count(>(t), v) / length(v) for t in (0.01, 0.05, 0.10, 0.25))
    out["dzmax_share_below"] = Dict(string(t) => count(<(-t), w) / length(w) for t in (0.01, 0.05, 0.10, 0.25))
end
let s = filter(isfinite, C("sill_min_minus_zmin")), below = filter(<(-1e-3), s), fb = filter(<(-1e-3), face_dz)
    sh(v, t) = 100 * count(<(-t), v) / max(length(v), 1)
    @printf("\nA(v) cells whose lowest sampled face point is below the table z_min by more than 1 mm / 1 cm / 5 cm / 10 cm / 25 cm: %.1f / %.1f / %.1f / %.1f / %.1f %% of %d cells; among those below by > 1 mm: median %.3f m, p95 %.3f m, max %.3f m\n",
            sh(s, 1e-3), sh(s, 0.01), sh(s, 0.05), sh(s, 0.10), sh(s, 0.25), length(s),
            isempty(below) ? 0.0 : median(-below), isempty(below) ? 0.0 : quantile(-below, 0.95), isempty(below) ? 0.0 : maximum(-below))
    @printf("      face sides (a face seen from one of its cells) with the sill below that cell's table z_min by more than 1 mm / 1 cm / 5 cm / 10 cm: %.1f / %.1f / %.1f / %.1f %% of %d; max %.3f m\n",
            sh(face_dz, 1e-3), sh(face_dz, 0.01), sh(face_dz, 0.05), sh(face_dz, 0.10), length(face_dz), isempty(fb) ? 0.0 : maximum(-fb))
    out["faces_below_zmin"] = Dict("cells_with_faces" => length(s), "cell_deficit_above_1mm" => FA.dist_summary(-below),
        "face_sides_total" => length(face_dz), "face_deficit_above_1mm" => FA.dist_summary(-fb),
        "cell_share_below" => Dict(string(t) => count(<(-t), s) / max(length(s), 1) for t in (1e-3, 0.01, 0.05, 0.10, 0.25)),
        "face_share_below" => Dict(string(t) => count(<(-t), face_dz) / max(length(face_dz), 1) for t in (1e-3, 0.01, 0.05, 0.10)))
end
err_table("A(vi) level interpolation alone: 100 uniform levels built from the dense samples, against the dense curve", "uni")
let a = absq(C("uni_mid_max")), b = absq(C("tab_mid_max"))
    @printf("\nlargest volume error at the mid-points between levels (relative, floor 1 cm): dense-sample table median %.3f %%, p95 %.3f %%, max %.2f %%; stored table median %.3f %%, p95 %.3f %%, max %.2f %%\n",
            100a[1], 100a[2], 100a[3], 100b[1], 100b[2], 100b[3])
    out["uni_mid_max"] = FA.dist_summary(C("uni_mid_max")); out["tab_mid_max"] = FA.dist_summary(C("tab_mid_max"))
end
const insub = [i in subset for i in 1:n]
err_table("A sampling convergence, $(count(insub .& ok)) cells: stored table (512 candidates) against dense", "tab"; sel = insub, key = "tab_subset")
err_table("A sampling convergence, same cells: table rebuilt from $NHI candidates against dense", "hi")
dist_table("extremes on the same cells", [
    ("512: table z_min − dense z_min (m)", C("dzmin")[insub], 1.0, "sub_dzmin"),
    ("$NHI: z_min − dense z_min (m)", C("hi_dzmin"), 1.0, "hi_dzmin"),
    ("512: table z_max − dense z_max (m)", C("dzmax")[insub], 1.0, "sub_dzmax"),
    ("$NHI: z_max − dense z_max (m)", C("hi_dzmax"), 1.0, "hi_dzmax"),
    ("$NHI: points inside the polygon", C("hi_in"), 1.0, "hi_in")])
err_table("reference convergence, $(length(conv_subset)) cells: dense grid at h against h/2", "conv")
dist_table("extremes, h against h/2", [("dense z_min(h) − z_min(h/2) (m)", C("conv_dzmin"), 1.0, "conv_dzmin"),
                                         ("dense z_max(h) − z_max(h/2) (m)", C("conv_dzmax"), 1.0, "conv_dzmax")])
if sp !== nothing
    err_table("strip: stored table against the analytic terrain (0.05 m grid)", "tab_an")
    err_table("strip: dense DEM reference against the analytic terrain", "dense_an")
    dist_table("strip extremes against the analytic terrain", [
        ("table z_min − analytic z_min (m)", C("an_dzmin"), 1.0, "an_dzmin"),
        ("table z_max − analytic z_max (m)", C("an_dzmax"), 1.0, "an_dzmax")])
end

# ── B. sub-basins ──────────────────────────────────────────────────────────
function basin_table(title, sel, key)
    ii = [i for i in 1:n if sel[i] && ok[i]]
    println("\n", title, " (", length(ii), " cells)")
    isempty(ii) && return
    println("| area > | depth > 0.10 m | 0.25 m | 0.50 m | two-sided 0.10 m | 0.25 m | 0.50 m | lowest point in such a basin 0.10 m | 0.25 m | 0.50 m |")
    println("|---|---|---|---|---|---|---|---|---|---|")
    d = Dict{String,Any}()
    for (ia, a) in enumerate(AREA_T)
        one = [100 * count(i -> trees[i].one[ia, id] > 0, ii) / length(ii) for id in eachindex(DEPTH_T)]
        two = [100 * count(i -> trees[i].two[ia, id] > 0, ii) / length(ii) for id in eachindex(DEPTH_T)]
        low = [100 * count(i -> !trees[i].lowest_in_main && trees[i].low_area > a && trees[i].low_depth > DEPTH_T[id], ii) / length(ii)
               for id in eachindex(DEPTH_T)]
        @printf("| %g %% | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |\n", 100a, one..., two..., low...)
        d[string(a)] = Dict("share_pct_one_sided" => one, "share_pct_two_sided" => two, "share_pct_lowest_point_in_it" => low)
    end
    nl = count(i -> !trees[i].lowest_in_main && trees[i].low_area > 0.01 && trees[i].low_depth > 0.02, ii)
    multi = [100 * count(i -> trees[i].one[1, 1] >= m, ii) / length(ii) for m in 1:3]
    @printf("lowest point outside the main (largest-area) basin, in a basin of > 1 %% of the area and > 2 cm deep: %d cells (%.2f %%). Cells with ≥ 1 / 2 / 3 basins at (1 %%, 0.10 m): %.2f / %.2f / %.2f %%\n",
            nl, 100nl / length(ii), multi...)
    d["lowest_outside_main_pct"] = 100nl / length(ii); d["n"] = length(ii)
    out[key] = d
end
@printf("\nB. working grid: %d × %d points across a cell's bounding box (median), spacing %.2f m; DEM pixels across a cell: %.1f\n",
        round(Int, median(trees[i].nx for i in 1:n if ok[i])), round(Int, median(trees[i].ny for i in 1:n if ok[i])),
        dense_h(median(AREA)), sqrt(median(AREA)) / PX_M)
basin_table("B. share of cells (%) with at least one basin beside the main one, all cells", trues(n), "basins_all")

# ── runs: wet cells, C. time step ──────────────────────────────────────────
norm_id(s) = lowercase(lstrip(String(s), '0'))
function run_order(r)
    pos = Dict(norm_id(c.id) => i for (i, c) in enumerate(mesh.cells))
    return [get(pos, norm_id(id), 0) for id in r.mesh.ids]       # run cell → mesh row
end
"4A/P from a 2 × nv lon/lat polygon and the cell's area (the model's cell_length)."
function cell_length(poly::AbstractMatrix, area)
    nv = size(poly, 2); P = 0.0
    for a in 1:nv
        b = a == nv ? 1 : a + 1
        P += A5Grid._haversine_m(poly[1, a], poly[2, a], poly[1, b], poly[2, b])
    end
    return P > 0 ? 4area / P : sqrt(area)
end
function courant_rows(title, r, defs, dt, key)
    ℓ = [cell_length(r.mesh.polygons[i], r.mesh.area[i]) for i in 1:length(r.mesh)]
    println("\n", title)
    @printf("mean step %.4f s (t_end %.0f s / %d steps); cell length 4A/P median %.2f m (min %.2f)\n",
            dt, r.t[end], r.info["summary"]["n_steps"], median(ℓ), minimum(ℓ))
    println("| depth measure | wet cells | depth median (m) | p95 | max | Courant median | p95 | p99 | max | share > 0.7 (%) | share > 1.0 (%) | step this depth requires at Courant 0.7 (s) |")
    println("|---|---|---|---|---|---|---|---|---|---|---|---|")
    d = Dict{String,Any}("dt_mean" => dt, "cell_length_median" => median(ℓ))
    for (name, h) in defs
        e = FA.courant_exposure(h, ℓ, dt)
        w = filter(x -> isfinite(x) && x > 1e-4, h)
        @printf("| %s | %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.3f |\n", name, e.n,
                median(w), quantile(w, 0.95), maximum(w), e.median, e.p95, e.p99, e.max, 100e.share_above_0p7, 100e.share_above_1, e.dt_required)
        d[name] = Dict("n" => e.n, "depth" => FA.dist_summary(w), "courant_median" => e.median, "courant_p95" => e.p95,
                       "courant_p99" => e.p99, "courant_max" => e.max, "share_above_0p7" => e.share_above_0p7,
                       "share_above_1" => e.share_above_1, "dt_required" => e.dt_required)
    end
    out[key] = d
end

if haskey(OPT, "run")
    r = FA.load_run(OPT["run"]); k = FA.nframes(r)
    perm = run_order(r)
    all(>(0), perm) || error("the run's cells are not those of the mesh")
    h = FA.depth(r, k); V = FA.volume(r, k); η = FA.wse(r, k)
    zmin = r.mesh.static["sgs_z_min"]; zmax = r.mesh.static["sgs_z_max"]
    wet = [isfinite(h[i]) && h[i] > 1e-4 for i in eachindex(h)]
    wet_mesh = falses(n); for i in eachindex(h); wet[i] && (wet_mesh[perm[i]] = true); end
    @printf("\nrun %s: last frame t = %.0f s, %d wet cells (depth > 1e-4 m), %.2f Mm³; water above z_max in %.1f %% of wet cells\n",
            basename(OPT["run"]), r.t[k], count(wet), sum(V[wet]) / 1e6, 100 * count(i -> wet[i] && η[i] > zmax[i], eachindex(h)) / count(wet))
    basin_table("B. the same, cells wet at the last frame of $(basename(OPT["run"]))", wet_mesh, "basins_wet")
    err_table("A(i, ii) stored table against dense, cells wet at the last frame", "tab"; sel = wet_mesh, key = "tab_wet")
    capped = [wet[i] ? min(h[i], max(0.0, zmax[i] - zmin[i])) : 0.0 for i in eachindex(h)]
    hyd = map(eachindex(h)) do i
        wet[i] || return 0.0
        m = perm[i]
        Aw = FA.sgs_table_wet_area(view(EB, :, m), view(AC, :, m), AREA[m], η[i])
        Aw > 0 ? V[i] / Aw : 0.0
    end
    dt = r.t[end] / r.info["summary"]["n_steps"]
    courant_rows("C. Courant numbers at the last frame of $(basename(OPT["run"])) (SGS)", r,
        [("capped: min(η − z_min, z_max − z_min) [the model's]", capped), ("η − z_min", [wet[i] ? h[i] : 0.0 for i in eachindex(h)]),
         ("V / A_wet", hyd)], dt, "courant_sgs")
    out["run"] = Dict("path" => OPT["run"], "t" => r.t[k], "wet_cells" => count(wet), "argv" => r.info["argv"])
end
if haskey(OPT, "compare-run")
    r = FA.load_run(OPT["compare-run"]); k = FA.nframes(r)
    h = FA.depth(r, k)
    courant_rows("C. Courant numbers at the last frame of $(basename(OPT["compare-run"])) ($(r.flow_model))", r,
        [("water depth", [isfinite(x) ? x : 0.0 for x in h])], r.t[end] / r.info["summary"]["n_steps"], "courant_compare")
end

# ── embankment: the cells on it ────────────────────────────────────────────
if haskey(OPT, "embankment-params")
    ep = JSON3.read(read(OPT["embankment-params"], String), Dict{String,Any})
    X(lon) = (lon - ep["lon0"]) / (ep["lon1"] - ep["lon0"]) * ep["domain_m_ew"]
    Y(lat) = (lat - ep["lat0"]) / (ep["lat1"] - ep["lat0"]) * ep["domain_m_ns"]
    xw = ep["emb_x_m"] - ep["emb_width_m"] / 2; xe = ep["emb_x_m"] + ep["emb_width_m"] / 2
    yn0 = ep["notch_y_centre_m"] - ep["notch_width_m"] / 2; yn1 = ep["notch_y_centre_m"] + ep["notch_width_m"] / 2
    rows = Any[]
    for i in 1:n
        ok[i] || continue
        b = mesh.cells[i].boundary
        x0, x1 = extrema(X(Float64(v[1])) for v in b); y0, y1 = extrema(Y(Float64(v[2])) for v in b)
        (x1 > xw && x0 < xe) || continue
        cls = (x0 < xw && x1 > xe) ? "spans" : x0 < xw ? "west face" : x1 > xe ? "east face" : "on crest"
        push!(rows, Dict("cell" => i, "id" => mesh.cells[i].id, "x0" => x0, "x1" => x1, "y0" => y0, "y1" => y1, "class" => cls,
            "notch" => (y1 > yn0 && y0 < yn1), "relief" => ZMAX[i] - ZMIN[i], "z_min" => ZMIN[i], "z_max" => ZMAX[i],
            "basins_1pct_0p1m" => trees[i].one[1, 1], "basins_5pct_0p25m" => trees[i].one[2, 2],
            "two_sided_1pct_0p1m" => trees[i].two[1, 1], "lowest_in_main" => !(!trees[i].lowest_in_main && trees[i].low_area > 0.01 && trees[i].low_depth > 0.02)))
    end
    @printf("\nB. embankment x = %.0f–%.0f m (crest %.1f m), notch y = %.0f–%.0f m: %d cells touch it\n", xw, xe, ep["emb_crest_z"], yn0, yn1, length(rows))
    println("| class | cells | of which on the notch | x-extent median (m) | relief median (m) | with a basin (1 %, 0.10 m) | (5 %, 0.25 m) | two-sided (1 %, 0.10 m) | lowest point outside main (basin > 1 %, > 2 cm) |")
    println("|---|---|---|---|---|---|---|---|---|")
    for cls in ("spans", "west face", "east face", "on crest")
        rr = filter(q -> q["class"] == cls, rows)
        isempty(rr) && (println("| $cls | 0 | | | | | | | |"); continue)
        @printf("| %s | %d | %d | %.1f | %.2f | %d | %d | %d | %d |\n", cls, length(rr), count(q -> q["notch"], rr),
                median(q["x1"] - q["x0"] for q in rr), median(q["relief"] for q in rr),
                count(q -> q["basins_1pct_0p1m"] > 0, rr), count(q -> q["basins_5pct_0p25m"] > 0, rr),
                count(q -> q["two_sided_1pct_0p1m"] > 0, rr), count(q -> !q["lowest_in_main"], rr))
    end
    nn = filter(q -> q["notch"], rows)
    println("cells on the notch: ", join(["#$(q["cell"]) x $(round(q["x0"], digits = 1))–$(round(q["x1"], digits = 1)) ($(q["class"])), z $(round(q["z_min"], digits = 2))–$(round(q["z_max"], digits = 2)), basins $(q["basins_1pct_0p1m"])" for q in nn], "; "))
    sb = filter(q -> q["basins_1pct_0p1m"] > 0, rows)
    println("cells on the embankment with a second basin (1 %, 0.10 m): ", isempty(sb) ? "none" :
            join(["#$(q["cell"]) x $(round(q["x0"], digits = 1))–$(round(q["x1"], digits = 1)) y $(round(q["y0"], digits = 1))–$(round(q["y1"], digits = 1)) ($(q["class"]))" for q in sb], "; "))
    @printf("largest x-extent of any cell in the mesh: %.1f m (embankment width %.0f m)\n",
            maximum(maximum(X(Float64(v[1])) for v in c.boundary) - minimum(X(Float64(v[1])) for v in c.boundary) for c in mesh.cells), xe - xw)
    out["embankment_cells"] = rows
end

mkpath(dirname(abspath(OUTF)))
open(OUTF, "w") do f; JSON3.write(f, out); end
println("\nwrote ", OUTF)
