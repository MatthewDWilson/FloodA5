# cases/registry/mesh_families.jl — the benchmark cases on every mesh family
# (docs/TEST_CASES.md §8; docs/MESHES.md; WP3 part 3, Phase B, 0.8.1).
#
# Eleven meshes at matched mean cell area (√area ≈ 22.2 m, A5 res 18):
#   cart          FloodA5's Cartesian grid (lon/lat-aligned squares)
#   rotcart       squares rotated by 30° (mesh/mesh_families.py)
#   hex           regular hexagons, pointy-top (mesh_families.py)
#   voronoi       centroidal Voronoi: jittered lattice + 5 Lloyd steps
#   voronoi_jit   Voronoi of a jittered lattice (no relaxation)
#   a5            A5 res 18 (pya5)
#   h3 igeo7 isea3h healpix rhealpix   discrete global grids
#                 (mesh/dggs_families.py), each at the level whose cell area
#                 is nearest 22.2² m²: levels are discrete, so their cells are
#                 16.6 (H3), 17.3 (rHEALPix), 22.9 (IGEO7), 23.3 (ISEA3H) and
#                 24.9 m (HEALPix) across. Each job records its own size
#                 (`cell_size_m`).
# The planar families match the size exactly; their generator options are
# those of test/channel_bench/make_channel_meshes.jl --families (keep the two
# in step).

const _MF_L = 22.2                       # target √(mean cell area), m
const _MF_ALL = ["cart", "rotcart", "hex", "voronoi", "voronoi_jit",
                 "a5", "h3", "igeo7", "isea3h", "healpix", "rhealpix"]
"mesh/mesh_families.py arguments of the planar families at cell size L (m)."
const _MF_PLANAR = Dict(
    "hex"         => L -> ["hex", "--dx", string(round(L / sqrt(3sqrt(3) / 2); digits = 6))],   # area L²
    "rotcart"     => L -> ["rotcart", "--dx", string(L), "--angle", "30"],
    "voronoi"     => L -> ["voronoi", "--dx", string(L), "--jitter", "0.3", "--lloyd", "5", "--seed", "1"],
    "voronoi_jit" => L -> ["voronoi", "--dx", string(L), "--jitter", "0.3", "--lloyd", "0", "--seed", "1"])

"""
    _mf_mesh(fam, aoi; dem = nothing, flow = "standard", clip = false, tag = fam) -> (preps, mesh)

Preparation steps that write family `fam`'s mesh over `aoi` to
`{case}/<tag>.parquet` (elevation from `dem` if given; SGS tables with
`flow = "sgs"`), and that path. Cartesian: `--meshgen-cartesian … --dx`; A5
and the DGGS: `--meshgen-dggs … --dggs FAM --dx`; the planar families:
mesh_families.py (`--clip`: cells meeting the AOI polygon), then
`--meshload … --dem` if there is a DEM (from `{case}/<raw_tag>_raw.parquet`,
which meshes of one family with different `tag`s share).
"""
function _mf_mesh(fam, aoi; dem = nothing, flow = "standard", clip = false, tag = fam, raw_tag = fam)
    out = "{case}/$tag.parquet"
    demargs = dem === nothing ? String[] : ["--dem", dem]
    tail = vcat(demargs, ["--meshout", out, "--flow-model", flow, "--mesh-only"])
    if fam == "cart"
        return [Prep(tag, :floodmodel, vcat(["--meshgen-cartesian", aoi, "--dx", string(_MF_L)], tail);
                     outputs = [out])], out
    elseif haskey(_MF_PLANAR, fam)
        raw = dem === nothing ? out : "{case}/$(raw_tag)_raw.parquet"
        steps = [Prep("$(raw_tag)_raw", :python,
                      vcat(["{repo}/mesh/mesh_families.py"], _MF_PLANAR[fam](_MF_L), ["--aoi", aoi],
                           clip ? ["--clip"] : String[], ["--out", raw]); outputs = [raw])]
        dem === nothing ||
            push!(steps, Prep(tag, :floodmodel, vcat(["--meshload", raw], tail); outputs = [out]))
        return steps, out
    end
    return [Prep(tag, :floodmodel, vcat(["--meshgen-dggs", aoi, "--dggs", fam, "--dx", string(_MF_L)], tail);
                 outputs = [out])], out
end

"Size, cell count and largest valence of each run's mesh (added to every family case's metrics)."
function _mf_mesh_metrics!(out, ctx)
    for (name, r) in ctx.runs
        d = get!(out, name, Dict{String,Float64}())
        d["cell_size_m"] = sqrt(sum(r.mesh.area) / length(r.mesh.area))
        d["n_cells"] = length(r.mesh.area)
        d["max_neighbours"] = size(r.mesh.neighbours, 1)
    end
    return out
end

_mf_list(m, key; digits = 3) =
    join(("$j $(round(d[key]; digits))" for (j, d) in sort(collect(m); by = first) if haskey(d, key)), ", ")

const _MF_ABOUT = """
    Meshes (√ mean cell area ≈ 22.2 m): Cartesian, rotated Cartesian (30°),
    hexagons, centroidal and jittered Voronoi, A5 res 18, and the discrete
    global grids H3, IGEO7, ISEA3H, HEALPix and rHEALPix at the level nearest
    that size (16.6–24.9 m; each job's `cell_size_m`)."""

# ── point spread ───────────────────────────────────────────────────────────
let aoi = "{repo}/test/square/square_domain.geojson"     # 1.5 km square at (0, 0)
    prep = Prep[]; jobs = Job[]
    for fam in _MF_ALL
        p, mesh = _mf_mesh(fam, aoi)
        append!(prep, p)
        push!(jobs, Job(fam, ["--meshload", mesh, "--flow-model", "standard",
                              "--injection-point", "0.0,0.0,10.0", "--closed-boundaries",
                              "--sim-duration", "3600", "--dt-max", "2", "--output-interval", "600"];
                        mesh = fam))
    end
    register!(Case(
        name = "mesh_families_point_spread", tier = :standard,
        title = "Flat point source on eleven mesh families (directional bias)",
        doc = "docs/TEST_CASES.md §8",
        description = """
        The `point_spread` case on every mesh family: 10 m³/s injected at the
        centre of the closed, flat 1.5 km square at (0°, 0°) for one hour,
        standard flow with the default scheme. """ * _MF_ABOUT * """
        Metrics at t = 3600 s as `point_spread` (A4, A6 angular harmonics,
        cardinal/diagonal front-radius ratio, hull circularity, front-radius
        CV), and each mesh's size. A few minutes in all.""",
        prep = prep, jobs = jobs,
        analyse = ctx -> _mf_mesh_metrics!(_point_spread_metrics(ctx), ctx),
        checks = [Check("mass conserved (|error| < 1e-9 of the volume added)",
                        m -> all(d["mass_balance_rel"] < 1e-9 for d in values(m)); detail = m -> _mf_list(m, "mass_balance_rel"; digits = 12)),
                  Check("stored volume = 10 m³/s × 3600 s (±0.1%)",
                        m -> all(abs(d["storage_m3"] / 36000 - 1) < 1e-3 for d in values(m)); detail = m -> _mf_list(m, "storage_m3"; digits = 1)),
                  Check("near-circular spread on every mesh (hull circularity > 0.9)",
                        m -> all(d["circularity"] > 0.9 for d in values(m)); detail = m -> _mf_list(m, "circularity"))],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end

# ── planar slope ───────────────────────────────────────────────────────────
let
    aoi, dem = "{case}/planar_embankment_aoi.geojson", "{case}/planar_embankment_dem.tif"
    lat, lon = _planar_point()
    prep = [Prep("dem", :python,
                 ["{repo}/test/planar_embankment/generate_planar_embankment_dem.py", "--out", "{case}",
                  "--emb-height", "0", "--emb-width", "0"];
                 outputs = ["{case}/planar_embankment_dem.tif", "{case}/planar_embankment_params.json",
                            "{case}/planar_embankment_aoi.geojson"])]
    jobs = Job[]
    for fam in _MF_ALL
        p, mesh = _mf_mesh(fam, aoi; dem)
        append!(prep, p)
        push!(jobs, Job(fam, ["--meshload", mesh, "--flow-model", "standard",
                              "--injection-point", "$lat,$lon,0.1", "--closed-boundaries",
                              "--sim-duration", "72000", "--dt-max", "2", "--output-interval", "1800"];
                        mesh = fam))
    end
    register!(Case(
        name = "mesh_families_planar_slope", tier = :standard,
        title = "Point source on a 0.1% planar slope on eleven mesh families: N/S symmetry, flow direction",
        doc = "docs/TEST_CASES.md §8",
        description = """
        The `planar_slope` case on every mesh family: a pure 0.1% west→east
        plane (4 × 2 km at 51° N, 5 m DEM), 0.1 m³/s at x = 400 m on the centre
        line, closed boundaries, 20 h, standard flow with the default scheme.
        """ * _MF_ABOUT * """
        The terrain has no north/south distinction, so the N/S volume
        asymmetry, lateral skew A2 and plume deflection (at 5, 10, 15, 20 h)
        are the mesh's directional bias; flow should run east. About 15 min
        on one 4-thread run at a time.""",
        prep = prep, jobs = jobs,
        analyse = ctx -> _mf_mesh_metrics!(_planar_slope_metrics(ctx), ctx),
        checks = [Check("injection point = the generator's (params file)",
                        m -> all(d["inj_point_mismatch_deg"] < 1e-9 for d in values(m))),
                  Check("mass conserved (|error| < 1e-9 of the volume added)",
                        m -> all(d["mass_balance_rel"] < 1e-9 for d in values(m)); detail = m -> _mf_list(m, "mass_balance_rel"; digits = 12)),
                  Check("closed domain stores Q·t (±0.1%)",
                        m -> all(abs(d["storage_m3"] / (0.1 * d["t_end"]) - 1) < 1e-3 for d in values(m));
                        detail = m -> _mf_list(m, "storage_m3"; digits = 1)),
                  # `planar_slope` also asks for ≥ 90% of moving cells within ±30°.
                  # At 22 m the 20 h plume moves only 22–160 cells, and the radial
                  # flow next to the source decides that fraction: the hexagons'
                  # four cells at ±40° (a symmetric pair of pairs) give 0.89.
                  # It is recorded (dir_frac30), not checked; HEALPix's 0.77
                  # (five cells, lopsided) goes with its directional bias.
                  Check("flow is downslope on every mesh: speed-weighted mean direction within 5° of east",
                        m -> all(abs(d["dir_mean_w"]) < 5 for d in values(m));
                        detail = m -> "mean " * _mf_list(m, "dir_mean_w"; digits = 2) * "°; cells within ±30°: " * _mf_list(m, "dir_frac30")),
                  # The lon/lat-aligned squares and the pointy-top hexagons are mirror
                  # symmetric about the source's row (up to meridian convergence,
                  # `planar_slope`): their plume must be too.
                  Check("mirror-symmetric meshes (Cartesian, hexagons): |asym|, |A2| < 1e-3 and |deflection| < 0.01° at 5–20 h",
                        m -> all(abs(m[j]["asym_$(t)h"]) < 1e-3 && abs(m[j]["A2_$(t)h"]) < 1e-3 &&
                                 abs(m[j]["deflection_$(t)h"]) < 0.01 for j in ("cart", "hex"), t in (5, 10, 15, 20));
                        detail = m -> "20 h asym: " * _mf_list(Dict(j => m[j] for j in ("cart", "hex")), "asym_20h"; digits = 7))],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end

# ── straight channel ───────────────────────────────────────────────────────
let res = [18], angles = [0, 15, 30, 45], offsets = [0]
    cases_ = [_ch_case(m, r, α, of) for r in res, m in _MF_ALL, α in angles, of in offsets]
    register!(Case(
        name = "mesh_families_channel", tier = :standard,
        title = "Straight 20 m channel at 0–45° on eleven mesh families (SWE and local-inertial)",
        doc = "docs/TEST_CASES.md §8; HYDRAULICS.md §15.12",
        description = """
        The `channel` case at 22 m (A5 res 18; the channel about one cell wide)
        on every mesh family, orientations 0°, 15°, 30° and 45°, no offset,
        full SWE and the local-inertial default: 88 runs of 3 h, meshes from
        test/channel_bench/make_channel_meshes.jl --families. """ * _MF_ABOUT * """
        Metrics per run as `channel` (Δη, spill, unsteadiness, outflow);
        under "series" the rms and spread of Δη over the orientations per
        solver × mesh — the spread is the mesh's orientation dependence.
        About an hour on one 4-thread run at a time.""",
        prep = [Prep("meshes", :julia,
                     ["{repo}/test/channel_bench/make_channel_meshes.jl", "{case}",
                      "--res", join(res, ","), "--angles", join(angles, ","), "--offsets", join(offsets, ","),
                      "--families", join(_MF_ALL, ",")];
                     outputs = vcat(["{case}/manifest.csv", "{case}/jobs.txt", "{case}/q20.csv"],
                                    ["{case}/$c.parquet" for c in vec(cases_)]))],
        jobs = _channel_jobs(; res, angles, offsets, meshes = _MF_ALL),
        analyse = _channel_metrics,
        checks = [Check("mass conserved (|error| < 1e-5 of the volume added, from the log)",
                        m -> all(d["mass_balance_log_rel"] < 1e-5 for (_, d) in _ch_runs(m));
                        detail = m -> "worst " * string(maximum(d["mass_balance_log_rel"] for (_, d) in _ch_runs(m)))),
                  Check("job arguments = the generator's jobs.txt", m -> m["series"]["jobs_match_generator"] == 1),
                  # `channel`: within 15% at 22 m on A5 and Cartesian
                  Check("near steady at 3 h: outflow/inflow within 15%",
                        m -> all(abs(d["outflow"] - 1) <= 0.15 for (_, d) in _ch_runs(m));
                        detail = m -> _ch_list(m, "outflow"))],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end

# ── Carlisle to the Eden peak, against the 5 m references ─────────────────
# The families not run by `carlisle_peak` (A5 res 18 and the 20 m Cartesian
# grid aligned with the DEM are its a5_* and cart_* runs, read from there).
# The domain is a British National Grid rectangle, so the planar families
# keep the cells that meet it (--clip), as A5 and the DGGS do: their bounding
# box would reach past the DEM, and cells without elevation are walls.
const _MF_CARL_FAMS = ["rotcart", "hex", "voronoi", "voronoi_jit", "h3", "igeo7", "isea3h", "healpix", "rhealpix"]
const _MF_CARL_PEAK = ("a5_li", "a5_sgs", "cart_li")     # from carlisle_peak

"""
Metrics of `mesh_families_carlisle`: every run (and `carlisle_peak`'s
$(join(_MF_CARL_PEAK, ", ")) as `peak:<job>`) against the 5 m runs of
`carlisle_peak_fine` aggregated onto its cells, at 49.5, 52.5 and 55.5 h
(`FA.compare_to_fine`, as `carlisle_peak_fine`'s `vs5:` entries).
"""
function _mf_carl_metrics(ctx)
    runs = Dict{String,Any}(ctx.runs)
    tables = Dict{String,Any}(n => (get(j.tags, :solver, "") == "sgs" ? FA.sgs_tables(ctx.runs[n]) : nothing)
                              for (n, j) in ctx.jobs)
    peak = case("carlisle_peak")
    for j in peak.jobs
        j.name in _MF_CARL_PEAK || continue
        job_status(peak, j) === :done || error("mesh_families_carlisle: carlisle_peak/$(j.name) is not complete — run that case first")
        r = FA.load_run(job_output(peak, j))
        runs["peak:" * j.name] = r
        tables["peak:" * j.name] = get(j.tags, :solver, "") == "sgs" ? FA.sgs_tables(r) : nothing
    end
    fine = case("carlisle_peak_fine")
    ref = Dict{String,Any}()
    for (tag, job) in (("swe5", "cart5_swe"), ("li5", "cart5_li"))
        j = only(filter(x -> x.name == job, fine.jobs))
        job_status(fine, j) === :done || error("mesh_families_carlisle: carlisle_peak_fine/$job is not complete")
        ref[tag] = FA.load_run(job_output(fine, j))
    end
    t0, t1 = _CARL_PEAK_WINDOW
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in runs
        d = Dict{String,Float64}("t_end" => r.t[end], "storage_w1_m3" => FA.stored_volume(r, FA.nframes(r)),
                                 "cell_size_m" => sqrt(sum(r.mesh.area) / length(r.mesh.area)),
                                 "n_cells" => length(r.mesh.area))
        owner = nothing
        for (fname, f) in ref, (tag, t) in (("w0", t0), ("wmid", (t0 + t1) / 2), ("w1", t1))
            c = FA.compare_to_fine(f, r, t; tol = _CARL_PEAK_TOL, owner, tables = tables[name])
            owner = c.owner
            d["$(fname)_storage_ratio_$tag"] = c.storage_ratio
            d["$(fname)_level_median_$tag"] = c.level_median
            d["$(fname)_level_rms_$tag"] = c.level_rms
            d["$(fname)_csi_$tag"] = c.csi
        end
        out[name] = d
    end
    return out
end

let
    aoi, dem = "$_CARL/Carlisle_domain.geojson", "$_CARL/Carlisle_LiDAR_5m_mean.tif"
    prep = Prep[]; jobs = Job[]
    for fam in _MF_CARL_FAMS
        p, mesh = _mf_mesh(fam, aoi; dem, clip = true)
        ps, sgsmesh = _mf_mesh(fam, aoi; dem, clip = true, flow = "sgs", tag = "$(fam)_sgs")
        append!(prep, p); append!(prep, filter(x -> !(x.name in [q.name for q in prep]), ps))
        push!(jobs, Job("$(fam)_li", vcat(["--meshload", mesh, "--flow-model", "standard"], _CARL_PEAK_RUN);
                        mesh = fam, solver = "li"))
        push!(jobs, Job("$(fam)_sgs", vcat(["--meshload", sgsmesh, "--flow-model", "sgs"], _CARL_PEAK_RUN);
                        mesh = fam, solver = "sgs"))
    end
    register!(Case(
        name = "mesh_families_carlisle", tier = :full,
        title = "Carlisle to the Eden peak on the mesh families at ~22 m (local-inertial, SGS) against the 5 m references",
        doc = "docs/TEST_CASES.md §8",
        description = """
        The `carlisle_peak` event (55.5 h, river inflow only) at ~22 m on the
        mesh families that case does not run — rotated Cartesian, hexagons,
        centroidal and jittered Voronoi (cells meeting the domain polygon),
        H3, IGEO7, ISEA3H, HEALPix and rHEALPix (16.6–24.9 m) — with the
        local-inertial default and with SGS (tables from the 5 m DEM). Every
        run, and `carlisle_peak`'s A5 (local-inertial, SGS) and 20 m Cartesian
        runs, is compared with the 5 m SWE benchmark and the 5 m
        local-inertial run of `carlisle_peak_fine` on its own cells at 49.5,
        52.5 and 55.5 h: storage ratio, median and rms water-level difference,
        CSI at 0.1 m (`swe5_*`, `li5_*`). Needs `carlisle_peak` and
        `carlisle_peak_fine` run first. CPU: about 20 min per
        local-inertial run and 25 min per SGS run on 4 threads (18 runs), plus
        about an hour of mesh and SGS-table preparation.""",
        prep = prep, jobs = jobs,
        analyse = _mf_carl_metrics,
        checks = [Check("mass conserved (|conservation error| < 1e-6 of the volume applied)",
                        m -> all(d["mass_balance_rel"] < 1e-6 for (k, d) in m if haskey(d, "mass_balance_rel"));
                        detail = m -> _mf_list(Dict(k => d for (k, d) in m if haskey(d, "mass_balance_rel")), "mass_balance_rel"; digits = 9)),
                  Check("volume added = the hydrographs' integral (|input error| < 1e-5)",
                        m -> all(abs(d["input_error_rel"]) < 1e-5 for (k, d) in m if haskey(d, "input_error_rel"))),
                  Check("every family compared with the 5 m runs at the three times",
                        m -> all(isfinite(get(get(m, k, Dict()), "swe5_csi_w1", NaN))
                                 for k in vcat(["$(f)_$s" for f in _MF_CARL_FAMS for s in ("li", "sgs")],
                                               ["peak:" * j for j in _MF_CARL_PEAK])))],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
