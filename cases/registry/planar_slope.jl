# cases/registry/planar_slope.jl — the planar N/S volume-symmetry benchmark
# (docs/TEST_CASES.md §2b, §5c). A point source on a uniform 0.1% west→east
# slope: the terrain has no north/south distinction, so any N/S asymmetry of
# the plume is a mesh/solver artefact.

_ps_all(f, m) = all(f(d) for d in values(m))
_ps_list(m, key; digits = 4) =
    join(("$j $(round(d[key]; digits))" for (j, d) in sort(collect(m); by = first) if haskey(d, key)), ", ")

const _PS_LAT, _PS_LON = _planar_point()
const _PS_Q = 0.1                       # m³/s, as the historical runs (HYDRAULICS.md §15.1)
const _PS_T = 72000.0                   # s
const _PS_TIMES = (18000.0, 36000.0, 54000.0, 72000.0)
# The N/S split excludes |y| ≤ 37.5 m about the source cell (FA.ns_split). Not
# the historical 50 m: on the 25 m Cartesian grid ±50 m falls exactly on cell
# centres, and rounding then puts one row "south" and its mirror image in the
# band (asym −0.11 at 20 h for a plume symmetric to A2 ≈ 3e-6); 37.5 m lies
# between rows.
const _PS_BAND = 37.5

function _planar_slope_metrics(ctx)
    p = _planar_params(ctx.dir)
    lat0, lon0 = Float64(p["injection_lat"]), Float64(p["injection_lon"])
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in ctx.runs
        d = Dict{String,Float64}("inj_point_mismatch_deg" => max(abs(lat0 - _PS_LAT), abs(lon0 - _PS_LON)))
        for t in _PS_TIMES
            k = FA.frame_at(r, t; tol = 30.0)   # frames fall on the first step past each output time
            pm = FA.planar_metrics(r, k; lon0, lat0, band = _PS_BAND)
            tag = string(round(Int, t / 3600)) * "h"
            d["asym_$tag"] = pm.asym
            d["A2_$tag"] = pm.A2
            d["deflection_$tag"] = pm.deflection
        end
        k = FA.nframes(r)
        fd = FA.flow_direction_metrics(r, k)
        d["dir_mean_w"] = fd.mean_dir_w
        d["dir_std_w"] = fd.circ_std_w
        d["dir_frac30"] = fd.frac_within_30
        d["storage_m3"] = FA.storage(r, k)
        d["t_end"] = r.t[k]
        out[name] = d
    end
    return out
end

let
    aoi, dem = "{case}/planar_embankment_aoi.geojson", "{case}/planar_embankment_dem.tif"
    common = ["--flow-model", "standard",
              "--injection-point", "$(_PS_LAT),$(_PS_LON),$(_PS_Q)", "--closed-boundaries",
              "--sim-duration", string(round(Int, _PS_T)), "--dt-max", "2", "--output-interval", "1800"]
    jobs = [Job("a5_res16", vcat(["--meshload", "{case}/a5_res16.parquet"], common); mesh = "a5", res = 16),
            Job("a5_res18", vcat(["--meshload", "{case}/a5_res18.parquet"], common); mesh = "a5", res = 18),
            Job("cart_dx25", vcat(["--meshload", "{case}/cart_dx25.parquet"], common); mesh = "cartesian", dx = 25)]
    prep = [Prep("dem", :python,
                 ["{repo}/test/planar_embankment/generate_planar_embankment_dem.py", "--out", "{case}",
                  "--emb-height", "0", "--emb-width", "0"];
                 outputs = ["{case}/planar_embankment_dem.tif", "{case}/planar_embankment_params.json",
                            "{case}/planar_embankment_aoi.geojson"]),
            [Prep("a5_res$res", :floodmodel,
                  ["--meshgen", aoi, "--meshres", "$res", "--dem", dem,
                   "--meshout", "{case}/a5_res$res.parquet", "--flow-model", "standard", "--mesh-only"];
                  outputs = ["{case}/a5_res$res.parquet"]) for res in (16, 18)]...,
            Prep("cart_dx25", :floodmodel,
                 ["--meshgen-cartesian", aoi, "--dx", "25", "--dem", dem,
                  "--meshout", "{case}/cart_dx25.parquet", "--mesh-only"];
                 outputs = ["{case}/cart_dx25.parquet"])]
    register!(Case(
        name = "planar_slope", tier = :standard,
        title = "Point source on a 0.1% planar slope: N/S volume symmetry",
        doc = "docs/TEST_CASES.md §2b",
        description = """
        A pure 0.1% west→east plane (4 × 2 km, 5 m DEM; the embankment generator
        with --emb-height 0 --emb-width 0), 0.1 m³/s injected at x = 400 m on the
        centre line, closed boundaries, 20 h (72 000 s, Δt ≤ 2 s), standard flow
        with the default scheme, on A5 res 16, 18 and a 25 m Cartesian grid. The
        terrain has no north/south distinction, so the N/S volume asymmetry
        (asym), the lateral skew A2 and the plume deflection are mesh/solver
        artefacts; flow directions should be downslope (east). Metrics at 5, 10,
        15 and 20 h (FA.planar_metrics, about the source cell; ±37.5 m band
        excluded) and flow-direction statistics at 20 h. The historical duration
        and Δt cap (HYDRAULICS.md §15.1) are kept: 36 000 steps cost 20 s (res 16),
        65 s (Cartesian) and 90 s (res 18) of simulation on 2–3 threads.""",
        prep = prep, jobs = jobs,
        analyse = _planar_slope_metrics,
        checks = [
            Check("injection point = the generator's (params file)",
                  m -> _ps_all(d -> d["inj_point_mismatch_deg"] < 1e-9, m)),
            Check("mass conserved (|error| < 1e-9 of the volume added)",
                  m -> _ps_all(d -> d["mass_balance_rel"] < 1e-9, m); detail = m -> _ps_list(m, "mass_balance_rel"; digits = 12)),
            Check("closed domain stores Q·t = $(_PS_Q * _PS_T) m³ (±0.1%)",
                  m -> _ps_all(d -> abs(d["storage_m3"] / (_PS_Q * d["t_end"]) - 1) < 1e-3, m);
                  detail = m -> _ps_list(m, "storage_m3"; digits = 1)),
            # The 25 m grid is mirror-symmetric about the source row up to
            # meridian convergence (a fixed dlon is ~6e-5 narrower 300 m north
            # than south) and N/S walls 12.5 m off-centre that the plume never
            # reaches; asym ~1e-5 is that. A scheme defect would show far above 1e-3.
            Check("Cartesian: symmetric about the source row (|asym|, |A2| < 1e-3, |deflection| < 0.01° at 5–20 h)",
                  m -> all(abs(m["cart_dx25"]["asym_$(t)h"]) < 1e-3 && abs(m["cart_dx25"]["A2_$(t)h"]) < 1e-3 &&
                           abs(m["cart_dx25"]["deflection_$(t)h"]) < 0.01 for t in (5, 10, 15, 20));
                  detail = m -> "20 h: asym $(m["cart_dx25"]["asym_20h"]), A2 $(m["cart_dx25"]["A2_20h"])"),
            # A5 has no mirror symmetry; its residual bias is the mesh's. Before WP1
            # (legacy gradient, component friction) the 20 h asymmetry was 0.68
            # (res 16) and 0.87 (res 18) with 14–15° deflection; with the diamond
            # gradient and vector friction of the 0.3.0 default it is ≤ 0.1 and
            # ≤ 3° here (≈ 0.19 recorded at res 18 for another source position,
            # HYDRAULICS.md §15.8). The bounds separate the two regimes.
            Check("A5: residual directional bias well below the pre-WP1 level (|asym| < 0.2, |deflection| < 5° at 5–20 h)",
                  m -> all(abs(m[j]["asym_$(t)h"]) < 0.2 && abs(m[j]["deflection_$(t)h"]) < 5
                           for j in ("a5_res16", "a5_res18"), t in (5, 10, 15, 20));
                  detail = m -> "20 h asym: " * _ps_list(Dict(j => m[j] for j in ("a5_res16", "a5_res18")), "asym_20h") *
                                "; deflection: " * _ps_list(Dict(j => m[j] for j in ("a5_res16", "a5_res18")), "deflection_20h"; digits = 2) * "°"),
            Check("flow is downslope: speed-weighted mean direction within 5° of east, ≥ 90% of moving cells within ±30°",
                  m -> _ps_all(d -> abs(d["dir_mean_w"]) < 5 && d["dir_frac30"] >= 0.9, m);
                  detail = m -> "mean " * _ps_list(m, "dir_mean_w"; digits = 2) * "°; within ±30°: " * _ps_list(m, "dir_frac30"; digits = 3)),
        ],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
