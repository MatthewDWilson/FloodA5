# cases/registry/sgs_pool.jl — a pool at rest on rough ground: does a 25 m
# model store the water the 5 m ground holds? (docs/TEST_CASES.md §2e;
# HYDRAULICS §15.16.)
#
# Why: at the Carlisle peak SGS stores 0.70–0.79 of the 5 m SWE benchmark.
# Part of that could be the volume representation — a 25 m cell's table
# against the 5 m pixels it stands for — rather than the hydraulics. Water at
# rest separates the two: its level depends on the geometry alone, and on
# ground where every pixel drains to one sink the exact level for a volume V
# is the η with Σ_p a_p·max(η − z_p, 0) = V over the 5 m pixels.
#
# Terrain (cases/scripts/sgs_pool_dem.jl): 500 × 400 m, 5 m pixels, British
# National Grid; a bowl with a tilt, ridges and hollows of 25–80 m, pixel
# noise and a one-pixel channel, every depression filled to its spill level
# towards the sink. The 25 m grid is the exact 5 × 5 block means
# (cases/scripts/block_dem.jl), so its cells are whole blocks of the pixels.
# 52 866 m³ (the pixels below the 45th percentile of the ground) enters at the
# sink in 2 h; then 4 h at rest; closed boundaries.
#
# Runs: local-inertial on the 5 m pixels (the reference: it must reach the
# exact level); on the 25 m grid SGS with tables from the pixels
# (`--sgs-sampling pixels`: must reach it exactly), SGS with the default
# bilinear tables, and local-inertial on the block means. Per run, at the end
# of the inflow and at the end: the stored volume, the median level of the
# cells deeper than 5 cm and its spread, the exact level for the stored volume
# (`eta_exact`), and the level the run's own geometry gives for it at rest
# (`eta_self`: its tables, or its cell-mean beds) — so `eta_self − eta_exact`
# is the representation's error, and `eta − eta_self` whether the run is at rest.

const _SP_T = 21600.0
const _SP_T_FILL = 7200.0
const _SP_WET = 0.05
const _SP_LAT, _SP_LON = 54.909001, -2.956895638          # the sink pixel's centre (sgs_pool_dem.jl)
const _SP_JOBS = ("px5_li" => "local-inertial, 5 m pixels",
                  "c25_sgs_px" => "SGS 25 m, tables from the pixels",
                  "c25_sgs_bl" => "SGS 25 m, bilinear tables (default)",
                  "c25_li" => "local-inertial 25 m, block means")

_sp_params(dir) = FA.FO.JSON3.read(read(joinpath(dir, "sgs_pool_params.json"), String), Dict{String,Any})

"The level η at which `vol(η) = V` (bisection; `vol` increasing from 0 at `lo`)."
function _sp_level(vol, V, lo, hi)
    for _ in 1:200
        m = 0.5 * (lo + hi); vol(m) < V ? (lo = m) : (hi = m)
    end
    return 0.5 * (lo + hi)
end

function _sgs_pool_metrics(ctx)
    p = _sp_params(ctx.dir)
    out = Dict{String,Dict{String,Float64}}()
    out["design"] = Dict("V0_m3" => Float64(p["V0_m3"]), "eta0" => Float64(p["eta0"]),
                         "inflow_pixel_ok" => p["inflow_pixel_ok"] ? 1.0 : 0.0,
                         "n_pixels_raised" => Float64(p["n_pixels_raised"]))
    r5 = ctx.runs["px5_li"]
    zp = r5.mesh.elevation; ap = r5.mesh.area
    pix(η) = sum(ap[i] * max(η - zp[i], 0.0) for i in eachindex(zp))
    lo, hi = minimum(zp), maximum(zp) + 10.0
    for (name, r) in ctx.runs
        own = if r.flow_model == "sgs"
            t = FA.sgs_tables(r)
            η -> sum(FA.sgs_volume(ti, η) for ti in t)
        else
            z = r.mesh.elevation; a = r.mesh.area
            η -> sum(a[i] * max(η - z[i], 0.0) for i in eachindex(z))
        end
        d = Dict{String,Float64}("t_end" => r.t[end])
        for (tag, k) in (("fill", FA.frame_at(r, _SP_T_FILL; tol = 600.0)), ("end", FA.nframes(r)))
            V = FA.storage(r, k); η = FA.wse(r, k); h = FA.depth(r, k)
            w = sort([η[i] for i in eachindex(h) if isfinite(h[i]) && h[i] > _SP_WET])
            d["stored_$(tag)_m3"] = V
            d["eta_$tag"] = w[cld(length(w), 2)]
            d["eta_spread_$tag"] = w[end] - w[1]
            d["eta_exact_$tag"] = _sp_level(pix, V, lo, hi)
            d["eta_self_$tag"] = _sp_level(own, V, lo - 1.0, hi)
            d["eta_err_$tag"] = d["eta_$tag"] - d["eta_exact_$tag"]
            d["repr_err_$tag"] = d["eta_self_$tag"] - d["eta_exact_$tag"]
        end
        out[name] = d
    end
    return out
end

_sp(m, j, k) = get(get(m, j, Dict{String,Float64}()), k, NaN)
_sp_mm(x) = string(round(1000 * x; digits = 2), " mm")
_sp_detail(m) = join(("$j: level − exact $(_sp_mm(_sp(m, j, "eta_err_end"))), own geometry − exact " *
                      _sp_mm(_sp(m, j, "repr_err_end")) for j in first.(_SP_JOBS)), "; ")

let
    D = "{case}"
    dem5 = "$D/sgs_pool_dem5.tif"; dem25 = "$D/sgs_pool_dem25.tif"
    mesh(n) = "$D/$n.parquet"
    prep = [Prep("dem", :julia, ["{repo}/cases/scripts/sgs_pool_dem.jl", D];
                 outputs = [dem5, "$D/sgs_pool_params.json", "$D/inflow.csv"]),
            Prep("dem25", :julia, ["{repo}/cases/scripts/block_dem.jl", dem5, dem25, "5"]; outputs = [dem25]),
            Prep("m5", :floodmodel, ["--meshgen-cartesian", dem5, "--meshout", mesh("m5"),
                                     "--flow-model", "standard", "--mesh-only"]; outputs = [mesh("m5")]),
            Prep("m25_li", :floodmodel, ["--meshgen-cartesian", dem25, "--meshout", mesh("m25_li"),
                                         "--flow-model", "standard", "--mesh-only"]; outputs = [mesh("m25_li")]),
            Prep("m25_sgs_px", :floodmodel, ["--meshgen-cartesian", dem25, "--dem", dem5, "--meshout", mesh("m25_sgs_px"),
                                             "--flow-model", "sgs", "--sgs-sampling", "pixels", "--mesh-only"];
                 outputs = [mesh("m25_sgs_px")]),
            Prep("m25_sgs_bl", :floodmodel, ["--meshgen-cartesian", dem25, "--dem", dem5, "--meshout", mesh("m25_sgs_bl"),
                                             "--flow-model", "sgs", "--mesh-only"]; outputs = [mesh("m25_sgs_bl")])]
    run(m, fm) = ["--meshload", mesh(m), "--flow-model", fm,
                  "--inflow-point", "$(_SP_LAT),$(_SP_LON),$D/inflow.csv", "--closed-boundaries",
                  "--sim-duration", string(round(Int, _SP_T)), "--dt-max", "2", "--output-interval", "600"]
    jobs = [Job("px5_li", run("m5", "standard"); threads = 2, mesh = "cart5", solver = "li"),
            Job("c25_sgs_px", run("m25_sgs_px", "sgs"); threads = 2, mesh = "cart25", solver = "sgs", tables = "pixels"),
            Job("c25_sgs_bl", run("m25_sgs_bl", "sgs"); threads = 2, mesh = "cart25", solver = "sgs", tables = "bilinear"),
            Job("c25_li", run("m25_li", "standard"); threads = 2, mesh = "cart25", solver = "li")]
    register!(Case(
        name = "sgs_pool", tier = :standard,
        title = "A pool at rest on rough 5 m ground: the volume a 25 m model stores against the pixels' exact volume",
        doc = "docs/TEST_CASES.md §2e",
        description = """
        A closed, rough 500 × 400 m basin on a 5 m grid (British National Grid;
        a bowl, 25–80 m ridges and hollows, pixel noise, a one-pixel channel),
        every depression filled to its spill level towards one sink, so that
        water at rest stands at one level over exactly the pixels below it.
        52 866 m³ enters at the sink in 2 h; then 4 h at rest. The exact level
        for the stored volume V is the η with Σ a_p·max(η − z_p, 0) = V over
        the pixels. Runs: local-inertial on the 5 m pixels; on the 25 m grid of
        exact 5 × 5 block means, SGS with tables from the pixels
        (`--sgs-sampling pixels`), SGS with the default bilinear tables, and
        local-inertial on the block means. Measures each run's level against
        the exact level and the error of its geometry alone (the level its own
        tables or beds give for its volume at rest). About a minute per run.""",
        prep = prep, jobs = jobs,
        analyse = _sgs_pool_metrics,
        checks = [
            Check("the inflow point lies in the sink pixel", m -> _sp(m, "design", "inflow_pixel_ok") == 1.0),
            Check("mass conserved in every run (|error| < 1e-9 of the volume added)",
                  m -> all(_sp(m, j, "mass_balance_rel") < 1e-9 for j in first.(_SP_JOBS))),
            Check("every run stores the inflow volume (±0.1%)",
                  m -> all(abs(_sp(m, j, "stored_end_m3") / _sp(m, "design", "V0_m3") - 1) < 1e-3 for j in first.(_SP_JOBS))),
            Check("every run is at rest at the end: levels within 5 mm of each other, and within 1 mm of the " *
                  "level its own geometry gives for its volume",
                  m -> all(_sp(m, j, "eta_spread_end") < 0.005 &&
                           abs(_sp(m, j, "eta_end") - _sp(m, j, "eta_self_end")) < 0.001 for j in first.(_SP_JOBS));
                  detail = m -> join(("$j spread $(_sp_mm(_sp(m, j, "eta_spread_end")))" for j in first.(_SP_JOBS)), ", ")),
            Check("the 5 m run reaches the exact level (±1 mm)",
                  m -> abs(_sp(m, "px5_li", "eta_err_end")) < 0.001; detail = _sp_detail),
            Check("SGS with tables from the pixels: its tables hold exactly the pixels' volume (±0.1 mm) and " *
                  "it reaches the exact level (±1 mm)",
                  m -> abs(_sp(m, "c25_sgs_px", "repr_err_end")) < 1e-4 && abs(_sp(m, "c25_sgs_px", "eta_err_end")) < 0.001;
                  detail = _sp_detail),
            Check("SGS with the default (bilinear) tables is within 1 cm of the exact level",
                  m -> abs(_sp(m, "c25_sgs_bl", "eta_err_end")) < 0.01; detail = _sp_detail),
        ],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
