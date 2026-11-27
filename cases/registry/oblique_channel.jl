# cases/registry/oblique_channel.jl — a straight channel at an angle to the grid:
# does a model convey it the same whatever the angle? (docs/TEST_CASES.md §2f;
# HYDRAULICS §15.16.)
#
# Why: at Carlisle the 5 m runs (the SWE benchmark among them) overtop the
# narrow Petteril, which runs obliquely to the grid, where 25 m SGS keeps it in
# bank. On a Cartesian grid an oblique channel is a staircase of pixels: is
# its conveyance the channel's, or the staircase's?
#
# Terrain (cases/scripts/oblique_channel_dem.jl): a channel 15 m wide (three
# pixels) and 2 m deep falling at 0.002 between floodplains, on 5 m pixels,
# at 0°, 26.565° (1:2) and 45° to the grid; walls beyond the floodplains, and a
# steep drop past s = 1300 m so that no outlet formulation reaches the
# measuring reach s = 400–1000 m. 25 m³/s, Manning n = 0.035, 3 h (steady
# after about 1 h). Analytic: the 2-D normal depth for the rasterised channel
# (its pixels' area per metre of reach, W_r — 15.0, 15.71 and 17.62 m), as
# every model here sees those pixels.
#
# Runs per angle: local-inertial and SWE on the 5 m pixels, SGS on the exact
# 25 m blocks (`--sgs-sampling pixels`), local-inertial on the block means.
# Per run: the mean, least and greatest level − analytic level over the
# channel cells of the reach (5 m: the channel pixels; 25 m: the cells within
# 12.5 m of the axis), and the storage.

const _OC_ANGLES = ("0" => 0.0, "26.565" => 26.565, "45" => 45.0)
const _OC_MODELS = ("li5" => ("m5", "standard"), "swe5" => ("m5", "swe"),
                    "sgs25" => ("m25_sgs_px", "sgs"), "li25" => ("m25", "standard"))
const _OC_T = 10800.0
const _OC_REACH = (400.0, 1000.0)

"The generator's inflow point for an angle (its arithmetic: the job arguments need it first)."
function _oc_inflow(θ; P = 5.0, s_in = 20.0)
    mlon = 6_371_000.0 * cosd(51.0) * π / 180; mlat = 6_371_000.0 * π / 180
    x = P / 2 + s_in * cosd(θ); y = P / 2 - s_in * sind(θ)
    return (round(51.0 + y / mlat; digits = 9), round(x / mlon; digits = 9))
end

function _oc_metrics(ctx)
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in ctx.runs
        a = ctx.jobs[name].tags[:angle]
        p = FA.FO.JSON3.read(read(joinpath(ctx.dir, "a$a", "params.json"), String), Dict{String,Any})
        θ = Float64(p["angle_deg"]); d = (cosd(θ), -sind(θ)); nrm = (sind(θ), cosd(θ))
        x = (r.mesh.lon .- Float64(p["lon0"])) .* Float64(p["m_per_deg_lon"])
        y = (r.mesh.lat .- Float64(p["lat0"])) .* Float64(p["m_per_deg_lat"])
        s = x .* d[1] .+ y .* d[2]; n = x .* nrm[1] .+ y .* nrm[2]
        half = endswith(ctx.jobs[name].tags[:model], "25") ? 12.5 : Float64(p["W"]) / 2
        sel = [i for i in eachindex(s) if _OC_REACH[1] <= s[i] <= _OC_REACH[2] && abs(n[i]) < half]
        k = FA.frame_at(r, _OC_T; tol = 300.0); η = FA.wse(r, k); h = FA.depth(r, k)
        hn = Float64(p["h_normal"]); S = Float64(p["S"]); Z0 = Float64(p["Z0"])
        err = [η[i] - (Z0 - S * s[i] + hn) for i in sel if isfinite(η[i]) && h[i] > 0]
        lat_j, lon_j = _oc_inflow(θ)
        out[name] = Dict("angle_deg" => θ, "W_raster" => Float64(p["W_raster"]), "h_normal" => hn,
                         "err_mean" => isempty(err) ? NaN : sum(err) / length(err),
                         "err_min" => isempty(err) ? NaN : minimum(err), "err_max" => isempty(err) ? NaN : maximum(err),
                         "n_cells" => length(err), "stored_m3" => FA.storage(r, k),
                         "inflow_mismatch_deg" => max(abs(Float64(p["inflow_lat"]) - lat_j),
                                                      abs(Float64(p["inflow_lon"]) - lon_j)))
    end
    return out
end

_oc(m, j, k) = get(get(m, j, Dict{String,Float64}()), k, NaN)
_oc_jobs(model) = ["$(model)_a$a" for (a, _) in _OC_ANGLES]
_oc_detail(model) = m -> join(("$(a)°: $(round(_oc(m, "$(model)_a$a", "err_mean"); digits = 3)) m" for (a, _) in _OC_ANGLES), ", ")

let
    prep = Prep[]; jobs = Job[]
    for (a, θ) in _OC_ANGLES
        D = "{case}/a$a"
        dem5 = "$D/dem5.tif"; dem25 = "$D/dem25.tif"; mesh(n) = "$D/$n.parquet"
        append!(prep, [
            Prep("dem_a$a", :julia, ["{repo}/cases/scripts/oblique_channel_dem.jl", D, a];
                 outputs = [dem5, "$D/params.json", "$D/inflow.csv"]),
            Prep("dem25_a$a", :julia, ["{repo}/cases/scripts/block_dem.jl", dem5, dem25, "5"]; outputs = [dem25]),
            Prep("m5_a$a", :floodmodel, ["--meshgen-cartesian", dem5, "--meshout", mesh("m5"),
                                         "--flow-model", "standard", "--mesh-only"]; outputs = [mesh("m5")]),
            Prep("m25_a$a", :floodmodel, ["--meshgen-cartesian", dem25, "--meshout", mesh("m25"),
                                          "--flow-model", "standard", "--mesh-only"]; outputs = [mesh("m25")]),
            Prep("m25_sgs_px_a$a", :floodmodel, ["--meshgen-cartesian", dem25, "--dem", dem5, "--meshout", mesh("m25_sgs_px"),
                                                 "--flow-model", "sgs", "--sgs-sampling", "pixels", "--mesh-only"];
                 outputs = [mesh("m25_sgs_px")])])
        lat, lon = _oc_inflow(θ)
        for (model, (m, fm)) in _OC_MODELS
            push!(jobs, Job("$(model)_a$a",
                            ["--meshload", mesh(m), "--flow-model", fm, "--inflow-point", "$lat,$lon,$D/inflow.csv",
                             "--manning-n", "0.035", "--sim-duration", string(round(Int, _OC_T)), "--dt-max", "2",
                             "--output-interval", "600"];
                            threads = 2, angle = a, model = model))
        end
    end
    register!(Case(
        name = "oblique_channel", tier = :standard,
        title = "A straight channel at 0°, 26.6° and 45° to the grid: in-bank levels against the analytic normal depth",
        doc = "docs/TEST_CASES.md §2f",
        description = """
        A channel 15 m wide (three 5 m pixels) and 2 m deep, falling at 0.002
        between floodplains, at 0°, 26.565° and 45° to the grid; 25 m³/s in bank,
        n = 0.035, 3 h. Per angle: local-inertial and SWE on the 5 m pixels, SGS
        on the exact 25 m blocks (tables from the pixels), local-inertial on the
        block means. Measures the steady level over the reach s = 400–1000 m
        against the analytic 2-D normal depth for the rasterised channel. On a
        Cartesian grid an oblique channel is a staircase of pixels: the case
        shows how much each model's conveyance depends on the angle. About
        40 min (the 5 m SWE runs, 4–13 min each).""",
        prep = prep, jobs = jobs,
        analyse = _oc_metrics,
        checks = [
            Check("each run's inflow point is the generator's", m -> all(_oc(m, j, "inflow_mismatch_deg") < 1e-8 for j in keys(m))),
            Check("mass conserved in every run (|error| < 1e-6 of the volume added)",
                  m -> all(_oc(m, j, "mass_balance_rel") < 1e-6 for j in keys(m))),
            Check("aligned with the grid (0°), the 5 m runs and SGS are within 5 cm of the analytic level",
                  m -> all(abs(_oc(m, "$(x)_a0", "err_mean")) < 0.05 for x in ("li5", "swe5", "sgs25"));
                  detail = m -> join(("$x $(round(_oc(m, "$(x)_a0", "err_mean"); digits = 3))" for x in ("li5", "swe5", "sgs25")), ", ")),
            Check("SGS (25 m, tables from the pixels) is within 10 cm of the analytic level at every angle",
                  m -> all(abs(_oc(m, j, "err_mean")) < 0.10 for j in _oc_jobs("sgs25")); detail = _oc_detail("sgs25")),
            Check("local-inertial on the 5 m pixels is within 15 cm of the analytic level at every angle",
                  m -> all(abs(_oc(m, j, "err_mean")) < 0.15 for j in _oc_jobs("li5")); detail = _oc_detail("li5")),
            Check("SWE on the 5 m pixels is within 15 cm of the analytic level at every angle",
                  m -> all(abs(_oc(m, j, "err_mean")) < 0.15 for j in _oc_jobs("swe5")); detail = _oc_detail("swe5"),
                  xfail = "the staircase of an oblique channel's rasterised banks costs SWE head (form drag at " *
                          "every step): +0.85 m at 26.565° and +0.22 m at 45° here, +0.72 / +0.12 m on 2.5 m " *
                          "pixels (HYDRAULICS §15.16) — a property of the 5 m reference, not of SGS"),
        ],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
