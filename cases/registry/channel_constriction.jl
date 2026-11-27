# cases/registry/channel_constriction.jl — an abrupt contraction and expansion in a
# straight channel aligned with the grid: the head each model loses there
# (docs/TEST_CASES.md §2g; HYDRAULICS §15.16).
#
# Why: at the Carlisle peak the 5 m SWE run stands 0.4–0.5 m above 5 m
# local-inertial over the Eden's upstream floodplain, the step at the
# constricted reach through the city, on the original grid and on one turned
# 26.565°. Local-inertial has no convective acceleration, so it cannot lose
# head where the flow contracts and expands; the shallow-water equations
# conserve momentum, so they do. Here the channel is aligned with the grid (no
# staircase), and each answer is 1-D (cases/scripts/channel_constriction_1d.jl).
#
# Terrain (cases/scripts/channel_constriction_dem.jl): a channel 65 m wide
# between walls, bed slope 0.0005, narrowing abruptly to 35 m for 25 m at
# s = 2000 m; a broad-crested sill at s = 3500 m holds the channel above it at
# normal depth (2.99 m for 300 m³/s, n = 0.03). 4 h (steady after about 2 h).
#
# Runs: standard flow (local-inertial) and SWE on the 5 m pixels; SGS on the
# exact 25 m blocks (`--sgs-sampling pixels`); standard flow on the block means.
# Per run: the mean level over the channel cells of s = 1500–1950 m (above the
# opening) less each 1-D answer, started from the run's own level at
# s = 2600 m (so the downstream control drops out): `err_li` (local-inertial:
# dη/ds = −S_f), `err_energy` (lossless), `err_borda` (Borda–Carnot loss at the
# expansion), `err_abrupt` (abrupt-transition coefficients 0.6 / 0.8).

include(joinpath(@__DIR__, "..", "scripts", "channel_constriction_1d.jl"))

const _CC_T = 14400.0
const _CC_SD = 2600.0
const _CC_UP = (1500.0, 1950.0)
const _CC_JOBS = ("li5" => ("m5", "standard"), "swe5" => ("m5", "swe"),
                  "sgs25" => ("m25_sgs_px", "sgs"), "li25" => ("m25", "standard"))

"The generator's inflow point (its arithmetic, at 0°; the job arguments need it first)."
function _cc_inflow(; P = 5.0, s_in = 20.0)
    mlon = 6_371_000.0 * cosd(51.0) * π / 180; mlat = 6_371_000.0 * π / 180
    return (round(51.0 + (P / 2) / mlat; digits = 9), round((P / 2 + s_in) / mlon; digits = 9))
end

function _cc_metrics(ctx)
    p = FA.FO.JSON3.read(read(joinpath(ctx.dir, "params.json"), String), Dict{String,Any})
    θ = Float64(p["angle_deg"]); d = (cosd(θ), -sind(θ)); nrm = (sind(θ), cosd(θ))
    Q, n, S, Z0 = Float64(p["Q"]), Float64(p["n"]), Float64(p["S"]), Float64(p["Z0"])
    SC, LC = Float64(p["SC"]), Float64(p["LC"]); W1r, W2r = Float64(p["W1_raster"]), Float64(p["W2_raster"])
    W(s) = SC <= s < SC + LC ? W2r : W1r
    zb(s) = Z0 - S * s
    hn = (Q * n / (W1r * sqrt(S)))^(3 / 5)
    lat_j, lon_j = _cc_inflow()
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in ctx.runs
        x = (r.mesh.lon .- Float64(p["lon0"])) .* Float64(p["m_per_deg_lon"])
        y = (r.mesh.lat .- Float64(p["lat0"])) .* Float64(p["m_per_deg_lat"])
        s = x .* d[1] .+ y .* d[2]; nn = x .* nrm[1] .+ y .* nrm[2]
        k = FA.frame_at(r, _CC_T; tol = 300.0); η = FA.wse(r, k); h = FA.depth(r, k)
        inch(i) = isfinite(η[i]) && h[i] > 0 && abs(nn[i]) < Float64(p["W1"]) / 2
        dn = [η[i] for i in eachindex(s) if inch(i) && abs(s[i] - _CC_SD) <= 30]
        ηd = sum(dn) / length(dn)
        up = [i for i in eachindex(s) if inch(i) && _CC_UP[1] <= s[i] <= _CC_UP[2]]
        dd = Dict{String,Float64}("level_sd_minus_normal" => ηd - (zb(_CC_SD) + hn), "n_cells" => length(up),
                                  "level_up_mean" => sum(η[up]) / length(up), "stored_m3" => FA.storage(r, k),
                                  "inflow_mismatch_deg" => max(abs(Float64(p["inflow_lat"]) - lat_j),
                                                               abs(Float64(p["inflow_lon"]) - lon_j)))
        for m in (:li, :energy, :borda, :abrupt)
            dd["err_$m"] = try
                ss, ηs = cc_profile(m, ηd, _CC_SD, _CC_UP[1] - 100; Q, n, W, zb, edges = (SC + LC, SC))
                η1(sv) = (j = findfirst(<=(sv), ss); j == 1 ? ηs[1] : ηs[j] + (sv - ss[j]) / (ss[j-1] - ss[j]) * (ηs[j-1] - ηs[j]))
                sum(η[i] - η1(s[i]) for i in up) / length(up)
            catch
                NaN                                     # the 1-D answer chokes in the opening
            end
        end
        out[name] = dd
    end
    # the head the local-inertial runs do not lose: SWE's level above the opening less theirs
    if haskey(out, "swe5")
        for j in ("li5", "sgs25")
            haskey(out, j) && (out[j]["swe5_minus_this_up"] = out["swe5"]["level_up_mean"] - out[j]["level_up_mean"])
        end
    end
    return out
end

_ccm(m, j, k) = get(get(m, j, Dict{String,Float64}()), k, NaN)
_cc_r(x) = round(x; digits = 3)
_cc_detail(j) = m -> "$j − 1-D: local-inertial $(_cc_r(_ccm(m, j, "err_li"))), lossless $(_cc_r(_ccm(m, j, "err_energy"))), " *
                     "Borda $(_cc_r(_ccm(m, j, "err_borda"))), abrupt coefficients $(_cc_r(_ccm(m, j, "err_abrupt"))) m"

let
    D = "{case}"
    dem5 = "$D/dem5.tif"; dem25 = "$D/dem25.tif"; mesh(n) = "$D/$n.parquet"
    prep = [Prep("dem", :julia, ["{repo}/cases/scripts/channel_constriction_dem.jl", D];
                 outputs = [dem5, "$D/params.json", "$D/inflow.csv"]),
            Prep("dem25", :julia, ["{repo}/cases/scripts/block_dem.jl", dem5, dem25, "5"]; outputs = [dem25]),
            Prep("m5", :floodmodel, ["--meshgen-cartesian", dem5, "--meshout", mesh("m5"),
                                     "--flow-model", "standard", "--mesh-only"]; outputs = [mesh("m5")]),
            Prep("m25", :floodmodel, ["--meshgen-cartesian", dem25, "--meshout", mesh("m25"),
                                      "--flow-model", "standard", "--mesh-only"]; outputs = [mesh("m25")]),
            Prep("m25_sgs_px", :floodmodel, ["--meshgen-cartesian", dem25, "--dem", dem5, "--meshout", mesh("m25_sgs_px"),
                                             "--flow-model", "sgs", "--sgs-sampling", "pixels", "--mesh-only"];
                 outputs = [mesh("m25_sgs_px")])]
    lat, lon = _cc_inflow()
    jobs = [Job(name, ["--meshload", mesh(m), "--flow-model", fm, "--inflow-point", "$lat,$lon,$D/inflow.csv",
                       "--manning-n", "0.03", "--sim-duration", string(round(Int, _CC_T)), "--dt-max", "2",
                       "--output-interval", "900"]; threads = 2, model = name)
            for (name, (m, fm)) in _CC_JOBS]
    register!(Case(
        name = "channel_constriction", tier = :standard,
        title = "An abrupt contraction and expansion in a channel aligned with the grid: the head each model loses there",
        doc = "docs/TEST_CASES.md §2g",
        description = """
        A channel 65 m wide between walls (bed slope 0.0005) narrows abruptly to
        35 m for 25 m; 300 m³/s, n = 0.03, a sill downstream holding normal depth
        (2.99 m); aligned with the 5 m grid, so there is no staircase. Runs:
        standard flow and SWE on the 5 m pixels, SGS on the exact 25 m blocks,
        standard flow on the block means. The level above the opening is
        compared with 1-D answers started from each run's own level below it:
        local-inertial (no convective term: the level is continuous through the
        opening), lossless energy, energy with the Borda–Carnot expansion loss,
        and the step-backwater coefficients for abrupt transitions (0.6 / 0.8).
        Local-inertial models lose no head at the opening; the shallow-water
        equations do. About 20 min.""",
        prep = prep, jobs = jobs,
        analyse = _cc_metrics,
        checks = [
            Check("each run's inflow point is the generator's", m -> all(_ccm(m, j, "inflow_mismatch_deg") < 1e-8 for j in keys(m))),
            Check("mass conserved in every run (|error| < 1e-6 of the volume added)",
                  m -> all(_ccm(m, j, "mass_balance_rel") < 1e-6 for j in keys(m))),
            Check("standard flow on the 5 m pixels follows the local-inertial 1-D answer above the opening (±5 cm)",
                  m -> abs(_ccm(m, "li5", "err_li")) < 0.05; detail = _cc_detail("li5")),
            Check("SGS (25 m, tables from the pixels) follows the same local-inertial answer (±5 cm), as standard flow at 5 m does",
                  m -> abs(_ccm(m, "sgs25", "err_li")) < 0.05; detail = _cc_detail("sgs25")),
            Check("SWE on the 5 m pixels loses at least the Borda–Carnot head and no more than the abrupt-transition " *
                  "coefficients give (±5 cm)",
                  m -> _ccm(m, "swe5", "err_borda") > -0.05 && _ccm(m, "swe5", "err_abrupt") < 0.05; detail = _cc_detail("swe5")),
        ],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
