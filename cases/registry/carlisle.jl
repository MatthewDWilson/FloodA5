# cases/registry/carlisle.jl — Carlisle 2005, first 6 h, the 0.3.0 local-inertial
# default against the full shallow-water reference on the same mesh
# (docs/TEST_CASES.md §3; HYDRAULICS.md §15.10–15.11).

const _CARL = "{repo}/test/carlisle"
const _CARL_RUN = ["--inflow-bci", "$_CARL/carlisle_wgs84.bci", "--rainfall", "10",
                   "--sim-duration", "21600", "--dt-max", "10", "--output-interval", "900"]
# The 10 mm/h uniform rain is a test forcing (HYDRAULICS §15.6–15.11), not the 2005
# event's rain: 60 mm in 6 h over the whole domain. The *_norain runs show what the
# river inflow alone floods in the same 6 h (early in the event: 100–230 m³/s).
const _CARL_RUN_NORAIN = (i = findfirst(==("--rainfall"), _CARL_RUN); vcat(_CARL_RUN[1:i-1], _CARL_RUN[i+2:end]))
"""
Area (m²) wetter than `depth` at frame k. Cell-mean runs: whole cells with
h > depth. SGS runs (0.7.0): the area of each cell whose local depth exceeds
`depth`, from its tables (`FA.wetted_fraction`) — its depth is measured from
the cell's lowest ground, so counting whole cells by it overstates the extent
(at 6 h 4.96 km² against 1.79 km²).
"""
function _carl_wet_area(r, k, depth)
    if r.flow_model == "sgs"
        f = FA.wetted_fraction(r, k; wet = depth)
        return sum(r.mesh.area[i] * f[i] for i in eachindex(f) if isfinite(f[i]); init = 0.0)
    end
    h = FA.depth(r, k)
    return sum(r.mesh.area[i] for i in eachindex(h) if isfinite(h[i]) && h[i] > depth; init = 0.0)
end

"""
Time (s) at which the volume leaving the domain first reaches `vol` (m³),
interpolated linearly between the step lines of a run log (NaN if never).
"""
function _carl_arrival(logpath, vol)
    e = isfile(logpath) ? FA.parse_log_volumes(logpath) : []
    for i in eachindex(e)
        e[i].vol_removed >= vol || continue
        i == 1 && return e[i].t
        a, b = e[i-1], e[i]
        return a.t + (vol - a.vol_removed) / (b.vol_removed - a.vol_removed) * (b.t - a.t)
    end
    return NaN
end

"""
    _case_log_mass_balance(logpath) -> (added, removed, domain, error, rel)

The mass balance at a run log's last step line (`domain_vol`, `vol_removed`,
`mb_err`; added = domain + removed + error). FloodModel's end-of-run summary
does not fill its volumes (0.4.1: always 0), so the cases take them from the log.
"""
function _case_log_mass_balance(logpath)
    last = nothing
    if isfile(logpath)
        for l in eachline(logpath)
            m = match(r"domain_vol=([-0-9.e+]+)m³\s+vol_removed=([-0-9.e+]+)m³\s+mb_err=([-0-9.e+]+)m³", l)
            m === nothing || (last = m)
        end
    end
    last === nothing && return (added = NaN, removed = NaN, domain = NaN, error = NaN, rel = NaN)
    d, r, e = (parse(Float64, x) for x in last.captures)
    a = d + r + e
    return (added = a, removed = r, domain = d, error = e, rel = abs(e) / max(a, eps()))
end

"Storage, outflow and arrival of every run; each local-inertial run against the SWE run on its mesh."
function _carl_metrics(ctx)
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in ctx.runs
        k = FA.nframes(r)
        log = replace(r.path, ".h5" => ".log")
        mb = _case_log_mass_balance(log)
        # Σ volume: equal to Σ h·A except with SGS, whose depth is above the cell's lowest point
        out[name] = Dict("t" => r.t[k], "storage_m3" => FA.stored_volume(r, k),
                         "vol_added_m3" => mb.added, "vol_removed_m3" => mb.removed,
                         "mass_balance_log_rel" => mb.rel, "arrival_s" => _carl_arrival(log, 50_000.0),
                         "wet_area_0p1_km2" => _carl_wet_area(r, k, 0.1) / 1e6)
    end
    for (name, j) in ctx.jobs
        ref = get(j.tags, :reference, nothing)
        ref === nothing && continue
        c = FA.compare(ctx.runs[ref], ctx.runs[name])
        c.t == c.t_ref || error("$name vs $ref: frames at different times ($(c.t), $(c.t_ref))")
        merge!(out[name], Dict("storage_ratio" => c.storage_ratio, "rms_dh" => c.rms,
                               "csi" => c.csi, "bias" => c.bias, "speed_ratio_deep" => c.speed_ratio))
    end
    return out
end

_carl_fmt(x) = string(round(x; sigdigits = 4))
_carl_list(m, key, jobs) = join(("$j $(_carl_fmt(m[j][key]))" for j in jobs), ", ")

# Recorded 6 h values (HYDRAULICS §15.10–15.11; runs test/wp1_runs/r1b, bcnb, 0.2.6/0.3.0).
# A5 re-recorded in 0.8.0 (DEM sampled at pixel centres, §15.14): SWE 2.690 →
# 2.641 Mm³, arrival 6 688 → 6 676 s, SGS 1.755 → 1.775 Mm³. The Cartesian mesh
# takes the raster's pixels and is unchanged.
const _CARL_SWE_STORAGE = Dict("a5_swe" => 2.641e6, "cart_swe" => 2.331e6)
const _CARL_LI_RATIO = Dict("a5_li" => 1.25, "cart_li" => 0.99)

const _CARL_CHECKS = [
    # Of ~5.9 Mm³ added, SWE loses ~20 m³ (3e-6) and the local-inertial default ~150 m³
    # (2.5e-5, as in the recorded 0.3.0 runs)
    Check("mass conserved (|error| < 1e-4 of the volume added, from the log)",
          m -> all(m[j]["mass_balance_log_rel"] < 1e-4 for j in keys(m));
          detail = m -> _carl_list(m, "mass_balance_log_rel", sort(collect(keys(m))))),
    # The reference is converged in its own parameters (Courant 0.7 vs 0.35, eddy
    # viscosity: storage equal to 3 digits; twin rms 0.3 mm): 1% allows for
    # round-off-level changes of the bed (mesh regeneration) and the code since 0.2.6.
    Check("SWE 6 h storage = recorded (A5 2.641, Cartesian 2.331 Mm³) ± 1%",
          m -> all(abs(m[j]["storage_m3"] / v - 1) <= 0.01 for (j, v) in _CARL_SWE_STORAGE);
          detail = m -> _carl_list(m, "storage_m3", ("a5_swe", "cart_swe"))),
    # Halving the Courant number changes the local-inertial 6 h storage by 0.2%
    # (A5) and 4% (Cartesian) (§15.11): ±0.05 on the ratio covers that.
    Check("local-inertial / SWE storage = recorded (A5 1.25, Cartesian 0.99) ± 0.05",
          m -> all(abs(m[j]["storage_ratio"] - v) <= 0.05 for (j, v) in _CARL_LI_RATIO);
          detail = m -> _carl_list(m, "storage_ratio", ("a5_li", "cart_li"))),
    # §15.10: the A5 excess is the step-bed representation of sub-cell channels,
    # not the dissipation; it should stay well above the Cartesian ratio.
    Check("local-inertial excess larger on A5 than on Cartesian (ratio difference > 0.15)",
          m -> m["a5_li"]["storage_ratio"] - m["cart_li"]["storage_ratio"] > 0.15;
          detail = m -> _carl_list(m, "storage_ratio", ("a5_li", "cart_li"))),
    # Recorded: CSI 0.62 / 0.69, rms Δh 0.29 / 0.13 m (A5 / Cartesian, FA.compare
    # on the recorded runs).
    Check("extent vs SWE at 0.1 m: CSI ≥ 0.55; rms Δh ≤ 0.35 m",
          m -> all(m[j]["csi"] >= 0.55 && m[j]["rms_dh"] <= 0.35 for j in ("a5_li", "cart_li"));
          detail = m -> _carl_list(m, "csi", ("a5_li", "cart_li")) * "; rms " * _carl_list(m, "rms_dh", ("a5_li", "cart_li"))),
    # Recorded SWE arrival (50 000 m³ out): 6 676 s (A5, 0.8.0), 5 736 s (Cartesian).
    Check("SWE arrival (50 000 m³ out) = recorded (A5 6 676, Cartesian 5 736 s) ± 2%",
          m -> abs(m["a5_swe"]["arrival_s"] / 6676 - 1) <= 0.02 && abs(m["cart_swe"]["arrival_s"] / 5736 - 1) <= 0.02;
          detail = m -> _carl_list(m, "arrival_s", ("a5_swe", "cart_swe"))),
    # SGS has no reference here: a regression pin on the recorded 0.8.0 run,
    # ±2% as the local-inertial time-step sensitivity. (0.7.x: 1.755 Mm³; until
    # 0.6.0: 4.30 Mm³, 1.6× the SWE reference — conveyance held back;
    # cases/expected/v0.6.0/.)
    Check("SGS 6 h storage = recorded 1.775 Mm³ ± 2%",
          m -> abs(m["a5_sgs"]["storage_m3"] / 1.775e6 - 1) <= 0.02;
          detail = m -> _carl_list(m, "storage_m3", ("a5_sgs",))),
    # The 10 mm/h rain is a test forcing: early in the event (inflow 100–230 m³/s)
    # the river alone wets little more than its channels, so most of the 6 h
    # flooded area is rain-on-grid (A5: 3.06 km² with rain, 1.25 without).
    Check("most of the 6 h flooded area (> 0.1 m) is rain-on-grid: without rain it is < 60% of the area with rain",
          m -> all(m["$(g)_li_norain"]["wet_area_0p1_km2"] < 0.6 * m["$(g)_li"]["wet_area_0p1_km2"] for g in ("a5", "cart"));
          detail = m -> join(("$g: $(_carl_fmt(m["$(g)_li"]["wet_area_0p1_km2"])) km² with rain, " *
                              "$(_carl_fmt(m["$(g)_li_norain"]["wet_area_0p1_km2"])) without" for g in ("a5", "cart")), "; ")),
]

let
    a5mesh = "{case}/carlisle_a5_res18.parquet"
    cart = "{case}/carlisle_cart20.parquet"
    sgsmesh = "{case}/carlisle_a5_res18_sgs.parquet"
    register!(Case(
        name = "carlisle", tier = :full,
        title = "Carlisle 6 h: local-inertial default vs SWE reference (A5 res 18, Cartesian 20 m)",
        doc = "docs/TEST_CASES.md §3; HYDRAULICS.md §15.10–15.11",
        description = """
        The first 6 h of the Carlisle 2005 event (the three inflow hydrographs,
        open boundaries, the river initially dry) plus 10 mm/h uniform rain — a
        test forcing from the damping series (HYDRAULICS §15.6–15.11), not the
        event's rain: 60 mm over the whole domain, which floods most of the
        6 h extent (the two *_norain local-inertial runs show the river inflow
        alone) — on A5 res 18 (~22 m cells) and on a 20 m
        Cartesian grid aligned with the 5 m LiDAR averaged to 20 m, each with
        the 0.3.0 local-inertial default and with the full SWE reference
        (`--flow-model swe`, order 2), and SGS on A5 res 18 (no SWE comparison
        here; storage, and its flooded area from its tables). Metrics at 6 h: storage (Σ volume), and from the log the
        volume added and removed, the mass-balance error and the arrival time
        (50 000 m³ out); each local-inertial run against the SWE run on its
        mesh (storage ratio, rms Δh over wet cells, CSI and bias of the extent
        at 0.1 m, mean speed ratio where the reference is ≥ 3 m deep). The
        checks are the §15.10–15.11 record (A5 re-recorded in 0.8.0, §15.14):
        SWE 2.641 / 2.331 Mm³ and arrival 6 676 / 5 736 s; local-inertial
        1.25× SWE on A5, 0.99× on Cartesian; SGS 1.775 Mm³. CPU, one 4-thread run at a time: about 50 min (SWE 21
        and 17 min, the others 1.5–2.5 min), plus 7 min of mesh generation.
        A5 results are sensitive to round-off in the bed: a regenerated A5 mesh
        that differed from the recorded one by ≤ 3e-9 m in elevation moved the
        local-inertial storage by 2e-4 (rms Δh 1.5 cm). The expected values
        hold for the mesh they were recorded with; re-record after
        regenerating it on another machine.""",
        prep = [Prep("a5_res18", :floodmodel,
                     ["--meshgen", "$_CARL/Carlisle_domain.geojson", "--meshres", "18",
                      "--dem", "$_CARL/Carlisle_LiDAR_5m_mean.tif", "--meshout", a5mesh,
                      "--flow-model", "standard", "--mesh-only"]; outputs = [a5mesh]),
                Prep("dem20", :julia,
                     ["{repo}/cases/scripts/average_dem.jl", "$_CARL/Carlisle_LiDAR_5m_mean.tif",
                      "{case}/carlisle_dem20.tif", "20"]; outputs = ["{case}/carlisle_dem20.tif"]),
                Prep("cart20", :floodmodel,
                     ["--meshgen-cartesian", "{case}/carlisle_dem20.tif", "--meshout", cart, "--mesh-only"];
                     outputs = [cart]),
                Prep("a5_res18_sgs", :floodmodel,
                     ["--meshgen", "$_CARL/Carlisle_domain.geojson", "--meshres", "18",
                      "--dem", "$_CARL/Carlisle_LiDAR_5m_mean.tif", "--meshout", sgsmesh,
                      "--flow-model", "sgs", "--mesh-only"]; outputs = [sgsmesh])],
        jobs = [Job("a5_swe", vcat(["--meshload", a5mesh, "--flow-model", "swe"], _CARL_RUN);
                    mesh = "a5", solver = "swe"),
                Job("a5_li", vcat(["--meshload", a5mesh, "--flow-model", "standard"], _CARL_RUN);
                    mesh = "a5", solver = "li", reference = "a5_swe"),
                Job("cart_swe", vcat(["--meshload", cart, "--flow-model", "swe"], _CARL_RUN);
                    mesh = "cartesian", solver = "swe"),
                Job("cart_li", vcat(["--meshload", cart, "--flow-model", "standard"], _CARL_RUN);
                    mesh = "cartesian", solver = "li", reference = "cart_swe"),
                Job("a5_sgs", vcat(["--meshload", sgsmesh, "--flow-model", "sgs"], _CARL_RUN);
                    mesh = "a5", solver = "sgs"),
                Job("a5_li_norain", vcat(["--meshload", a5mesh, "--flow-model", "standard"], _CARL_RUN_NORAIN);
                    mesh = "a5", solver = "li", rain = false),
                Job("cart_li_norain", vcat(["--meshload", cart, "--flow-model", "standard"], _CARL_RUN_NORAIN);
                    mesh = "cartesian", solver = "li", rain = false)],
        analyse = _carl_metrics,
        checks = _CARL_CHECKS,
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end

# ── carlisle_peak: the 6 h before the River Eden's peak ─────────────────────
# The Eden ("main" in carlisle_wgs84.bdy) peaks at 1273 m³/s at t = 55.5 h; the
# comparison window is 49.5–55.5 h. The model has no start-time option, so the
# runs simulate the event from t = 0 (river initially dry; the first hours fill
# it). River inflow only: the 10 mm/h test rain of `carlisle` would add 555 mm
# over 55 h.
const _CARL_PEAK_T = 55.5 * 3600.0
const _CARL_PEAK_WINDOW = (49.5 * 3600.0, _CARL_PEAK_T)
const _CARL_PEAK_RUN = ["--inflow-bci", "$_CARL/carlisle_wgs84.bci",
                        "--sim-duration", string(Int(_CARL_PEAK_T)), "--dt-max", "10",
                        "--output-interval", "1800"]

# Output frames drift past the nominal times (the next output is due 1800 s after the
# last one was written, so each frame is late by up to one step, cumulatively): at
# 49.5 h the local-inertial frames are ~100 s late, the SWE ones ~30 s. Frames are
# taken within 300 s of the nominal time.
const _CARL_PEAK_TOL = 300.0

function _carl_peak_metrics(ctx)
    out = Dict{String,Dict{String,Float64}}()
    t0, t1 = _CARL_PEAK_WINDOW
    for (name, r) in ctx.runs
        d = Dict{String,Float64}("t_end" => r.t[end])
        for (tag, t) in (("w0", t0), ("wmid", (t0 + t1) / 2), ("w1", t1))
            k = FA.frame_at(r, t; tol = _CARL_PEAK_TOL)
            d["storage_$(tag)_m3"] = FA.stored_volume(r, k)
            d["wet_area_0p1_$(tag)_km2"] = _carl_wet_area(r, k, 0.1) / 1e6
        end
        out[name] = d
    end
    for (name, j) in ctx.jobs
        ref = get(j.tags, :reference, nothing)
        ref === nothing && continue
        for (tag, t) in (("w0", t0), ("wmid", (t0 + t1) / 2), ("w1", t1))
            kr = FA.frame_at(ctx.runs[ref], t; tol = _CARL_PEAK_TOL)
            if get(j.tags, :solver, "") == "sgs"
                merge!(out[name], _carl_sgs_compare(ctx.runs[ref], ctx.runs[name], kr, t, tag))
                continue
            end
            c = FA.compare(ctx.runs[ref], ctx.runs[name]; k_ref = kr)
            abs(c.t - c.t_ref) <= _CARL_PEAK_TOL || error("$name vs $ref: frames at different times")
            merge!(out[name], Dict("storage_ratio_$tag" => c.storage_ratio, "rms_dh_$tag" => c.rms,
                                   "csi_$tag" => c.csi, "bias_$tag" => c.bias))
        end
    end
    return out
end

"""
SGS against a reference on the same cells (0.7.0). An SGS depth is measured
from the cell's lowest ground, a cell-mean depth from its mean bed, so the two
are compared through quantities that mean the same in both: the stored volume,
the mean depth volume/area (wet: ≥ 0.1 m, for CSI and bias) and the water
level (rms over cells wet in both — for two runs on one bed this is the rms
depth difference, `rms_dh` of the local-inertial runs).
"""
function _carl_sgs_compare(ref, r, kr, t, tag)
    k = FA.frame_at(r, t; tol = _CARL_PEAK_TOL)
    abs(r.t[k] - ref.t[kr]) <= _CARL_PEAK_TOL || error("SGS vs reference: frames at different times")
    length(r.mesh.area) == length(ref.mesh.area) || error("SGS and its reference are on different meshes")
    A = r.mesh.area
    V = FA.volume(r, k); Vr = FA.volume(ref, kr); η = FA.wse(r, k); ηr = FA.wse(ref, kr)
    ok(i) = isfinite(V[i]) && isfinite(Vr[i]) && A[i] > 0
    wet = [ok(i) && V[i] / A[i] >= 0.1 for i in eachindex(A)]
    wetr = [ok(i) && Vr[i] / A[i] >= 0.1 for i in eachindex(A)]
    both = wet .& wetr
    hit = sum(A[both]); un = sum(A[wet .| wetr])
    dη = [η[i] - ηr[i] for i in eachindex(A) if both[i] && isfinite(η[i]) && isfinite(ηr[i])]
    return Dict("storage_ratio_$tag" => FA.stored_volume(r, k) / FA.stored_volume(ref, kr),
                "rms_dh_$tag" => isempty(dη) ? NaN : sqrt(sum(abs2, dη) / length(dη)),
                "mean_dh_$tag" => isempty(dη) ? NaN : sum(dη) / length(dη),
                "csi_$tag" => un > 0 ? hit / un : NaN, "bias_$tag" => sum(A[wet]) / max(sum(A[wetr]), eps()))
end

let
    a5mesh = "{case}/../carlisle/carlisle_a5_res18.parquet"      # the `carlisle` case's meshes
    cart = "{case}/../carlisle/carlisle_cart20.parquet"
    sgsmesh = "{case}/../carlisle/carlisle_a5_res18_sgs.parquet"
    register!(Case(
        name = "carlisle_peak", tier = :full,
        title = "Carlisle, the 6 h before the Eden peak (49.5–55.5 h): local-inertial vs SWE, river inflow only",
        doc = "docs/TEST_CASES.md §3",
        description = """
        The Carlisle 2005 event from t = 0 to the River Eden's peak (1273 m³/s at
        55.5 h), with the three inflow hydrographs only (no rain), open
        boundaries and the river initially dry. Compared over the 6 h before the
        peak (49.5, 52.5 and 55.5 h): each local-inertial run against the SWE
        reference on its mesh (storage ratio, rms Δh, CSI and bias of the extent
        at 0.1 m), on A5 res 18 and Cartesian 20 m; SGS on A5 (with and without
        the critical-flow limit of 0.8.0) against the A5 SWE run (storage, water
        level where both are wet, extent by mean depth V/A), for information:
        the 22 m SWE run is not a reference. The benchmark is the 5 m SWE run
        of `carlisle_peak_fine`, whose checks hold the agreed SGS bar.
        Uses the `carlisle` case's meshes (run `cases/run.jl prep carlisle`
        first).""",
        prep = Prep[],
        jobs = [Job("a5_swe", vcat(["--meshload", a5mesh, "--flow-model", "swe"], _CARL_PEAK_RUN);
                    threads = 6, mesh = "a5", solver = "swe"),
                Job("a5_li", vcat(["--meshload", a5mesh, "--flow-model", "standard"], _CARL_PEAK_RUN);
                    mesh = "a5", solver = "li", reference = "a5_swe"),
                Job("cart_swe", vcat(["--meshload", cart, "--flow-model", "swe"], _CARL_PEAK_RUN);
                    threads = 6, mesh = "cartesian", solver = "swe"),
                Job("cart_li", vcat(["--meshload", cart, "--flow-model", "standard"], _CARL_PEAK_RUN);
                    mesh = "cartesian", solver = "li", reference = "cart_swe"),
                Job("a5_sgs", vcat(["--meshload", sgsmesh, "--flow-model", "sgs"], _CARL_PEAK_RUN);
                    mesh = "a5", solver = "sgs", reference = "a5_swe"),
                # the 0.7.x SGS hydraulics (no critical-flow limit) on the same tables
                Job("a5_sgs_nocrit", vcat(["--meshload", sgsmesh, "--flow-model", "sgs", "--sgs-critical", "off"],
                                          _CARL_PEAK_RUN);
                    mesh = "a5", solver = "sgs", reference = "a5_swe")],
        analyse = _carl_peak_metrics,
        checks = [Check("mass conserved (|conservation error| < 1e-6 of the volume applied)",
                        m -> all(d["mass_balance_rel"] < 1e-6 for (_, d) in m);
                        detail = m -> _carl_list(m, "mass_balance_rel", sort(collect(keys(m))))),
                  # the inputs are the hydrograph integral (157.63 Mm³); the step-start
                  # evaluation of the rates is first order (FloodModel.jl run_simulation!)
                  Check("volume added = the hydrographs' integral (|input error| < 1e-5)",
                        m -> all(abs(d["input_error_rel"]) < 1e-5 for (_, d) in m);
                        detail = m -> _carl_list(m, "input_error_rel", sort(collect(keys(m))))),
                  # At high flow the local-inertial default conveys more than the full
                  # SWE (faster in deep water, lets more out) — the opposite sign to the
                  # 6 h rain test (`carlisle`: 1.25× on A5). Recorded 0.74–0.76 (A5),
                  # 0.85–0.86 (Cartesian).
                  Check("near the peak local-inertial stores less than SWE on both meshes (ratio < 1 at 49.5, 52.5, 55.5 h)",
                        m -> all(m[j]["storage_ratio_$w"] < 1 for j in ("a5_li", "cart_li"), w in ("w0", "wmid", "w1"));
                        detail = m -> _carl_list(m, "storage_ratio_w1", ("a5_li", "cart_li")) * " at 55.5 h"),
                  Check("extent vs SWE at 0.1 m through the window: CSI ≥ 0.80, bias in [0.80, 1.00]; rms Δh ≤ 0.65 m (A5), ≤ 0.45 m (Cartesian)",
                        m -> all(m[j]["csi_$w"] >= 0.80 && 0.80 <= m[j]["bias_$w"] <= 1.0 for j in ("a5_li", "cart_li"), w in ("w0", "wmid", "w1")) &&
                             all(m["a5_li"]["rms_dh_$w"] <= 0.65 && m["cart_li"]["rms_dh_$w"] <= 0.45 for w in ("w0", "wmid", "w1"));
                        detail = m -> "CSI " * _carl_list(m, "csi_w1", ("a5_li", "cart_li")) * "; rms " *
                                      _carl_list(m, "rms_dh_w1", ("a5_li", "cart_li")) * " at 55.5 h"),
                  # the two LI runs nearly coincide (1%); the references differ by 12%
                  Check("the A5 and Cartesian local-inertial runs store within 5% of each other at 55.5 h",
                        m -> abs(m["a5_li"]["storage_w1_m3"] / m["cart_li"]["storage_w1_m3"] - 1) < 0.05;
                        detail = m -> _carl_list(m, "storage_w1_m3", ("a5_li", "cart_li")))],
        # The SGS bar (agreed 2026-11-16) was checked here against the 22 m SWE run
        # until 0.7.1 (an expected failure); since 0.8.0 it is checked against the
        # 5 m SWE benchmark in `carlisle_peak_fine`. Until 0.7.0 SGS stored 2.5×
        # SWE here, the Eden ~3 m high: cases/expected/v0.6.0/.
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end

# ── carlisle_peak_fine: 5 m references for the peak (0.7.0) ────────────────
# The 22 m runs of `carlisle_peak` cannot say which of them is nearer the
# truth: their reference is SWE on cell-mean beds at the same resolution, and
# an SGS cell with water in a sub-grid channel is not the same thing as a
# cell-mean cell with that volume spread over it. A 5 m Cartesian grid
# resolves the banks: local-inertial and SWE on it are what the 22 m models
# (standard, SWE and SGS) should each be judged against (plan of 2026-11,
# §9.3). GPU runs: 582 012 cells for 55.5 h is hours for local-inertial and
# about a day for SWE; they are not part of `carlisle_peak`, which must stay
# runnable on a CPU in under an hour. The mesh is the benchmark suite's
# (`julia --project=. benchmarks/run.jl`, test/case_runs/perf/).
"The 22 m runs of `carlisle_peak` compared with the 5 m runs (`_carl_fine_metrics`)."
const _CARL_FINE_VS = ("a5_swe", "a5_li", "a5_sgs", "a5_sgs_nocrit", "cart_swe", "cart_li")

"""
Metrics of `carlisle_peak_fine`: each 5 m run's storage, and for each 22 m run
of `carlisle_peak` (entries `vs5:<job>`) its comparison with each 5 m run
aggregated onto its cells (`FA.compare_to_fine`; prefixes `swe5_`, `li5_`) at
49.5, 52.5 and 55.5 h (`_w0`, `_wmid`, `_w1`): storage ratio, the median and
rms of the water-level difference where both are at least half wet, and the
CSI of the wetted area at 0.1 m local depth (SGS: from its tables).
"""
function _carl_fine_metrics(ctx)
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in ctx.runs
        out[name] = Dict("t_end" => r.t[end], "storage_w1_m3" => FA.stored_volume(r, FA.nframes(r)))
    end
    peak = case("carlisle_peak")
    t0, t1 = _CARL_PEAK_WINDOW
    owners = Dict{String,Any}()                    # coarse mesh → fine-cell owners
    for j in peak.jobs
        j.name in _CARL_FINE_VS || continue
        job_status(peak, j) === :done || (@warn "carlisle_peak_fine: carlisle_peak/$(j.name) is not complete"; continue)
        r = FA.load_run(job_output(peak, j))
        mesh = j.args[findfirst(==("--meshload"), j.args) + 1]
        tables = get(j.tags, :solver, "") == "sgs" ? FA.sgs_tables(r) : nothing
        d = Dict{String,Float64}()
        for (fname, fjob) in (("swe5", "cart5_swe"), ("li5", "cart5_li"))
            for (tag, t) in (("w0", t0), ("wmid", (t0 + t1) / 2), ("w1", t1))
                c = FA.compare_to_fine(ctx.runs[fjob], r, t; tol = _CARL_PEAK_TOL,
                                       owner = get(owners, mesh, nothing), tables)
                owners[mesh] = c.owner
                d["$(fname)_storage_ratio_$tag"] = c.storage_ratio
                d["$(fname)_level_median_$tag"] = c.level_median
                d["$(fname)_level_rms_$tag"] = c.level_rms
                d["$(fname)_csi_$tag"] = c.csi
            end
        end
        out["vs5:" * j.name] = d
    end
    return out
end

"SGS (`s`) no further from the 5 m run `ref` (SWE) than local-inertial (`l`) at the three times."
_carl_bar(m, s, l; ref = "swe5") =
    all(abs(m[s]["$(ref)_storage_ratio_$w"] - 1) <= abs(m[l]["$(ref)_storage_ratio_$w"] - 1) &&
        m[s]["$(ref)_level_rms_$w"] <= m[l]["$(ref)_level_rms_$w"] &&
        m[s]["$(ref)_csi_$w"] >= m[l]["$(ref)_csi_$w"] for w in ("w0", "wmid", "w1"))
function _carl_bar_detail(m, s, l; ref = "swe5")
    r(k, w) = round(m[s]["$(ref)_$(k)_$w"]; digits = 3)
    rl(k, w) = round(m[l]["$(ref)_$(k)_$w"]; digits = 3)
    join(("$w: storage ÷ 5 m $(ref == "swe5" ? "SWE" : "LI") $(r("storage_ratio", w)) " *
          "(LI $(rl("storage_ratio", w))), rms Δη $(r("level_rms", w)) m (LI $(rl("level_rms", w))), " *
          "CSI $(r("csi", w)) (LI $(rl("csi", w)))" for w in ("w0", "wmid", "w1")), "; ")
end

let
    cart5 = "{case}/../perf/carlisle_cart5.parquet"
    gpu(prec) = ["--backend", "gpu", "--precision", prec]
    register!(Case(
        name = "carlisle_peak_fine", tier = :full,
        title = "Carlisle to the Eden peak on a 5 m Cartesian grid (GPU): the fine reference for carlisle_peak",
        doc = "docs/TEST_CASES.md §3",
        description = """
        The `carlisle_peak` event (55.5 h, river inflow only) on the 5 m Cartesian
        mesh of the same domain (582 012 cells), local-inertial (GPU, Float64) and
        SWE (GPU, Float32), output every 30 min. The 5 m SWE run is the
        benchmark for the 22 m runs of `carlisle_peak` — standard, SWE and SGS
        (decided 2026-11-23) — which this case's analysis compares with both
        5 m runs on their own cells (run `carlisle_peak` first); it holds the
        agreed SGS bar. The 5 m mesh takes its elevations from the raster's
        pixels, so the DEM registration fix of 0.8.0 does not change it.
        Needs a CUDA GPU and the 5 m mesh of the benchmark suite; hours
        (local-inertial) to about a day (SWE).""",
        prep = Prep[],
        jobs = [Job("cart5_li", vcat(["--meshload", cart5, "--flow-model", "standard"], _CARL_PEAK_RUN, gpu("float64"));
                    threads = 2, mesh = "cartesian5", solver = "li"),
                Job("cart5_swe", vcat(["--meshload", cart5, "--flow-model", "swe"], _CARL_PEAK_RUN, gpu("float32"));
                    threads = 2, mesh = "cartesian5", solver = "swe")],
        analyse = _carl_fine_metrics,
        checks = [Check("mass conserved (|conservation error| < 1e-6 of the volume applied)",
                        m -> all(d["mass_balance_rel"] < 1e-6 for (k, d) in m if haskey(d, "mass_balance_rel"))),
                  Check("the 22 m runs of carlisle_peak are all compared (run that case first)",
                        m -> all(haskey(m, "vs5:" * j) for j in _CARL_FINE_VS)),
                  # The bar for SGS (plan of 2026-11, agreed 2026-11-16; benchmark the
                  # 5 m SWE run since 2026-11-23): if a sub-grid model is no closer to
                  # the full shallow-water equations resolved at 5 m than the cell-mean
                  # local-inertial model on the same mesh, there is little evidence of
                  # its benefit.
                  Check("SGS on A5 is no further from the 5 m SWE benchmark than local-inertial on the same mesh, " *
                        "at 49.5, 52.5 and 55.5 h: storage, rms water-level difference and extent (CSI at 0.1 m)",
                        m -> _carl_bar(m, "vs5:a5_sgs", "vs5:a5_li");
                        detail = m -> _carl_bar_detail(m, "vs5:a5_sgs", "vs5:a5_li") *
                                      " | against the 5 m local-inertial run: " *
                                      _carl_bar_detail(m, "vs5:a5_sgs", "vs5:a5_li"; ref = "li5"),
                        xfail = "the hydraulics, not the SGS geometry (HYDRAULICS §15.14): at 5 m the " *
                                "local-inertial equations store 0.84 of SWE's volume, and SGS — local-inertial " *
                                "on its tables — stores 0.87 of the 5 m local-inertial run, closer to it in " *
                                "level and extent than the cell-mean run; the cell-mean run matches SWE's " *
                                "storage by compensating errors (its beds store ~1.15× the 5 m local-inertial " *
                                "run). Needs a hydraulic formulation with momentum on the SGS geometry"),
                  # The same bar against the 5 m standard-flow (local-inertial) run (added
                  # 0.9.1): the like-for-like reference for SGS as it is — local-inertial
                  # hydraulics on its tables — free of the head the shallow-water equations
                  # lose at the Eden's constricted reach (case channel_constriction).
                  Check("SGS on A5 is no further from the 5 m standard-flow run than local-inertial (cell-mean) on " *
                        "the same mesh, at 49.5, 52.5 and 55.5 h: storage, rms water-level difference and extent (CSI at 0.1 m)",
                        m -> _carl_bar(m, "vs5:a5_sgs", "vs5:a5_li"; ref = "li5");
                        detail = m -> _carl_bar_detail(m, "vs5:a5_sgs", "vs5:a5_li"; ref = "li5"),
                        xfail = "storage only, at 55.5 h, by 0.01 (|0.872 − 1| against |1.118 − 1|; SGS is closer in " *
                                "level and extent at every time): the 5 m run's Petteril overtops into the city and " *
                                "the defended compartments, a threshold that moves with the grid's alignment — on a " *
                                "5 m grid turned 26.565° standard flow stores 6.6% less (HYDRAULICS §15.16; Nguyen et " *
                                "al. 2025, doi:10.1029/2024WR038919)")],
    ))
end
