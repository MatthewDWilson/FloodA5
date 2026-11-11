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
_carl_wet_area(r, k, depth) = (h = FA.depth(r, k);
    sum(r.mesh.area[i] for i in eachindex(h) if isfinite(h[i]) && h[i] > depth; init = 0.0))

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

# Recorded 6 h values (HYDRAULICS §15.10–15.11; runs test/wp1_runs/r1b, bcnb, 0.2.6/0.3.0)
const _CARL_SWE_STORAGE = Dict("a5_swe" => 2.690e6, "cart_swe" => 2.331e6)
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
    Check("SWE 6 h storage = recorded (A5 2.690, Cartesian 2.331 Mm³) ± 1%",
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
    # Recorded SWE arrival (50 000 m³ out): 6 688 s (A5), 5 736 s (Cartesian).
    Check("SWE arrival (50 000 m³ out) = recorded (A5 6 688, Cartesian 5 736 s) ± 2%",
          m -> abs(m["a5_swe"]["arrival_s"] / 6688 - 1) <= 0.02 && abs(m["cart_swe"]["arrival_s"] / 5736 - 1) <= 0.02;
          detail = m -> _carl_list(m, "arrival_s", ("a5_swe", "cart_swe"))),
    # SGS has no reference: the recorded 0.3.0 run (§15.11, 4.299 Mm³; Σ volume
    # at 21 600 s 4.300). ±2%, as the local-inertial time-step sensitivity.
    Check("SGS 6 h storage = recorded 4.30 Mm³ ± 2%",
          m -> abs(m["a5_sgs"]["storage_m3"] / 4.30e6 - 1) <= 0.02;
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
        (`--flow-model swe`, order 2), and SGS on A5 res 18 (no reference;
        storage only). Metrics at 6 h: storage (Σ volume), and from the log the
        volume added and removed, the mass-balance error and the arrival time
        (50 000 m³ out); each local-inertial run against the SWE run on its
        mesh (storage ratio, rms Δh over wet cells, CSI and bias of the extent
        at 0.1 m, mean speed ratio where the reference is ≥ 3 m deep). The
        checks are the §15.10–15.11 record: SWE 2.690 / 2.331 Mm³ and arrival
        6 688 / 5 736 s; local-inertial 1.25× SWE on A5, 0.99× on Cartesian;
        SGS 4.30 Mm³. CPU, one 4-thread run at a time: about 50 min (SWE 21
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
            c = FA.compare(ctx.runs[ref], ctx.runs[name]; k_ref = FA.frame_at(ctx.runs[ref], t; tol = _CARL_PEAK_TOL))
            abs(c.t - c.t_ref) <= _CARL_PEAK_TOL || error("$name vs $ref: frames at different times")
            merge!(out[name], Dict("storage_ratio_$tag" => c.storage_ratio, "rms_dh_$tag" => c.rms,
                                   "csi_$tag" => c.csi, "bias_$tag" => c.bias))
        end
    end
    return out
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
        at 0.1 m), on A5 res 18 and Cartesian 20 m; SGS on A5 (no reference).
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
                    mesh = "a5", solver = "sgs")],
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
        # no SGS check: its 6 h-before-peak storage (2.5× SWE; Eden ~3 m high) is
        # under investigation (PROJECT_STATE) and a check would enshrine it
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
