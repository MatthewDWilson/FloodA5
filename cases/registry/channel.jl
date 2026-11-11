# cases/registry/channel.jl — the straight-channel representation series
# (docs/TEST_CASES.md §7; HYDRAULICS.md §15.12): a 20 m channel at several
# orientations and cross-channel offsets on A5 and equal-area Cartesian
# meshes, full SWE and the local-inertial default, 3 h to a steady state.
#
# The meshes, inflow file and manifest come from
# test/channel_bench/make_channel_meshes.jl (a preparation step). The jobs
# are defined here, statically, with the same names and arguments as the
# generator's jobs.txt; the metric `jobs_match_generator` (under "series")
# checks that they agree with the jobs.txt the preparation wrote.

# The generator's constants (make_channel_meshes.jl; FA.CHANNEL)
const _CH_L, _CH_LAT0, _CH_LON0, _CH_R = 2000.0, 51.0, 0.0, 6_371_000.0   # R = A5Grid._EARTH_R
"Nominal cell size (m) of each A5 resolution over the channel domain: √(mean cell area), as the generator computes it."
const _CH_CELL = Dict(17 => 44.40893177272335, 18 => 22.204465875421874,
                      19 => 11.102232940083269, 20 => 5.551116469880078)
const _CH_RUN = ["--sim-duration", "10800", "--dt-max", "10", "--output-interval", "900", "--manning-n", "0.03"]
const _CH_SOLVER = Dict("swe" => ["--flow-model", "swe"], "li" => ["--flow-model", "standard"])

"The generator's case name, e.g. `a5_r18_a45_o00`."
_ch_case(mesh, res, angle, of) = string(mesh, "_r", res, "_a", lpad(round(Int, angle), 2, '0'),
                                        "_o", lpad(round(Int, 100of), 2, '0'))

"The inflow point (lat, lon), with the generator's arithmetic and rounding."
function _ch_inflow(angle, off)
    dlon(x) = rad2deg(x / (_CH_R * cosd(_CH_LAT0))); dlat(y) = rad2deg(y / _CH_R)
    s0 = -0.45_CH_L
    x = s0 * cosd(angle) - off * sind(angle); y = s0 * sind(angle) + off * cosd(angle)
    return round(_CH_LAT0 + dlat(y), digits = 9), round(_CH_LON0 + dlon(x), digits = 9)
end

"Jobs `chan_<case>_<solver>` for every resolution × angle × offset (fraction of a cell) × mesh × solver."
function _channel_jobs(; res, angles, offsets, meshes = ("a5", "cart"), solvers = ("swe", "li"))
    jobs = Job[]
    for r in res, mesh in meshes, α in angles, of in offsets
        c = _ch_case(mesh, r, α, of)
        off = of * _CH_CELL[r]
        lat, lon = _ch_inflow(Float64(α), off)
        src = ["--meshload", "{case}/$c.parquet", "--inflow-point", "$lat,$lon,{case}/q20.csv"]
        for s in solvers
            push!(jobs, Job("chan_$(c)_$s", vcat(src, _CH_RUN, _CH_SOLVER[s]);
                            mesh, res = r, angle = Float64(α), offset = of, solver = s))
        end
    end
    return jobs
end

"The generator run that prepares `res × angles × offsets` (both meshes) in `{case}`."
function _channel_prep(; res, angles, offsets)
    cases = [_ch_case(m, r, α, of) for r in res, m in ("a5", "cart"), α in angles, of in offsets]
    Prep("meshes", :julia,
         ["{repo}/test/channel_bench/make_channel_meshes.jl", "{case}",
          "--res", join(res, ","), "--angles", join(angles, ","), "--offsets", join(offsets, ",")];
         outputs = vcat(["{case}/manifest.csv", "{case}/jobs.txt", "{case}/q20.csv"],
                        ["{case}/$c.parquet" for c in vec(cases)]))
end

_ch_norm(s) = replace(s, '\\' => '/')
"1 if every job's arguments are those of the generator's jobs.txt (paths compared with / separators), else 0."
function _ch_jobs_match(ctx)
    gen = Dict{String,String}()
    for l in eachline(joinpath(ctx.dir, "jobs.txt"))
        n, a = split(l, '|'; limit = 2)
        gen[n] = _ch_norm(a)
    end
    bad = String[]
    for (name, j) in ctx.jobs
        a = _ch_norm(join(_subst(j.args, ctx.case), " "))
        get(gen, name, "") == a || push!(bad, name)
    end
    isempty(bad) || @warn "jobs that differ from the generator's jobs.txt: $(join(sort(bad), ", "))"
    return isempty(bad) ? 1.0 : 0.0
end

# (`_case_log_mass_balance` is in carlisle.jl, the same module.)
"Per-run channel metrics (FA.channel_series_metrics over the case's manifest) and series-level summaries."
function _channel_metrics(ctx)
    rows = FA.channel_series_metrics(ctx.dir, joinpath(ctx.dir, "runs"); solvers = ("swe", "li"),
                                     csv = joinpath(ctx.dir, "channel_metrics.csv"), io = nothing)
    out = Dict{String,Dict{String,Float64}}()
    for r in rows
        out["chan_$(r.case)_$(r.solver)"] = Dict("deta" => r.deta, "spill" => r.spill,
                                                 "unsteady" => r.unsteady, "outflow" => r.outflow,
                                                 "mass_balance_log_rel" => _case_log_mass_balance(joinpath(ctx.dir, "runs", "chan_$(r.case)_$(r.solver).log")).rel)
    end
    out["series"] = merge(_ch_series(rows), Dict("jobs_match_generator" => _ch_jobs_match(ctx)))
    return out
end

"""
HYDRAULICS §15.12's table from metric rows (`solver`, `mesh`, `res`, `angle`,
`deta`): per solver × mesh × resolution, the rms and spread (max − min) of Δη
over the orientations, each orientation's Δη the mean over the offsets.
"""
function _ch_series(rows)
    series = Dict{String,Float64}()
    for s in unique(r.solver for r in rows), m in unique(r.mesh for r in rows), res in unique(r.res for r in rows)
        sel = [r for r in rows if r.solver == s && r.mesh == m && r.res == res]
        isempty(sel) && continue
        per_angle = [sum(r.deta for r in sel if r.angle == a) / count(r -> r.angle == a, sel)
                     for a in unique(r.angle for r in sel)]
        series["$(s)_$(m)_r$(res)_rms"] = sqrt(sum(abs2, per_angle) / length(per_angle))
        series["$(s)_$(m)_r$(res)_spread"] = maximum(per_angle) - minimum(per_angle)
    end
    return series
end

# ── checks ──
_ch_runs(m) = [(j, d) for (j, d) in m if startswith(j, "chan_")]
_ch_job(m, mesh, res, α, of, s) = m["chan_$(_ch_case(mesh, res, α, of))_$s"]
_ch_list(m, key; filt = j -> true) = join(("$(replace(j, "chan_" => "")) $(round(d[key]; digits = 3))"
                                           for (j, d) in sort(_ch_runs(m); by = first) if filt(j)), ", ")
_ch_res(j) = parse(Int, match(r"_r(\d+)_", j)[1])

function _channel_checks(; res, angles, offsets)
    fine, coarse = maximum(res), minimum(res)
    pairs = [(mesh, α, of, s) for mesh in ("a5", "cart"), α in angles, of in offsets, s in ("swe", "li")]
    [Check("mass conserved (|error| < 1e-5 of the volume added, from the log)",
           m -> all(d["mass_balance_log_rel"] < 1e-5 for (_, d) in _ch_runs(m));
           detail = m -> "worst " * string(maximum(d["mass_balance_log_rel"] for (_, d) in _ch_runs(m)))),
     Check("job arguments = the generator's jobs.txt", m -> m["series"]["jobs_match_generator"] == 1),
     Check("steady at 3 h: outflow/inflow within 2% at $(round(_CH_CELL[fine]; digits = 1)) m, 15% at coarser cells",
           m -> all(abs(d["outflow"] - 1) <= (_ch_res(j) == fine ? 0.02 : 0.15) for (j, d) in _ch_runs(m));
           detail = m -> _ch_list(m, "outflow")),
     Check("Δη falls by ≥ 0.2 m from $(round(_CH_CELL[coarse]; digits = 1)) to $(round(_CH_CELL[fine]; digits = 1)) m (every mesh, angle, offset, solver)",
           m -> all(_ch_job(m, p[1], fine, p[2], p[3], p[4])["deta"] <=
                    _ch_job(m, p[1], coarse, p[2], p[3], p[4])["deta"] - 0.2 for p in pairs);
           detail = m -> _ch_list(m, "deta")),
     Check("local-inertial |Δη| ≤ SWE |Δη| + 0.15 m at $(round(_CH_CELL[fine]; digits = 1)) m",
           m -> all(abs(_ch_job(m, mesh, fine, α, of, "li")["deta"]) <=
                    abs(_ch_job(m, mesh, fine, α, of, "swe")["deta"]) + 0.15
                    for mesh in ("a5", "cart"), α in angles, of in offsets);
           detail = m -> _ch_list(m, "deta"; filt = j -> _ch_res(j) == fine)),
     Check("SWE stays in bank (spill ≤ 10 m); nothing spills at $(round(_CH_CELL[fine]; digits = 1)) m",
           m -> all((endswith(j, "_swe") ? d["spill"] <= 10 : true) && (_ch_res(j) == fine ? d["spill"] == 0 : true)
                    for (j, d) in _ch_runs(m));
           detail = m -> _ch_list(m, "spill"))]
end

"The full series' checks, on HYDRAULICS §15.12's table (rms of Δη over the orientations, each the mean over the offsets)."
function _series_checks(; res)
    rms(m, s, mesh, r) = m["series"]["$(s)_$(mesh)_r$(r)_rms"]
    combos = [(s, mesh) for s in ("swe", "li"), mesh in ("a5", "cart")]
    table(m) = join(("$s $mesh " * join((round(rms(m, s, mesh, r); digits = 2) for r in res), "/")
                     for (s, mesh) in combos), "; ")
    [Check("mass conserved (|error| < 1e-5 of the volume added, from the log)",
           m -> all(d["mass_balance_log_rel"] < 1e-5 for (_, d) in _ch_runs(m));
           detail = m -> "worst " * string(maximum(d["mass_balance_log_rel"] for (_, d) in _ch_runs(m)))),
     Check("job arguments = the generator's jobs.txt", m -> m["series"]["jobs_match_generator"] == 1),
     Check("steady at 3 h: outflow/inflow within 1% from 11 m, 15% at 44 and 22 m",
           m -> all(abs(d["outflow"] - 1) <= (_ch_res(j) >= 19 ? 0.01 : 0.15) for (j, d) in _ch_runs(m));
           detail = m -> "off by > 0.5%: " * join(("$(j) $(round(d["outflow"]; digits = 3))" for (j, d) in _ch_runs(m)
                                          if abs(d["outflow"] - 1) > 0.005), ", ")),
     Check("rms Δη over orientations falls at every refinement (each solver × mesh)",
           m -> all(rms(m, s, mesh, res[i+1]) < rms(m, s, mesh, res[i]) for (s, mesh) in combos, i in 1:length(res)-1);
           detail = table),
     Check("rms Δη ≤ 0.2 m at 5.6 m cells (§15.12: 0.06–0.15 m)",
           m -> all(rms(m, s, mesh, 20) <= 0.2 for (s, mesh) in combos); detail = table),
     Check("local-inertial rms Δη ≤ SWE's + 0.05 m at 11 and 5.6 m",
           m -> all(rms(m, "li", mesh, r) <= rms(m, "swe", mesh, r) + 0.05 for mesh in ("a5", "cart"), r in (19, 20));
           detail = table),
     Check("SWE stays in bank (spill ≤ 10 m); nothing spills from 11 m",
           m -> all((endswith(j, "_swe") ? d["spill"] <= 10 : true) && (_ch_res(j) >= 19 ? d["spill"] == 0 : true)
                    for (j, d) in _ch_runs(m));
           detail = m -> "spilling: " * join(("$(j) $(round(d["spill"]; digits = 1))" for (j, d) in _ch_runs(m)
                                              if d["spill"] > 0), ", "))]
end

let res = [18, 19], angles = [0, 45], offsets = [0]
    register!(Case(
        name = "channel", tier = :standard,
        title = "Straight channel, 22 and 11 m cells, 0° and 45° (A5 and Cartesian, SWE and local-inertial)",
        doc = "docs/TEST_CASES.md §7; HYDRAULICS.md §15.12",
        description = """
        A subset of the channel series for routine use: res 18 and 19 (cells
        22 and 11 m, i.e. the 20 m channel ~1 and ~2 cells wide), orientations
        0° and 45°, no cross-channel offset, A5 and Cartesian, full SWE and the
        local-inertial default — 16 runs of 3 h. Metrics per run (FA channel
        metrics): Δη `deta` (mean in-channel stage error against normal flow,
        m), `spill` (wet floodplain per metre of reach, m), `unsteady` (max
        stage change over the last 900 s, m) and `outflow` (fraction of the
        inflow leaving over the last 900 s, from the log); under "series", the
        rms and spread of Δη over the orientations per solver × mesh × res.
        The checks are HYDRAULICS §15.12's findings: Δη falls with resolution
        (0.7–1.3 m at 22 m to 0.07–0.7 m at 11 m), the runs are steady, the
        local-inertial default is as close to normal flow as SWE where the
        channel is resolved, and SWE does not spill. Res 17 is left out
        because its Δη does not fall towards res 18 on Cartesian at 0°
        (0.74 m at both). 30 min on one 4-thread run at a time (70–270 s per run,
        about a minute of it Julia start-up), ~15 min with -P 2 -T 4.""",
        prep = [_channel_prep(; res, angles, offsets)],
        jobs = _channel_jobs(; res, angles, offsets),
        analyse = _channel_metrics,
        checks = _channel_checks(; res, angles, offsets),
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end

let res = [17, 18, 19, 20], angles = [0, 15, 30, 45, 60, 75, 90], offsets = [0, 0.5]
    register!(Case(
        name = "channel_series", tier = :full,
        title = "Straight-channel series: 4 resolutions × 7 orientations × 2 offsets, A5 and Cartesian, SWE and LI",
        doc = "docs/TEST_CASES.md §7; HYDRAULICS.md §15.12",
        description = """
        The whole series of HYDRAULICS §15.12: res 17–20 (cells 44, 22, 11,
        5.6 m), orientations 0–90° in 15° steps, offsets 0 and ½ cell, A5 and
        Cartesian, full SWE and the local-inertial default: 224 runs of 3 h,
        named and run as test/channel_bench/make_channel_meshes.jl's jobs.txt
        (the series' `limd` runs, `--momentum-damping 0.05`, a sensitivity of
        an option that is off by default, are not part of the case). Metrics
        as the `channel` case; "series" holds §15.12's table (rms and spread
        of Δη over the orientations). Preparation 11 min (0.46 GB of meshes).
        Many hours on CPU: on 4 threads res 18 runs take ~1 min, res 19 1–4.5
        min (SWE the longest); res 20 (130 000 cells) is estimated at ~5 min
        local-inertial and ~30 min SWE — the recorded series ran res 20 on a
        GPU (`--backend gpu`). The checks are provisional: derived from, and
        passing on, the recorded series' metrics (0.3.0). No expected values
        are recorded until the case has been run in full.""",
        prep = [_channel_prep(; res, angles, offsets)],
        jobs = _channel_jobs(; res, angles, offsets),
        analyse = _channel_metrics,
        checks = _series_checks(; res),
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
