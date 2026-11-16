# cases/registry/sgs_strip.jl — the SGS strip: a sub-cell channel in a long
# sloping flume, where every answer is analytic (docs/TEST_CASES.md §2c).
#
# Why: SGS runs store far more water than local-inertial and SWE on Carlisle
# (26.7 against 8.0 and 10.5 Mm³ at the peak) and are deeper at the inflow of
# the planar embankment. Is that how storage is measured (SGS keeps the water
# below the cell-mean bed that a standard grid smooths away), or is it lost
# conveyance? Here the geometry is simple enough to answer both from first
# principles.
#
# Terrain (cases/scripts/sgs_strip_dem.jl, 0.5 m DEM; geometry and analytic
# hydraulics in cases/scripts/sgs_strip_geometry.jl): a 4000 × 60 m flume
# sloping S = 0.001 along x (east), a 5 m × 1 m rectangular channel on its
# centre line, no walls in the DEM (`--wall 0`: the floodplain plane carries
# on beyond the strip); the mesh edges are closed (`--bc-file`) except the
# downstream end. Below
# x = 3800 m the bed steepens to 0.021: the flow passes through critical
# depth at the break, so no outflow boundary formulation (each model has
# its own) can reach the measuring reach x = 1000–2000 m. Manning
# n = 0.03 everywhere. Two discharges: 3 m³/s (in bank: the 2-D normal depth
# in the channel is 0.73 m) and 20 m³/s (overbank: 0.37 m on the floodplain).
#
# Meshes: Cartesian 20 m (3 rows; the channel is a quarter of the middle row)
# with SGS tables and with cell-mean beds; A5 res 18 (~22 m pentagons, a
# ragged strip ≈ 95 m wide on average: every cell the AOI touches) likewise; a Cartesian 2.5 m
# reference (the channel two cells wide) with cell-mean beds.
#
# Analytic answers (no simulation), level δ above the floodplain:
#   2-D       Σ strips h^(5/3): what a fine 2-D model converges to (no lateral
#             exchange at uniform flow); δ = −0.267 / +0.372 m
#   cells     rows of 20 m with their mean bed (a standard 20 m grid):
#             δ = +0.045 / +0.409 m
#   lumped    rows of 20 m with their true section, A^(5/3)/P_w^(2/3) — the
#             SGS face formula with one hydraulic radius (--sgs-conveyance
#             lumped): δ = −0.254 / +0.409. The default (2d) sums h^(5/3)
#             over the face profile and should give the 2-D level.
#   storage   per metre of strip, true water volume at δ; above the cell-mean
#             bed it is the same once δ ≥ 0 (the bed's mean is only moved
#             below the floodplain by the channel's own volume)
#   arrival   kinematic shock: t(x) = 300 s + A_n·(x − x_in)/Q, A_n the
#             storage per metre at normal depth (the hydrograph ramps over
#             600 s, so its volume lags a step by 300 s)
# The SGS tables are compared with the terrain directly
# (cases/scripts/sgs_strip_tables.jl → {case}/tables_*.json).

const _SS_Q = (3.0, 20.0)
const _SS_N = 0.03
const _SS_T = 28800.0                    # coarse runs: 8 h (a 20 m cell-mean grid at 3 m³/s fills slowly: a thin sheet)
const _SS_T_FINE = 14400.0               # 2.5 m reference: steady after ~2 h
const _SS_DT_OUT = 300.0
const _SS_REACH = (1000.0, 2000.0)       # measuring reach (m from the upstream end)
const _SS_STATIONS = (1000.0, 2000.0)    # arrival-time stations; window ±50 m
const _SS_XIN, _SS_YIN = 31.25, 1.25     # inflow point (m): in the channel, the second 20 m cell
const _SS_LON0, _SS_LAT0, _SS_R = 0.0, 51.0, 6_371_000.0
const _SS_LAT = round(_SS_LAT0 + rad2deg(_SS_YIN / _SS_R); digits = 9)
const _SS_LON = round(_SS_LON0 + rad2deg(_SS_XIN / (_SS_R * cosd(_SS_LAT0))); digits = 9)
"The case's DEM parameters (`sgs_strip_dem.jl` arguments; its defaults)."
const _SS_DEM = (L = 4000.0, B = 60.0, W = 5.0, D = 1.0, S = 0.001, P = 0.5, H = 0.0, YC = 0.0, SX = 0.0, Z0 = 10.0,
                 XB = 3800.0, S2 = 0.02)
const _SS_WET = 1e-3                     # m: a cell is wet when volume/area exceeds this

include(joinpath(@__DIR__, "..", "scripts", "sgs_strip_geometry.jl"))

_ss_params_file(dir) = _ss_params(FA.FO.JSON3.read(read(joinpath(dir, "sgs_strip_params.json"), String), Dict{String,Any}))
_ss_tables(dir, mesh) = FA.FO.JSON3.read(read(joinpath(dir, "tables_$mesh.json"), String), Dict{String,Any})
_ss_qtag(Q) = "q" * string(round(Int, Q))

"""
    _ss_analytic(p) -> Dict(metric => value)

The analytic answers for both discharges (see the file header): normal level δ,
channel depth, storage per metre, arrival times, for the 2-D, `cells` and
`lumped` (20 m rows) sections; with `width_a5` (the A5 strip's mean width,
m) also `a5w`, the 2-D answer for a strip of that width.
"""
function _ss_analytic(p; width_a5 = NaN)
    d = Dict{String,Float64}()
    for Q in _SS_Q
        t = _ss_qtag(Q)
        secs = Any[("2d", :twod, (;)), ("cells20", :cells, (rows = 20.0,)), ("lumped20", :lumped, (rows = 20.0,))]
        isfinite(width_a5) && push!(secs, ("a5w", :twod, (ylim = width_a5 / 2,)))
        for (lab, mode, kw) in secs
            δ = _ss_normal_level(p, Q, _SS_N, mode; kw...)
            st = mode === :cells ? _ss_storage(p, δ, :cells; kw...) :
                 _ss_storage(p, δ; (haskey(kw, :ylim) ? (ylim = kw.ylim,) : (;))...)
            d["delta_$(lab)_$t"] = δ
            d["store_$(lab)_$t"] = st
            for xs in _SS_STATIONS
                d["t_arrive_$(lab)_$(t)_x$(round(Int, xs))"] = 300.0 + st * (xs - _SS_XIN) / Q
            end
        end
        d["depth_channel_2d_$t"] = d["delta_2d_$t"] + p.D
    end
    return d
end

"Cell x, y (m) in the strip's frame."
_ss_xy(r) = FA.local_xy(r.mesh, _SS_LON0, _SS_LAT0; R = _SS_R)

"Channel cells: polygons overlapping the channel's central half (|y − YC| < W/4)."
_ss_channel_cells(r, p) = map(r.mesh.polygons) do P
    y = (view(P, 2, :) .- _SS_LAT0) .* (π / 180 * _SS_R)
    minimum(y) < p.YC + p.W / 4 && maximum(y) > p.YC - p.W / 4
end

"""
    _ss_profile(r, k, p) -> NamedTuple

Per cell at frame k: x, y, the water level above the local floodplain
δ = η − z_f(x) (NaN where dry), volume, the volume above the cell-mean bed
area·max(η − z̄, 0), and the channel-cell mask.
"""
function _ss_profile(r, k, p)
    x, y = _ss_xy(r)
    V = FA.volume(r, k); η = FA.wse(r, k); A = r.mesh.area
    wet = [isfinite(V[i]) && V[i] > _SS_WET * A[i] for i in eachindex(V)]
    δ = [wet[i] ? η[i] - _ss_zf(p, x[i]) : NaN for i in eachindex(V)]
    above = [wet[i] ? A[i] * max(η[i] - r.mesh.elevation[i], 0.0) : 0.0 for i in eachindex(V)]
    return (; x, y, δ, V = [isfinite(v) ? v : 0.0 for v in V], above, wet, chan = _ss_channel_cells(r, p))
end

"Least-squares slope of v against x (per km)."
function _ss_slope(x, v)
    isempty(x) && return NaN
    xm = sum(x) / length(x); vm = sum(v) / length(v)
    den = sum((xi - xm)^2 for xi in x)
    return den == 0 ? NaN : 1000 * sum((x[i] - xm) * (v[i] - vm) for i in eachindex(x)) / den
end

function _ss_run_metrics(r, p, Q)
    n = FA.nframes(r)
    pr = _ss_profile(r, n, p)
    a, b = _SS_REACH
    inr = [a <= pr.x[i] < b for i in eachindex(pr.x)]
    ch = [i for i in eachindex(pr.x) if inr[i] && pr.chan[i] && pr.wet[i]]
    fp = [i for i in eachindex(pr.x) if inr[i] && !pr.chan[i] && pr.wet[i]]
    wetr = [i for i in eachindex(pr.x) if inr[i] && pr.wet[i]]
    mean(v) = isempty(v) ? NaN : sum(v) / length(v)
    d = Dict{String,Float64}(
        "t_end" => r.t[n], "stored_m3" => sum(pr.V),
        "delta_channel" => mean(pr.δ[ch]),
        "delta_floodplain" => mean(pr.δ[fp]),
        "delta_spread" => isempty(wetr) ? NaN : maximum(pr.δ[wetr]) - minimum(pr.δ[wetr]),
        "delta_slope_per_km" => _ss_slope(pr.x[ch], pr.δ[ch]),
        "store_per_m" => sum(pr.V[inr]) / (b - a),
        "store_above_mean_per_m" => sum(pr.above[inr]) / (b - a),
        "n_wet_reach" => length(wetr))
    # steadiness: storage change over the last hour against the inflow
    k1 = FA.frame_at(r, r.t[n] - 3600.0; tol = _SS_DT_OUT)
    d["dSdt_rel_end"] = (FA.stored_volume(r, n) - FA.stored_volume(r, k1)) / (r.t[n] - r.t[k1]) / Q
    # arrival: the ±50 m window around each station holds half the 2-D normal storage
    st2d = _ss_storage(p, _ss_normal_level(p, Q, _SS_N, :twod))
    for xs in _SS_STATIONS
        win = [i for i in eachindex(pr.x) if xs - 50 <= pr.x[i] < xs + 50]
        t_half = Inf
        for k in 1:n
            v = FA.volume(r, k)
            if sum(isfinite(v[i]) ? v[i] : 0.0 for i in win) / 100.0 >= 0.5st2d
                t_half = r.t[k]; break
            end
        end
        d["t_arrive_x$(round(Int, xs))"] = t_half
    end
    return d
end

function _ss_table_metrics(t)
    v(s, δ, stat) = Float64(get(t[s][string(Float64(δ))], stat, NaN))
    e(c, δ, stat) = Float64(get(t["edge_area_ratio"][c][string(Float64(δ))], stat, NaN))
    k(c, δ, stat) = Float64(get(t["edge_conveyance_ratio"][c][string(Float64(δ))], stat, NaN))
    d = Dict{String,Float64}(
        "volume_err_max" => maximum(v("volume_rel_err", δ, "max_abs") for δ in t["levels"]),
        "volume_err_median_max" => maximum(abs(v("volume_rel_err", δ, "median")) for δ in t["levels"]),
        "edge_A_ratio_channel_inbank" => e("channel", -0.25, "median"),
        "edge_A_ratio_channel_overbank" => e("channel", 0.4, "median"),
        "edge_A_ratio_floodplain_overbank" => e("floodplain", 0.4, "median"),
        "edge_A_ratio_channel_1m" => e("channel", 1.0, "median"),
        "edge_A_ratio_floodplain_1m" => e("floodplain", 1.0, "median"),
        "edge_K_ratio_channel_inbank" => k("channel", -0.25, "median"),
        "edge_K_ratio_channel_overbank" => k("channel", 0.4, "median"),
        "edge_K_ratio_floodplain_overbank" => k("floodplain", 0.4, "median"),
        "edge_K_ratio_channel_1m" => k("channel", 1.0, "median"))
    u = t["uniform_flow"]
    u === nothing || (d["uniform_Q_max"] = Float64(u["Q_max"]); d["uniform_delta_at_Q_max"] = Float64(u["delta_at_Q_max"]))
    return d
end

function _sgs_strip_metrics(ctx)
    p = _ss_params_file(ctx.dir)
    out = Dict{String,Dict{String,Float64}}()
    # the A5 strip's mean width over the reach (its edges are ragged): Σ cell area / reach length
    ra = ctx.runs["a5r18_std_q3"]; xa, _ = _ss_xy(ra)
    w_a5 = sum(ra.mesh.area[i] for i in eachindex(xa) if _SS_REACH[1] <= xa[i] < _SS_REACH[2]) / (_SS_REACH[2] - _SS_REACH[1])
    out["analytic"] = _ss_analytic(p; width_a5 = w_a5)
    out["analytic"]["width_a5"] = w_a5
    for mesh in ("cart20", "a5r18")
        out["tables_$mesh"] = _ss_table_metrics(_ss_tables(ctx.dir, mesh))
    end
    for (name, r) in ctx.runs
        out[name] = _ss_run_metrics(r, p, Float64(ctx.jobs[name].tags[:Q]))
    end
    return out
end

# ── checks ──
_ss(m, j, k) = get(get(m, j, Dict{String,Float64}()), k, NaN)
_ss_r(x; d = 3) = round(x; digits = d)
_ss_jobs(m, pre) = sort([j for j in keys(m) if startswith(j, pre)])
"Check that a run's channel level is within tol (m) of an analytic level and its storage within rtol."
_ss_level_ok(m, job, ref; tol, rtol) =
    abs(_ss(m, job, "delta_channel") - _ss(m, "analytic", "delta_$ref")) <= tol &&
    abs(_ss(m, job, "store_per_m") / _ss(m, "analytic", "store_$ref") - 1) <= rtol
_ss_level_detail(m, job, ref) =
    "$job: δ $(_ss_r(_ss(m, job, "delta_channel"))) m (analytic $(ref) $(_ss_r(_ss(m, "analytic", "delta_$ref")))), " *
    "storage $(_ss_r(_ss(m, job, "store_per_m"); d = 2)) m²/m (analytic $(_ss_r(_ss(m, "analytic", "store_$ref"); d = 2)))"
_ss_steady(m, job) = abs(_ss(m, job, "dSdt_rel_end")) <= 0.01
"Arrival time at station xs: of a run, or (job = \"analytic\") the 2-D kinematic-shock time."
_ss_t(m, job, q, xs) = job == "analytic" ? _ss(m, "analytic", "t_arrive_2d_$(q)_x$(round(Int, xs))") :
                                           _ss(m, job, "t_arrive_x$(round(Int, xs))")
_ss_arrive_ok(m, job, ref, q; rtol = 0.15) =
    all(abs(_ss_t(m, job, q, xs) / _ss_t(m, ref, q, xs) - 1) <= rtol for xs in _SS_STATIONS)
_ss_fmt_t(t) = isfinite(t) ? string(round(Int, t)) : "never"
_ss_arrive_detail(m, job, ref, q) = "$job: " * join(("x=$(round(Int, xs)) $(_ss_fmt_t(_ss_t(m, job, q, xs))) s " *
    "($ref $(_ss_fmt_t(_ss_t(m, ref, q, xs))))" for xs in _SS_STATIONS), ", ")

# Checks of the correct SGS behaviour that the model still fails: kept
# unweakened and marked as expected failures (matched by name prefix). Until
# 0.7.0 there were eight (frozen overbank face tables, the legacy donor cap:
# cases/expected/v0.6.0/); the face tables, positivity limiting and the
# face-area time-step limit of 0.7.0 fixed them.
# (The eighth, A5 in bank, was 0.054 m high with the face tables alone: on A5
# the channel zig-zags through the cells, and the kernel took the slope over
# the centre distance and the faces' full oblique cut. The flow-path geometry
# of 0.7.0, _sgs_path, takes both along the line between the wetted areas.)
const _SS_XFAIL = Pair{String,String}[]
function _ss_mark_xfail(chs)
    map(chs) do ch
        i = findfirst(p -> startswith(ch.name, first(p)), _SS_XFAIL)
        i === nothing ? ch : Check(ch.name, ch.test, ch.detail, last(_SS_XFAIL[i]))
    end
end

_ss_checks() = _ss_mark_xfail(_ss_checks_all())

function _ss_checks_all()
    q3, q20 = _ss_qtag(3.0), _ss_qtag(20.0)
    runs(m) = [j for j in keys(m) if !startswith(j, "analytic") && !startswith(j, "tables_")]
    [
     Check("mass conserved in every run (|error| < 1e-9 of the volume added)",
           m -> all(_ss(m, j, "mass_balance_rel") < 1e-9 for j in runs(m));
           detail = m -> "worst " * string(maximum(_ss(m, j, "mass_balance_rel") for j in runs(m)))),
     # — SGS tables against the terrain (no simulation) —
     Check("SGS tables, Cartesian 20 m: stored volume V(η) within 1% of the analytic terrain at every level " *
           "(δ = −0.75 … +1 m; error relative to max(V, 5 cm over the cell))",
           m -> _ss(m, "tables_cart20", "volume_err_max") <= 0.01;
           detail = m -> "max $(_ss_r(_ss(m, "tables_cart20", "volume_err_max"); d = 4))"),
     Check("SGS tables, A5 res 18: median V(η) error within 1% at every level (cells holding a sliver of the channel carry Halton error)",
           m -> _ss(m, "tables_a5r18", "volume_err_median_max") <= 0.01;
           detail = m -> "worst median $(_ss_r(_ss(m, "tables_a5r18", "volume_err_median_max"); d = 4)), " *
                         "max $(_ss_r(_ss(m, "tables_a5r18", "volume_err_max"); d = 3))"),
     Check("SGS tables, Cartesian 20 m: in-bank (δ = −0.25 m) channel-face flow area within 10% of the face's analytic area",
           m -> abs(_ss(m, "tables_cart20", "edge_A_ratio_channel_inbank") - 1) <= 0.1;
           detail = m -> "A_edge/A = $(_ss_r(_ss(m, "tables_cart20", "edge_A_ratio_channel_inbank")))"),
     Check("SGS tables, Cartesian 20 m: overbank (δ = +0.4 m) face flow areas within 10% of the analytic areas " *
           "(channel and floodplain faces)",
           m -> abs(_ss(m, "tables_cart20", "edge_A_ratio_channel_overbank") - 1) <= 0.1 &&
                abs(_ss(m, "tables_cart20", "edge_A_ratio_floodplain_overbank") - 1) <= 0.1;
           detail = m -> "A_edge/A: channel faces $(_ss_r(_ss(m, "tables_cart20", "edge_A_ratio_channel_overbank"))), " *
                         "floodplain faces $(_ss_r(_ss(m, "tables_cart20", "edge_A_ratio_floodplain_overbank"))) " *
                         "(at δ = +1 m: $(_ss_r(_ss(m, "tables_cart20", "edge_A_ratio_channel_1m"))), " *
                         "$(_ss_r(_ss(m, "tables_cart20", "edge_A_ratio_floodplain_1m"))))"),
     Check("SGS tables, A5 res 18: overbank (δ = +0.4 m) face flow areas within 10% of the analytic areas (medians)",
           m -> abs(_ss(m, "tables_a5r18", "edge_A_ratio_channel_overbank") - 1) <= 0.1 &&
                abs(_ss(m, "tables_a5r18", "edge_A_ratio_floodplain_overbank") - 1) <= 0.1;
           detail = m -> "A_edge/A: channel faces $(_ss_r(_ss(m, "tables_a5r18", "edge_A_ratio_channel_overbank"))), " *
                         "floodplain faces $(_ss_r(_ss(m, "tables_a5r18", "edge_A_ratio_floodplain_overbank")))"),
     Check("SGS tables, Cartesian 20 m: face conveyance Σ h^(5/3) within 5% of the face's analytic value in bank, " *
           "overbank (δ = +0.4 m) and at δ = +1 m (above the tabulated profile: the continuation)",
           m -> all(abs(_ss(m, "tables_cart20", "edge_K_ratio_$c") - 1) <= 0.05
                    for c in ("channel_inbank", "channel_overbank", "floodplain_overbank", "channel_1m"));
           detail = m -> "K/K_analytic: channel faces $(_ss_r(_ss(m, "tables_cart20", "edge_K_ratio_channel_inbank"))) in bank, " *
                         "$(_ss_r(_ss(m, "tables_cart20", "edge_K_ratio_channel_overbank"))) overbank, " *
                         "$(_ss_r(_ss(m, "tables_cart20", "edge_K_ratio_channel_1m"))) at +1 m; floodplain faces " *
                         "$(_ss_r(_ss(m, "tables_cart20", "edge_K_ratio_floodplain_overbank")))"),
     Check("SGS tables, Cartesian 20 m: at the bed slope they convey 20 m³/s (uniform flow through one cross-section)",
           m -> _ss(m, "tables_cart20", "uniform_Q_max") >= 20.0;
           detail = m -> "at most $(_ss_r(_ss(m, "tables_cart20", "uniform_Q_max"); d = 2)) m³/s " *
                         "(at δ = $(_ss_r(_ss(m, "tables_cart20", "uniform_delta_at_Q_max"); d = 2)) m)"),
     # — references —
     Check("2.5 m references (local-inertial and SWE): steady, channel level within 0.03 m and storage within 3% of the 2-D analytic",
           m -> all(_ss_steady(m, "fine_$(s)_$q") && _ss_level_ok(m, "fine_$(s)_$q", "2d_$q"; tol = 0.03, rtol = 0.03)
                    for s in ("std", "swe"), q in (q3, q20));
           detail = m -> join((_ss_level_detail(m, "fine_$(s)_$q", "2d_$q") for s in ("std", "swe"), q in (q3, q20)), "; ")),
     Check("2.5 m references, in bank: arrival at x = 1000 and 2000 m within 10% of the kinematic-shock time " *
           "(overbank the floodplain fills behind a faster in-bank front, so the single-shock time is only a lower " *
           "bound; the 2.5 m local-inertial run is the reference there)",
           m -> all(_ss_arrive_ok(m, "fine_$(s)_$q3", "analytic", q3; rtol = 0.10) for s in ("std", "swe"));
           detail = m -> join((_ss_arrive_detail(m, "fine_$(s)_$q", "analytic", q) for s in ("std", "swe"), q in (q3, q20)), "; ")),
     Check("Cartesian 20 m local-inertial and SWE (cell-mean beds): steady, level within 0.03 m and storage within 3% " *
           "of the cell-mean analytic (in bank: +0.045 m, 2.2× the true volume — the channel smeared over its row)",
           m -> all(_ss_steady(m, "cart20_$(s)_$q") && _ss_level_ok(m, "cart20_$(s)_$q", "cells20_$q"; tol = 0.03, rtol = 0.03)
                    for s in ("std", "swe"), q in (q3, q20));
           detail = m -> join((_ss_level_detail(m, "cart20_$(s)_$q", "cells20_$q") for s in ("std", "swe"), q in (q3, q20)), "; ")),
     # — SGS: the correct behaviour —
     Check("SGS Cartesian 20 m, in bank: steady, channel level within 0.03 m of the 2-D analytic " *
           "and storage within 5% of the true volume",
           m -> _ss_steady(m, "cart20_sgs_$q3") && _ss_level_ok(m, "cart20_sgs_$q3", "2d_$q3"; tol = 0.03, rtol = 0.05);
           detail = m -> _ss_level_detail(m, "cart20_sgs_$q3", "2d_$q3")),
     Check("SGS Cartesian 20 m, overbank: steady, level within 0.03 m of the 2-D analytic and storage within 5%",
           m -> _ss_steady(m, "cart20_sgs_$q20") && _ss_level_ok(m, "cart20_sgs_$q20", "2d_$q20"; tol = 0.03, rtol = 0.05);
           detail = m -> _ss_level_detail(m, "cart20_sgs_$q20", "2d_$q20") *
                         ", dS/dt/Q $(_ss_r(_ss(m, "cart20_sgs_$q20", "dSdt_rel_end")))"),
     Check("SGS Cartesian 20 m, --sgs-conveyance lumped (one hydraulic radius per face): steady, level within 0.03 m " *
           "and storage within 5% of the lumped analytic at both discharges (overbank +0.04 m, +8% against 2-D: " *
           "the level of a cell-mean grid)",
           m -> all(_ss_steady(m, "cart20_sgslumped_$q") &&
                    _ss_level_ok(m, "cart20_sgslumped_$q", "lumped20_$q"; tol = 0.03, rtol = 0.05) for q in (q3, q20));
           detail = m -> join((_ss_level_detail(m, "cart20_sgslumped_$q", "lumped20_$q") for q in (q3, q20)), "; ")),
     Check("SGS Cartesian 20 m: arrival at x = 1000 and 2000 m within 15% of the 2.5 m local-inertial reference (both discharges)",
           m -> all(_ss_arrive_ok(m, "cart20_sgs_$q", "fine_std_$q", q) for q in (q3, q20));
           detail = m -> join((_ss_arrive_detail(m, "cart20_sgs_$q", "fine_std_$q", q) for q in (q3, q20)), "; ")),
     Check("A5 res 18 local-inertial (cell-mean beds): steady at both discharges; overbank level within 0.15 m " *
           "of the 2-D analytic for the A5 strip's mean width (its ragged edges hold slack water)",
           m -> _ss_steady(m, "a5r18_std_$q3") && _ss_steady(m, "a5r18_std_$q20") &&
                abs(_ss(m, "a5r18_std_$q20", "delta_channel") - _ss(m, "analytic", "delta_a5w_$q20")) <= 0.15;
           detail = m -> _ss_level_detail(m, "a5r18_std_$q3", "a5w_$q3") * "; " * _ss_level_detail(m, "a5r18_std_$q20", "a5w_$q20") *
                         " (mean width $(_ss_r(_ss(m, "analytic", "width_a5"); d = 1)) m)"),
     Check("SGS A5 res 18, in bank: steady, channel level within 0.05 m of the 2-D analytic and storage within 10%",
           m -> _ss_steady(m, "a5r18_sgs_$q3") && _ss_level_ok(m, "a5r18_sgs_$q3", "2d_$q3"; tol = 0.05, rtol = 0.10);
           detail = m -> _ss_level_detail(m, "a5r18_sgs_$q3", "2d_$q3") * ", dS/dt/Q $(_ss_r(_ss(m, "a5r18_sgs_$q3", "dSdt_rel_end")))"),
     Check("SGS A5 res 18, overbank: steady; level and storage no further from the 2-D analytic for the A5 strip's " *
           "mean width than local-inertial on the same mesh (its ragged edges hold slack water in both)",
           m -> _ss_steady(m, "a5r18_sgs_$q20") &&
                abs(_ss(m, "a5r18_sgs_$q20", "delta_channel") - _ss(m, "analytic", "delta_a5w_$q20")) <=
                    abs(_ss(m, "a5r18_std_$q20", "delta_channel") - _ss(m, "analytic", "delta_a5w_$q20")) &&
                abs(_ss(m, "a5r18_sgs_$q20", "store_per_m") - _ss(m, "analytic", "store_a5w_$q20")) <=
                    abs(_ss(m, "a5r18_std_$q20", "store_per_m") - _ss(m, "analytic", "store_a5w_$q20"));
           detail = m -> "SGS δ $(_ss_r(_ss(m, "a5r18_sgs_$q20", "delta_channel"))) m, storage $(_ss_r(_ss(m, "a5r18_sgs_$q20", "store_per_m"); d = 2)) m²/m, " *
                         "dS/dt/Q $(_ss_r(_ss(m, "a5r18_sgs_$q20", "dSdt_rel_end"))); LI δ $(_ss_r(_ss(m, "a5r18_std_$q20", "delta_channel"))) m, " *
                         "storage $(_ss_r(_ss(m, "a5r18_std_$q20", "store_per_m"); d = 2)) m²/m; analytic " *
                         "$(_ss_r(_ss(m, "analytic", "delta_a5w_$q20"))) m, $(_ss_r(_ss(m, "analytic", "store_a5w_$q20"); d = 2)) m²/m"),
     # — the storage measurement question —
     Check("measurement: overbank, the volume above the cell-mean bed equals the water volume (within 1%) — " *
           "SGS and the 2.5 m reference",
           m -> all(abs(_ss(m, j, "store_above_mean_per_m") / _ss(m, j, "store_per_m") - 1) <= 0.01
                    for j in ("cart20_sgs_$q20", "cart20_sgslumped_$q20", "fine_std_$q20"));
           detail = m -> join(("$j $(_ss_r(_ss(m, j, "store_above_mean_per_m"); d = 2))/$(_ss_r(_ss(m, j, "store_per_m"); d = 2)) m²/m"
                               for j in ("cart20_sgs_$q20", "cart20_sgslumped_$q20", "fine_std_$q20")), ", ")),
    ]
end

let
    dem = ["{repo}/cases/scripts/sgs_strip_dem.jl", "{case}",
           "--length", string(_SS_DEM.L), "--width", string(_SS_DEM.B), "--channel-width", string(_SS_DEM.W),
           "--channel-depth", string(_SS_DEM.D), "--slope", string(_SS_DEM.S), "--pixel", string(_SS_DEM.P),
           "--wall", string(_SS_DEM.H), "--channel-offset", string(_SS_DEM.YC), "--cross-slope", string(_SS_DEM.SX),
           "--z0", string(_SS_DEM.Z0), "--break-x", string(_SS_DEM.XB), "--break-slope", string(_SS_DEM.S2)]
    aoi, demf = "{case}/sgs_strip_aoi.geojson", "{case}/sgs_strip_dem.tif"
    cart(dx, fm, name) = Prep(name, :floodmodel,
        ["--meshgen-cartesian", aoi, "--dx", string(dx), "--cartesian-lat-ref", string(_SS_LAT0), "--dem", demf,
         "--meshout", "{case}/$name.parquet", "--flow-model", fm, "--mesh-only"]; outputs = ["{case}/$name.parquet"])
    a5(fm, name) = Prep(name, :floodmodel,
        ["--meshgen", aoi, "--meshres", "18", "--dem", demf, "--meshout", "{case}/$name.parquet",
         "--flow-model", fm, "--mesh-only"]; outputs = ["{case}/$name.parquet"])
    tables(mesh) = Prep("tables_$mesh", :julia,
        ["{repo}/cases/scripts/sgs_strip_tables.jl", "{case}/$(mesh)_sgs.parquet", "{case}/sgs_strip_params.json",
         "{case}/tables_$mesh.json"]; outputs = ["{case}/tables_$mesh.json"])
    hydro(Q) = Prep("hydrograph_$(_ss_qtag(Q))", :write,
        ["{case}/$(_ss_qtag(Q)).csv", "t_s,Q_m3s\n0,0\n600,$Q\n1000000,$Q\n"]; outputs = ["{case}/$(_ss_qtag(Q)).csv"])
    prep = [Prep("dem", :julia, dem; outputs = ["{case}/sgs_strip_dem.tif", "{case}/sgs_strip_params.json",
                                                 "{case}/sgs_strip_aoi.geojson", "{case}/sgs_strip_bc.geojson"]),
            cart(20, "sgs", "cart20_sgs"), cart(20, "standard", "cart20_standard"),
            cart(2.5, "standard", "cart2p5_standard"),
            a5("sgs", "a5r18_sgs"), a5("standard", "a5r18_standard"),
            hydro(3.0), hydro(20.0), tables("cart20"), tables("a5r18")]
    args(mesh, fm, Q, T; extra = String[]) =
        vcat(["--meshload", "{case}/$mesh.parquet", "--flow-model", fm,
              "--inflow-point", "$(_SS_LAT),$(_SS_LON),{case}/$(_ss_qtag(Q)).csv",
              "--bc-file", "{case}/sgs_strip_bc.geojson", "--manning-n", string(_SS_N),
              "--sim-duration", string(round(Int, T)), "--dt-max", "10",
              "--output-interval", string(round(Int, _SS_DT_OUT))], extra)
    jobs = Job[]
    for Q in _SS_Q
        q = _ss_qtag(Q)
        push!(jobs, Job("cart20_sgs_$q", args("cart20_sgs", "sgs", Q, _SS_T); threads = 2, mesh = "cart20", model = "sgs", Q),
                    Job("cart20_sgslumped_$q", args("cart20_sgs", "sgs", Q, _SS_T; extra = ["--sgs-conveyance", "lumped"]);
                        threads = 2, mesh = "cart20", model = "sgs", conveyance = "lumped", Q),
                    Job("cart20_std_$q", args("cart20_standard", "standard", Q, _SS_T); threads = 2, mesh = "cart20", model = "standard", Q),
                    Job("cart20_swe_$q", args("cart20_standard", "swe", Q, _SS_T); threads = 2, mesh = "cart20", model = "swe", Q),
                    Job("a5r18_sgs_$q", args("a5r18_sgs", "sgs", Q, _SS_T); threads = 2, mesh = "a5r18", model = "sgs", Q),
                    Job("a5r18_std_$q", args("a5r18_standard", "standard", Q, _SS_T); threads = 2, mesh = "a5r18", model = "standard", Q),
                    Job("fine_std_$q", args("cart2p5_standard", "standard", Q, _SS_T_FINE); threads = 4, mesh = "cart2p5", model = "standard", Q),
                    Job("fine_swe_$q", args("cart2p5_standard", "swe", Q, _SS_T_FINE); threads = 4, mesh = "cart2p5", model = "swe", Q))
    end
    register!(Case(
        name = "sgs_strip", tier = :standard,
        title = "SGS strip: a sub-cell channel in a sloping flume, against analytic normal flow, storage and arrival",
        doc = "docs/TEST_CASES.md §2c",
        description = """
        A 4000 × 60 m strip on a 0.1% slope (0.5 m DEM, no walls) with a 5 m
        wide, 1 m deep rectangular channel on its centre line; the mesh edges
        closed except the downstream end, behind a steep (2.1%) last 200 m;
        Manning n = 0.03; a constant inflow (10-minute ramp) of 3 m³/s (in bank) or 20 m³/s (overbank) in the channel 31 m from
        the upstream end. Meshes: Cartesian 20 m (3 rows; the channel a quarter
        of a cell) with SGS tables or cell-mean beds, A5 res 18 likewise, and a
        Cartesian 2.5 m reference with the channel two cells wide. Runs: SGS
        (default, and --sgs-conveyance lumped), local-inertial and SWE on
        20 m; SGS and local-inertial on A5; local-inertial and SWE on 2.5 m — 8 h
        (2.5 m: 4 h), output every 5 min. Analytic answers (file header of
        cases/registry/sgs_strip.jl): the 2-D normal level (δ above the
        floodplain −0.267 / +0.372 m), the cell-mean and exact-table (lumped)
        20 m levels, storage per metre and kinematic-shock arrival times. The
        SGS tables are also checked against the terrain directly
        (tables_*.json, no simulation). Metrics per run over x = 1000–2000 m at
        the last frame: channel level δ, its spread and slope, storage per
        metre (water, and above the cell-mean bed), steadiness (storage change
        over the last hour / Q), arrival times at x = 1000 and 2000 m. The
        SGS checks state the correct behaviour. Until 0.7.0 eight of them
        failed (frozen overbank face tables, the legacy donor cap; record in
        cases/expected/v0.6.0/); all pass since 0.7.0. Preparation ~5 min; runs ~1 min each
        on 2 threads at 20 m and A5, 1.5–25 min on 4 threads at 2.5 m (SWE at
        20 m³/s the longest); ~35 min with -P 2.""",
        prep = prep, jobs = jobs,
        analyse = _sgs_strip_metrics,
        checks = _ss_checks(),
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
