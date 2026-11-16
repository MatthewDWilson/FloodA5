# cases/registry/planar_wall.jl — a wall with no gap across a planar slope: is a
# barrier that crosses the inside of a cell a barrier? (docs/TEST_CASES.md §2d;
# SGS conveyance fix plan §2.3, Phase 0 item 1.)
#
# Why: SGS stores water with a hypsometric curve per cell (one level pool) and
# conveys it through each face with the DEM profile on the face line. Nothing
# between the faces is used. A wall that crosses the inside of a cell lowers
# the cell's stored volume but should not stop water: it enters over the
# upstream faces, is stored at the cell's lowest ground — on a slope, the
# downstream side of the wall — and leaves through the downstream faces, whose
# profiles are natural ground. planar_embankment cannot show this (water may
# use its notch); here there is no notch, so in a correct model no water passes
# while the upstream level is below the crest.
#
# Terrain (generate_planar_embankment_dem.py --notch-width 0, one DEM per wall
# in {case}/wide and {case}/thin): the planar_embankment plane — 2 × 0.5 km,
# 0.5% west→east, z = 10 m at the west edge, 2 m DEM — with a wall across it,
# 3 m high above the ground on its centre line (flat top), no notch:
#   wide   40 m, x = 980–1020 m (20 pixels), crest 8.0 m: wider than a
#          res-18 A5 cell (~22 m) and two Cartesian 20 m cells
#   thin    4 m, x = 988–992 m (2 pixels), crest 8.05 m: narrower than any
#          cell, and inside the Cartesian 20 m column x = 981.1–1001.1 m, 7 m
#          or more from its faces. (At x = 998–1002 m it would lie on the
#          Cartesian face at 1001.1 m: a wall along a face, a different test —
#          docs/TEST_CASES.md §2d records that run.)
#
# Source and design. As planar_embankment: 10 m³/s (10-minute ramp) at
# x = 200 m on the centre line, closed boundaries, 6 h: 213 000 m³. Held as a
# level pool against the wall's west face (x_w, ground z_w there) that volume
# is a wedge
#     V(η) = ½ (η − z_w)² / S · W = 5·10⁴ (η − z_w)² m³   (S = 0.005, W = 500 m)
# so η − z_w = √(213 000 / 5·10⁴) = 2.064 m at the end: η = 7.164 m against
# the wide wall (z_w = 5.10 m) and 7.124 m against the thin one (5.06 m), 0.84
# and 0.93 m below the crests. (An upper bound: water still on its way down the
# slope is not in the pool. The wedge is 413 m long, inside the 780 m between
# the source and either wall.) The crest is never reached: any water east of the
# wall has crossed a closed barrier.
#
# Runs: SGS and standard flow (cell-mean beds), on A5 res 18 and Cartesian
# 20 m, for each wall (8). Standard flow is the comparison: it sees the wide
# wall at (nearly) its crest, and the thin wall as a low cell-mean ridge
# (≈ 3 m · 4 m / 22 m ≈ 0.5 m), which it overtops early — a different, known
# limitation of cell-mean beds, measured here but not what the case is for.

const _PW_DOMAIN_KM = (2.0, 0.5)
const _PW_SLOPE_PCT = 0.5
const _PW_HEIGHT = 3.0
const _PW_WALLS = ("wide" => (width = 40.0, x = 1000.0),    # wall => width, centre line (m)
                   "thin" => (width = 4.0, x = 990.0))
const _PW_PIXEL = 2.0
const _PW_LAT, _PW_LON = _planar_point(; domain_km = _PW_DOMAIN_KM)   # x = 200 m, centre line
const _PW_CENTRE_LAT = 51.0                 # the generator's default centre (Cartesian frame)
const _PW_Q = 10.0
const _PW_T = 21600.0
const _PW_DT_OUT = 120.0
const _PW_DOWN_THRESHOLD = 1.0          # m³ east of the wall = "water has passed"
const _PW_CELL_V = 0.01                 # m³: a cell downstream of the wall "holds water"
const _PW_WET = 0.01                    # m: a cell is wet above this depth
const _PW_MARGIN = 0.5                  # m: design freeboard required below the crest
const _PW_MESHES = ("a5" => "A5 res 18", "cart" => "Cartesian 20 m")
const _PW_MODELS = ("sgs" => "SGS", "std" => "standard")
_pw_job(mesh, model, wall) = "$(mesh)_$(model)_$(wall)"

"""
    _pw_design(wall) ->(; crest, z_foot, wse_end, freeboard)

The design of the case (file header), from its constants: the crest, the
ground at the wall's upstream face, the level-pool (wedge) level when the
whole hydrograph volume is held against the wall, and the freeboard left.
"""
function _pw_design(wall)
    S = _PW_SLOPE_PCT / 100; Lx, W = _PW_DOMAIN_KM .* 1000.0
    crest = S * (Lx - wall.x) + _PW_HEIGHT
    z_foot = S * (Lx - wall.x + wall.width / 2)
    wse_end = z_foot + sqrt(2 * _PW_Q * (_PW_T - 300.0) * S / W)
    return (; crest, z_foot, wse_end, freeboard = crest - wse_end)
end

"""
    _pw_run_metrics(r, p) -> Dict

One run against its wall (`p`: the generator's parameters). Cells are sorted
by their polygons: upstream (entirely west of the wall), downstream (entirely
east), straddling (the rest); `foot` the upstream cells within a cell of the
wall, `row` the downstream cells within a cell of it. y is from the south edge.
"""
function _pw_run_metrics(r, p)
    mesh = r.mesh
    crest = Float64(p["emb_crest_z"])
    xw = Float64(p["emb_x_m"]) - Float64(p["emb_width_m"]) / 2
    xe = Float64(p["emb_x_m"]) + Float64(p["emb_width_m"]) / 2
    W = Float64(p["domain_m_ns"])
    xmin, xmax = _pe_xrange(mesh, p)
    ymin, ymax = _planar_yrange(mesh, p)
    yc = (ymin .+ ymax) ./ 2
    cell = sqrt(sum(mesh.area) / length(mesh.area))
    up = xmax .<= xw
    down = xmin .>= xe
    strad = .!up .& .!down
    foot = up .& (xmax .>= xw - cell)
    row = down .& (xmin .<= xe + cell)
    n = FA.nframes(r)
    wse_up = fill(NaN, n); v_down = zeros(n); v_strad = zeros(n); row_wet = zeros(n)
    fin(V, i) = isfinite(V[i]) ? V[i] : 0.0
    for k in 1:n
        h = FA.depth(r, k); η = FA.wse(r, k); V = FA.volume(r, k)
        w = [i for i in findall(foot) if isfinite(h[i]) && h[i] > _PW_WET]
        isempty(w) || (wse_up[k] = maximum(η[w]))
        v_down[k] = sum(fin(V, i) for i in findall(down); init = 0.0)
        v_strad[k] = sum(fin(V, i) for i in findall(strad); init = 0.0)
        row_wet[k] = count(i -> fin(V, i) > _PW_CELL_V, findall(row)) / max(count(row), 1)
    end
    kf = findfirst(>(_PW_DOWN_THRESHOLD), v_down)
    ks = findfirst(>(_PW_DOWN_THRESHOLD), v_strad)
    stored = FA.stored_volume(r, n)
    d = Dict{String,Float64}(
        "inj_point_mismatch_deg" => max(abs(Float64(p["injection_lat"]) - _PW_LAT),
                                        abs(Float64(p["injection_lon"]) - _PW_LON)),
        "notch_width_m" => Float64(p["notch_width_m"]),
        "crest_z" => crest, "wall_width_m" => xe - xw,
        "z_foot" => Float64(p["z_west"]) - Float64(p["slope_pct"]) / 100 * xw,
        "cell_size_m" => cell,
        "n_up_foot" => count(foot), "n_straddle" => count(strad), "n_down" => count(down),
        "n_down_row" => count(row),
        "z_bneck_mean" => _pe_bottleneck(mesh.elevation, mesh.neighbours, up, down),
        "t_first_down" => kf === nothing ? Inf : r.t[kf],
        "t_first_straddle" => ks === nothing ? Inf : r.t[ks],
        "wse_up_at_first_down" => kf === nothing ? NaN : wse_up[kf],
        "wse_up_before_first_down" => (kf === nothing || kf == 1) ? NaN : wse_up[kf-1],
        "wse_up_end" => wse_up[n],
        "wse_up_max" => maximum(filter(isfinite, wse_up); init = -Inf),
        "v_down_end" => v_down[n], "v_straddle_end" => v_strad[n],
        "frac_down_end" => v_down[n] / stored,
        "row_wet_frac_end" => row_wet[n],
        "stored_m3" => stored, "t_end" => r.t[n])
    d["freeboard_min"] = crest - d["wse_up_max"]
    haskey(mesh.static, "sgs_z_min") &&
        (d["z_bneck_zmin"] = _pe_bottleneck(mesh.static["sgs_z_min"], mesh.neighbours, up, down))
    for th in 1:6
        k = FA.frame_at(r, 3600.0 * th; tol = 90.0)
        d["wse_up_$(th)h"] = wse_up[k]
        d["v_down_$(th)h"] = v_down[k]
        d["row_wet_frac_$(th)h"] = row_wet[k]
    end
    # where the water crossed: the downstream row's cells that hold water at the
    # first output with water downstream, against the stretch of the wall the
    # upstream water touches then (foot cells holding water) and the straddling cells
    span(I) = isempty(I) ? (NaN, NaN) : (minimum(ymin[I]), maximum(ymax[I]))
    if kf !== nothing
        V = FA.volume(r, kf)
        cross = [i for i in findall(row) if fin(V, i) > _PW_CELL_V]
        touch = [i for i in findall(foot) if fin(V, i) > _PW_CELL_V]
        sw = [i for i in findall(strad) if fin(V, i) > _PW_CELL_V]
        (c0, c1), (t0, t1) = span(cross), span(touch)
        d["cross_n_cells"] = length(cross)
        d["cross_frac_row"] = length(cross) / max(count(row), 1)
        d["cross_y_min"] = c0; d["cross_y_max"] = c1
        d["cross_y_mean"] = isempty(cross) ? NaN : sum(yc[cross]) / length(cross)
        d["cross_span_frac"] = (c1 - c0) / W
        d["touch_n_cells"] = length(touch)
        d["touch_y_min"] = t0; d["touch_y_max"] = t1
        d["cross_n_straddle_wet"] = length(sw)
        # crossing cells per wet foot cell: ≈ 1 if water crosses wherever it touches the wall
        d["cross_per_touch"] = isempty(touch) ? NaN : length(cross) / length(touch)
        d["v_down_at_first_down"] = v_down[kf]
        # the highest water level at the wall then: foot cells and the cells the wall crosses
        η = FA.wse(r, kf); at = vcat(touch, sw)
        d["wse_wall_at_first_down"] = isempty(at) ? NaN : maximum(η[at])
    end
    return d
end

function _planar_wall_metrics(ctx)
    out = Dict{String,Dict{String,Float64}}()
    for (wall, geom) in _PW_WALLS
        p = _planar_params(joinpath(ctx.dir, wall))
        g = _pw_design(geom)
        out["design_$wall"] = Dict("crest_z" => g.crest, "z_foot" => g.z_foot, "wse_end" => g.wse_end,
                                   "freeboard" => g.freeboard,
                                   "crest_z_dem" => Float64(p["emb_crest_z"]))
        for (name, r) in ctx.runs
            ctx.jobs[name].tags[:wall] == wall || continue
            out[name] = _pw_run_metrics(r, p)
        end
    end
    return out
end

# ── checks ──
_pw(m, job, key) = get(get(m, job, Dict{String,Float64}()), key, NaN)
_pw_r(x; d = 3) = round(x; digits = d)
_pw_runs(m) = sort([j for j in keys(m) if !startswith(j, "design_")])
_pw_fmt_t(t) = isfinite(t) ? "t = $(round(Int, t)) s" : "never"
"The outcome of one run in words: downstream volume, first crossing, levels, where."
function _pw_detail(m, j)
    s = "$(_pw_r(_pw(m, j, "v_down_end"); d = 1)) m³ downstream at the end " *
        "($(_pw_r(100 * _pw(m, j, "frac_down_end"); d = 2))% of the inflow), first " * _pw_fmt_t(_pw(m, j, "t_first_down"))
    isfinite(_pw(m, j, "t_first_down")) || return s * "; upstream WSE at most $(_pw_r(_pw(m, j, "wse_up_max"))) m (crest $(_pw(m, j, "crest_z")))"
    return s * " with the water at the wall no higher than $(_pw_r(_pw(m, j, "wse_wall_at_first_down"))) m (crest $(_pw(m, j, "crest_z")), " *
           "upstream WSE at the end $(_pw_r(_pw(m, j, "wse_up_end")))); crossing cells then: $(Int(_pw(m, j, "cross_n_cells"))) of " *
           "$(Int(_pw(m, j, "n_down_row"))) in the downstream row, y = $(_pw_r(_pw(m, j, "cross_y_min"); d = 0))–" *
           "$(_pw_r(_pw(m, j, "cross_y_max"); d = 0)) m (upstream water touches the wall over y = " *
           "$(_pw_r(_pw(m, j, "touch_y_min"); d = 0))–$(_pw_r(_pw(m, j, "touch_y_max"); d = 0)) m); " *
           "row wet at the end $(_pw_r(100 * _pw(m, j, "row_wet_frac_end"); d = 0))%"
end

# Checks of the correct behaviour that the model fails (matched by name
# prefix): none since 0.7.0. Until then SGS conveyed through each face from
# the DEM profile on the face line only, so a barrier inside a cell was not
# seen: both thin-wall SGS runs passed 58–64% of the inflow, from the moment
# the water reached the wall. The activation levels of 0.7.0 (connectivity,
# HYDRAULICS §9.3) close the faces behind the wall up to its crest.
const _PW_XFAIL = Pair{String,String}[]
function _pw_mark_xfail(chs)
    map(chs) do ch
        i = findfirst(p -> startswith(ch.name, first(p)), _PW_XFAIL)
        i === nothing ? ch : Check(ch.name, ch.test, ch.detail, last(_PW_XFAIL[i]))
    end
end

function _pw_checks_all()
    walls = first.(_PW_WALLS)
    chs = Check[
        Check("inflow point = the generator's injection point; the DEMs have no notch",
              m -> all(_pw(m, j, "inj_point_mismatch_deg") < 1e-9 && _pw(m, j, "notch_width_m") == 0 for j in _pw_runs(m))),
        Check("mass conserved in every run (|error| < 1e-9 of the volume added)",
              m -> length(_pw_runs(m)) == 8 && all(_pw(m, j, "mass_balance_rel") < 1e-9 for j in _pw_runs(m));
              detail = m -> "worst " * string(maximum(_pw(m, j, "mass_balance_rel") for j in _pw_runs(m)))),
        Check("closed domain stores the hydrograph volume, $(round(Int, _PW_Q * (_PW_T - 300))) m³ (±0.1%), in every run",
              m -> all(abs(_pw(m, j, "stored_m3") / (_PW_Q * (_pw(m, j, "t_end") - 300)) - 1) < 1e-3 for j in _pw_runs(m))),
        Check("design: the whole inflow held as a level pool against the wall stays ≥ $(_PW_MARGIN) m below the crest " *
              "(wedge storage), and the DEM's crest is the design crest",
              m -> all(_pw(m, "design_$w", "freeboard") >= _PW_MARGIN &&
                       abs(_pw(m, "design_$w", "crest_z_dem") - _pw(m, "design_$w", "crest_z")) < 1e-9 for w in walls);
              detail = m -> join(("$w: level $(_pw_r(_pw(m, "design_$w", "wse_end"))) m, crest $(_pw(m, "design_$w", "crest_z")) m"
                                  for w in walls), "; ")),
        Check("the upstream water level at the wall's foot stays below the crest in every run, and no higher than " *
              "the level-pool level (+0.05 m)",
              m -> all(_pw(m, j, "wse_up_max") < _pw(m, j, "crest_z") &&
                       _pw(m, j, "wse_up_max") <= _pw(m, "design_" * last(split(j, '_')), "wse_end") + 0.05 for j in _pw_runs(m));
              detail = m -> join(("$j $(_pw_r(_pw(m, j, "wse_up_max")))" for j in _pw_runs(m)), ", ") * " m"),
    ]
    for (wall, geom) in _PW_WALLS, (mesh, meshname) in _PW_MESHES, (model, modelname) in _PW_MODELS
        j = _pw_job(mesh, model, wall)
        label = "$modelname, $meshname, $wall wall ($(round(Int, geom.width)) m)"
        if model == "std" && wall == "thin"
            # cell-mean beds: a wall narrower than a cell is a low ridge in the cell means,
            # overtopped early — the limit of cell-mean beds that SGS exists to remove, not
            # the defect this case measures. Checked for what it is.
            push!(chs, Check("$label: cell-mean beds see the wall as a ridge more than 1 m below its crest; no water " *
                             "downstream until the upstream WSE reaches that ridge (the lowest cell-mean route, within 0.05 m)",
                             m -> _pw(m, j, "z_bneck_mean") < _pw(m, j, "crest_z") - 1.0 &&
                                  (!isfinite(_pw(m, j, "t_first_down")) ||
                                   _pw(m, j, "wse_up_at_first_down") >= _pw(m, j, "z_bneck_mean") - 0.05);
                             detail = m -> "cell-mean ridge $(_pw_r(_pw(m, j, "z_bneck_mean"))) m (ground at the wall's foot " *
                                           "$(_pw_r(_pw(m, j, "z_foot"))) m); " * _pw_detail(m, j)))
        else
            push!(chs, Check("$label: no water downstream of the wall " *
                             "(< $(_PW_DOWN_THRESHOLD) m³ at the end) — the upstream level never reaches the crest",
                             m -> _pw(m, j, "v_down_end") < _PW_DOWN_THRESHOLD;
                             detail = m -> _pw_detail(m, j)))
        end
    end
    return chs
end
_pw_checks() = _pw_mark_xfail(_pw_checks_all())

let
    hydro = "t_s,Q_m3s\n0,0\n600,$(_PW_Q)\n1000000,$(_PW_Q)\n"
    gen = "{repo}/test/planar_embankment/generate_planar_embankment_dem.py"
    aoi(w) = "{case}/$w/planar_embankment_aoi.geojson"
    dem(w) = "{case}/$w/planar_embankment_dem.tif"
    fm(model) = model == "std" ? "standard" : model
    meshfile(mesh, model, w) = "{case}/$w/$(mesh)_$(fm(model)).parquet"
    prep = Prep[]
    for (w, geom) in _PW_WALLS
        push!(prep, Prep("dem_$w", :python,
                         [gen, "--out", "{case}/$w",
                          "--domain-km", string(_PW_DOMAIN_KM[1]), string(_PW_DOMAIN_KM[2]),
                          "--pixel-size", string(_PW_PIXEL), "--slope-pct", string(_PW_SLOPE_PCT),
                          "--emb-height", string(_PW_HEIGHT), "--emb-width", string(geom.width),
                          "--emb-pos", string(geom.x / (_PW_DOMAIN_KM[1] * 1000.0)), "--notch-width", "0"];
                         outputs = [dem(w), "{case}/$w/planar_embankment_params.json", aoi(w)]))
        for (model, _) in _PW_MODELS
            push!(prep, Prep("a5_$(model)_$w", :floodmodel,
                             ["--meshgen", aoi(w), "--meshres", "18", "--dem", dem(w),
                              "--meshout", meshfile("a5", model, w), "--flow-model", fm(model), "--mesh-only"];
                             outputs = [meshfile("a5", model, w)]),
                        Prep("cart_$(model)_$w", :floodmodel,
                             ["--meshgen-cartesian", aoi(w), "--dx", "20", "--cartesian-lat-ref", string(_PW_CENTRE_LAT),
                              "--dem", dem(w), "--meshout", meshfile("cart", model, w),
                              "--flow-model", fm(model), "--mesh-only"];
                             outputs = [meshfile("cart", model, w)]))
        end
    end
    push!(prep, Prep("hydrograph", :write, ["{case}/inflow_q10.csv", hydro]; outputs = ["{case}/inflow_q10.csv"]))
    run_args(mesh, model, w) =
        ["--meshload", meshfile(mesh, model, w), "--flow-model", fm(model),
         "--inflow-point", "$(_PW_LAT),$(_PW_LON),{case}/inflow_q10.csv",
         "--closed-boundaries", "--sim-duration", string(round(Int, _PW_T)),
         "--dt-max", "2", "--output-interval", string(round(Int, _PW_DT_OUT))]
    jobs = [Job(_pw_job(mesh, model, w), run_args(mesh, model, w); threads = 2,
                mesh = mesh, flow_model = fm(model), wall = w)
            for (w, _) in _PW_WALLS for (mesh, _) in _PW_MESHES for (model, _) in _PW_MODELS]
    register!(Case(
        name = "planar_wall", tier = :standard,
        title = "A wall with no gap across a planar slope: no water may pass below the crest (SGS hops it)",
        doc = "docs/TEST_CASES.md §2d",
        description = """
        The planar_embankment plane (2 × 0.5 km, 0.5% west→east, 2 m DEM) with a
        3 m-high wall across it and no notch: `wide`, 40 m at x = 980–1020 m
        (crest 8.0 m; wider than a res-18 A5 cell, ~22 m) and `thin`, 4 m at
        x = 988–992 m (crest 8.05 m; 2 pixels, narrower than any cell and inside
        the cells of both meshes), each with its own DEM and meshes. 10 m³/s
        enters at x = 200 m on the centre line (10-minute ramp), closed
        boundaries, 6 h, output every 2 min: 213 000 m³, which as a level pool
        against the wall stands 2.06 m above the ground at its upstream face
        (wedge storage), 0.84 m (wide) and 0.93 m (thin) below the crest. So in
        a correct model no water passes. Runs: SGS and standard flow on A5
        res 18 and Cartesian 20 m, for each wall (8). Per run: the volume east
        of the wall through time (cells entirely east of it), its final value
        and fraction of the inflow, the time it first exceeds 1 m³, the
        upstream water level at the wall's foot then and at the end, and
        where the water crossed — the cells of the first row east of the wall
        that hold water at that time (count, y-range), against the stretch of
        the wall the upstream water touches. The checks state the correct
        behaviour; the SGS thin-wall ones fail (SGS conveys through a face from
        the DEM profile on the face line only: a barrier inside a cell is not
        seen) and are expected failures until that is fixed; standard flow on
        the thin wall is checked against the cell-mean ridge it sees. The wide
        wall holds in every run. Preparation ~10 min (2 DEMs, 8 meshes); runs
        ~1 min each on 2 threads.""",
        prep = prep, jobs = jobs,
        analyse = _planar_wall_metrics,
        checks = _pw_checks(),
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
