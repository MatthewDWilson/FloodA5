# cases/registry/planar_embankment.jl — SGS sub-cell routing through a notch in
# an embankment (docs/TEST_CASES.md §2a). Replaces
# test/planar_embankment/test_planar_embankment.jl (broken by design: shared
# DEM regenerated without the embankment, an upstream WSE taken over dry cells,
# and a source too small to reach the notch; PROJECT_STATE.md "Current plan" 3).
#
# Terrain (generate_planar_embankment_dem.py, own DEM in {case}): a 2 × 0.5 km
# plane sloping 0.5% west→east (z = 10 m at the west edge, 0 at the east), a
# 40 m-wide embankment across it at x = 1000 m (natural ground there 5.0 m),
# flat top at the crest z = 7.0 m (2 m high), cut by a 4 m notch on the centre
# line whose floor is the natural ground at the embankment centre line, the
# sill z = 5.0 m. On the 2 m DEM the notch is exactly 2 pixels (the generator
# takes pixels with |y − y_c| ≤ w/2) and the embankment 20 pixels.
#
# A res-18 A5 cell (~490 m², ~22 m) is 5× wider than the notch. Standard flow
# uses cell-mean elevations: the notch only lowers the mean of the cells it
# crosses by ≈ 2 m · 4 m/22 m ≈ 0.4 m, so water is held back until the upstream
# water surface nears the crest. SGS keeps each cell's sub-grid elevation
# distribution and the minimum along every shared edge, so it sees the notch
# and passes water at the sill level, ~2 m below the crest.

const _PE_DOMAIN_KM = (2.0, 0.5)
const _PE_SLOPE_PCT = 0.5
const _PE_HEIGHT = 2.0
const _PE_WIDTH = 40.0
const _PE_NOTCH = 4.0
const _PE_PIXEL = 2.0
const _PE_LAT, _PE_LON = _planar_point(; domain_km = _PE_DOMAIN_KM)   # x = 200 m, centre line
#
# Source. Upstream storage to a water level η, against the embankment's west
# face (x_w = 980 m, ground there z_w = 5.1 m), is a wedge
#     V(η) = ½ (η − z_w)² / S · W  = 5·10⁴ (η − z_w)² m³   (S = 0.005, W = 500 m),
# 1.8·10⁵ m³ to the crest (η = 7 m) and ≈ 1.2·10⁵ m³ to the cell-mean crest
# of the standard mesh (≈ 6.6 m). 10 m³/s (after a 10-minute ramp) therefore
# overtops the standard mesh after ≈ 3.5–5 h, inside the 6 h run. Through the
# 4 m notch a broad-crested weir passes 10 m³/s at a head of
# (10 / (1.7 · 4))^⅔ ≈ 1.3 m above the sill, so SGS should settle ≈ 0.7 m below
# the crest; its first downstream water comes when the plume reaches the notch.
const _PE_Q = 10.0
const _PE_T = 21600.0
const _PE_DT_OUT = 120.0
const _PE_DOWN_THRESHOLD = 1.0          # m³ downstream = "water has passed"
const _PE_WET = 0.01                    # m: a cell is wet above this depth

"Cells' x extent (m from the DEM's west edge, generator convention) from their polygons."
function _pe_xrange(mesh, p)
    m_lon = 111320.0 * cos(deg2rad(Float64(p["centre_lat"])))
    lon0 = Float64(p["lon0"])
    xmin = [minimum(view(P, 1, :)) for P in mesh.polygons]
    xmax = [maximum(view(P, 1, :)) for P in mesh.polygons]
    return (xmin .- lon0) .* m_lon, (xmax .- lon0) .* m_lon
end

"""
    _pe_bottleneck(z, nb, up, down) -> Float64

The lowest level at which the `up` cells connect to the `down` cells through
neighbouring cells, i.e. min over paths of the max of `z` along the path
(cells added in order of `z`, union–find). For cell-mean elevations this is
the crest standard flow sees; for `sgs_z_min` the lowest sub-grid route.
"""
function _pe_bottleneck(z::AbstractVector, nb::AbstractMatrix, up::AbstractVector{Bool},
                        down::AbstractVector{Bool})
    n = length(z)
    parent = collect(1:n)
    hasu = copy(up); hasd = copy(down)
    added = falses(n)
    find(i) = (while parent[i] != i; parent[i] = parent[parent[i]]; i = parent[i]; end; i)
    for i in sortperm(z)
        added[i] = true
        for s in axes(nb, 1)
            j = Int(nb[s, i]); (j == 0 || !added[j]) && continue
            a, b = find(i), find(j)
            a == b && continue
            parent[b] = a
            hasu[a] |= hasu[b]; hasd[a] |= hasd[b]
        end
        r = find(i)
        hasu[r] && hasd[r] && return Float64(z[i])
    end
    return NaN
end

function _planar_embankment_metrics(ctx)
    p = _planar_params(ctx.dir)
    crest, sill = Float64(p["emb_crest_z"]), Float64(p["notch_sill_z"])
    xw = Float64(p["emb_x_m"]) - Float64(p["emb_width_m"]) / 2
    xe = Float64(p["emb_x_m"]) + Float64(p["emb_width_m"]) / 2
    out = Dict{String,Dict{String,Float64}}()
    series = Dict{String,Any}()
    for (name, r) in ctx.runs
        mesh = r.mesh
        xmin, xmax = _pe_xrange(mesh, p)
        cell = sqrt(sum(mesh.area) / length(mesh.area))   # mean cell size, ≈ 22 m
        up = xmax .<= xw                          # entirely west of the embankment
        down = xmin .>= xe                        # entirely east of it
        adj = up .& (xmax .>= xw - cell)          # upstream cells along its foot
        n = FA.nframes(r)
        wse_up = fill(NaN, n); v_down = zeros(n)
        for k in 1:n
            h = FA.depth(r, k); η = FA.wse(r, k); V = FA.volume(r, k)
            w = [i for i in findall(adj) if isfinite(h[i]) && h[i] > _PE_WET]
            isempty(w) || (wse_up[k] = maximum(η[w]))
            v_down[k] = sum(V[i] for i in findall(down) if isfinite(V[i]); init = 0.0)
        end
        kf = findfirst(>(_PE_DOWN_THRESHOLD), v_down)
        d = Dict{String,Float64}(
            "inj_point_mismatch_deg" => max(abs(Float64(p["injection_lat"]) - _PE_LAT),
                                            abs(Float64(p["injection_lon"]) - _PE_LON)),
            "n_up_adjacent" => count(adj), "n_down" => count(down),
            "crest_z" => crest, "sill_z" => sill,
            "z_foot" => sill + Float64(p["slope_pct"]) / 100 * Float64(p["emb_width_m"]) / 2,
            "z_bneck_mean" => _pe_bottleneck(mesh.elevation, mesh.neighbours, up, down),
            "t_first_down" => kf === nothing ? Inf : r.t[kf],
            "wse_up_at_first_down" => kf === nothing ? NaN : wse_up[kf],
            "wse_up_before_first_down" => (kf === nothing || kf == 1) ? NaN : wse_up[kf-1],
            "wse_up_max" => maximum(filter(isfinite, wse_up); init = -Inf),
            "v_down_end" => v_down[end],
            # Σ volume, not Σ depth·area: SGS depth is measured from the cell minimum
            "stored_m3" => FA.stored_volume(r, n), "t_end" => r.t[n])
        haskey(mesh.static, "sgs_z_min") &&
            (d["z_bneck_zmin"] = _pe_bottleneck(mesh.static["sgs_z_min"], mesh.neighbours, up, down))
        for th in 1:6
            k = FA.frame_at(r, 3600.0 * th; tol = 30.0)
            d["wse_up_$(th)h"] = wse_up[k]
            d["v_down_$(th)h"] = v_down[k]
        end
        out[name] = d
        series[name] = (; t = r.t, wse_up, v_down)
    end
    # SGS against standard at matched times between SGS's first downstream
    # water and standard's (before standard overtops)
    s, t = series["sgs"], series["standard"]
    t0, t1 = out["sgs"]["t_first_down"], out["standard"]["t_first_down"]
    ks = [(ka, kb) for (ka, kb) in FA.matched_frames(ctx.runs["sgs"], ctx.runs["standard"]; tol = 30.0)
          if t0 <= s.t[ka] < t1]
    out["sgs_vs_standard"] = Dict(
        "n_frames_compared" => length(ks),
        "frac_sgs_more_down" => isempty(ks) ? NaN : count(((a, b),) -> s.v_down[a] > t.v_down[b], ks) / length(ks),
        "v_down_sgs_before_std_overtops" => isempty(ks) ? NaN : s.v_down[last(ks)[1]],
        "v_down_std_before_std_overtops" => isempty(ks) ? NaN : t.v_down[last(ks)[2]])
    return out
end

_pe(m, job, key) = get(get(m, job, Dict{String,Float64}()), key, NaN)
_pe_r(x) = round(x; digits = 3)

let
    aoi, dem = "{case}/planar_embankment_aoi.geojson", "{case}/planar_embankment_dem.tif"
    hydro = "t_s,Q_m3s\n0,0\n600,$(_PE_Q)\n1000000,$(_PE_Q)\n"
    run_args(mesh, fm) = ["--meshload", "{case}/$mesh.parquet", "--flow-model", fm,
                          "--inflow-point", "$(_PE_LAT),$(_PE_LON),{case}/inflow_q10.csv",
                          "--closed-boundaries", "--sim-duration", string(round(Int, _PE_T)),
                          "--dt-max", "2", "--output-interval", string(round(Int, _PE_DT_OUT))]
    prep = [Prep("dem", :python,
                 ["{repo}/test/planar_embankment/generate_planar_embankment_dem.py", "--out", "{case}",
                  "--domain-km", string(_PE_DOMAIN_KM[1]), string(_PE_DOMAIN_KM[2]),
                  "--pixel-size", string(_PE_PIXEL), "--slope-pct", string(_PE_SLOPE_PCT),
                  "--emb-height", string(_PE_HEIGHT), "--emb-width", string(_PE_WIDTH),
                  "--notch-width", string(_PE_NOTCH)];
                 outputs = ["{case}/planar_embankment_dem.tif", "{case}/planar_embankment_params.json",
                            "{case}/planar_embankment_aoi.geojson"]),
            [Prep("a5_res18_$fm", :floodmodel,
                  ["--meshgen", aoi, "--meshres", "18", "--dem", dem,
                   "--meshout", "{case}/a5_res18_$fm.parquet", "--flow-model", fm, "--mesh-only"];
                  outputs = ["{case}/a5_res18_$fm.parquet"]) for fm in ("sgs", "standard")]...,
            Prep("hydrograph", :write, ["{case}/inflow_q10.csv", hydro]; outputs = ["{case}/inflow_q10.csv"])]
    register!(Case(
        name = "planar_embankment", tier = :standard,
        title = "Sub-cell notch in an embankment: SGS passes water at the sill, standard at the crest",
        doc = "docs/TEST_CASES.md §2a",
        description = """
        A 2 × 0.5 km plane sloping 0.5% west→east (2 m DEM) with a 40 m-wide,
        2 m-high embankment across it at x = 1000 m (crest 7.0 m), cut by a 4 m
        notch down to the natural ground (sill 5.0 m) on the centre line — a
        fifth of a res-18 A5 cell. 10 m³/s enters at x = 200 m on the centre line
        (10-minute ramp), closed boundaries, 6 h, output every 2 min, on res-18
        SGS and standard meshes. Per output time: the upstream water surface (the
        highest WSE of the wet cells along the embankment's upstream foot) and the
        volume stored east of the embankment. Standard flow (cell-mean
        elevations) should hold the water until the upstream WSE reaches its
        cell-mean crest; SGS should pass water through the notch near the sill,
        far below the crest, and so route more water downstream until standard
        overtops. z_bneck_mean / z_bneck_zmin: the lowest level at which the
        upstream cells connect to the downstream ones through cell-mean
        elevations / cell minimum elevations (a mesh property). Each run is
        10 800 steps, 30–40 s of simulation on 2 threads (2189 cells).""",
        prep = prep,
        jobs = [Job("sgs", run_args("a5_res18_sgs", "sgs"); mesh = "a5", res = 18, flow_model = "sgs"),
                Job("standard", run_args("a5_res18_standard", "standard"); mesh = "a5", res = 18, flow_model = "standard")],
        analyse = _planar_embankment_metrics,
        checks = [
            Check("inflow point = the generator's injection point (params file)",
                  m -> _pe(m, "sgs", "inj_point_mismatch_deg") < 1e-9),
            Check("mass conserved (|error| < 1e-9 of the volume added)",
                  m -> _pe(m, "sgs", "mass_balance_rel") < 1e-9 && _pe(m, "standard", "mass_balance_rel") < 1e-9;
                  detail = m -> "sgs $(_pe(m, "sgs", "mass_balance_rel")), standard $(_pe(m, "standard", "mass_balance_rel"))"),
            Check("closed domain stores the hydrograph volume, $(round(Int, _PE_Q * (_PE_T - 300))) m³ (±0.1%)",
                  m -> all(abs(_pe(m, j, "stored_m3") / (_PE_Q * (_pe(m, j, "t_end") - 300)) - 1) < 1e-3 for j in ("sgs", "standard"));
                  detail = m -> "sgs $(round(_pe(m, "sgs", "stored_m3"); digits = 1)), standard $(round(_pe(m, "standard", "stored_m3"); digits = 1)) m³"),
            Check("standard: the 2 m embankment is a cell-mean crest within 0.5 m of the true crest " *
                  "(the 4 m notch lowers the cells it crosses by ≈ 0.4 m)",
                  m -> _pe(m, "standard", "z_bneck_mean") >= _pe(m, "standard", "crest_z") - 0.5;
                  detail = m -> "cell-mean crest $(_pe_r(_pe(m, "standard", "z_bneck_mean"))) m, crest $(_pe(m, "standard", "crest_z")) m"),
            Check("standard: no water downstream until the upstream WSE reaches the cell-mean crest " *
                  "(within 0.05 m) and so is within 0.5 m of the crest; it overtops within the run",
                  m -> isfinite(_pe(m, "standard", "t_first_down")) &&
                       _pe(m, "standard", "wse_up_at_first_down") >= _pe(m, "standard", "z_bneck_mean") - 0.05 &&
                       _pe(m, "standard", "wse_up_at_first_down") >= _pe(m, "standard", "crest_z") - 0.5;
                  detail = m -> "first downstream water at t = $(_pe(m, "standard", "t_first_down")) s, upstream WSE " *
                                "$(_pe_r(_pe(m, "standard", "wse_up_at_first_down"))) m (before: " *
                                "$(_pe_r(_pe(m, "standard", "wse_up_before_first_down"))) m)"),
            Check("SGS: water passes the notch while the upstream WSE is ≥ 0.5 m below the crest " *
                  "and within 0.5 m of the sill",
                  m -> _pe(m, "sgs", "wse_up_at_first_down") <= _pe(m, "sgs", "crest_z") - 0.5 &&
                       _pe(m, "sgs", "wse_up_at_first_down") <= _pe(m, "sgs", "sill_z") + 0.5;
                  detail = m -> "first downstream water at t = $(_pe(m, "sgs", "t_first_down")) s, upstream WSE " *
                                "$(_pe_r(_pe(m, "sgs", "wse_up_at_first_down"))) m (sill $(_pe(m, "sgs", "sill_z")), " *
                                "crest $(_pe(m, "sgs", "crest_z")))"),
            Check("SGS mesh: the notch is an open route — the lowest route through cell minima is set by " *
                  "the ground at the embankment's upstream foot (within 0.05 m), not by the embankment",
                  m -> _pe(m, "sgs", "z_bneck_zmin") <= _pe(m, "sgs", "z_foot") + 0.05;
                  detail = m -> "z_bneck_zmin $(_pe_r(_pe(m, "sgs", "z_bneck_zmin"))) m, foot " *
                                "$(_pe_r(_pe(m, "sgs", "z_foot"))) m, sill $(_pe(m, "sgs", "sill_z")) m"),
            Check("SGS routes more water downstream than standard at every matched time before standard overtops",
                  m -> _pe(m, "sgs_vs_standard", "n_frames_compared") >= 10 &&
                       _pe(m, "sgs_vs_standard", "frac_sgs_more_down") == 1.0;
                  detail = m -> "$(Int(_pe(m, "sgs_vs_standard", "n_frames_compared"))) frames; just before standard " *
                                "overtops: SGS $(round(_pe(m, "sgs_vs_standard", "v_down_sgs_before_std_overtops"); digits = 0)) m³, " *
                                "standard $(round(_pe(m, "sgs_vs_standard", "v_down_std_before_std_overtops"); digits = 1)) m³"),
        ],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
