# visualisation/plots/core.jl — the Makie-free part of FloodA5Plots.
#
# House style (colours as hex strings, which Makie accepts anywhere a colour
# goes), colour-range choice, frame selection, extent classification, map
# coordinates and run sources. Pure Julia on top of FA: CI tests these without
# CairoMakie (test/test_plots.jl).

# ───────────────────────────── house style ─────────────────────────────────
# Chart ink and surfaces (light mode), as test/channel_bench/plot_channel.jl.

"Primary ink (titles, labels)."
const INK = "#0b0b0b"
"Secondary ink (tick labels, legend text)."
const INK2 = "#52514e"
"Muted ink (axes, ticks, reference lines)."
const MUTED = "#898781"
"Gridline hairline."
const GRID = "#e1e0d9"
"Chart surface (figure and axis background)."
const SURF = "#fcfcfb"

"""
Categorical series colours, in fixed order (never cycled): blue, orange,
aqua, yellow, magenta, green, violet, red. Slots 1–3 are distinguishable for
colour-vision deficiencies in any pairing; for more than three series use a
legend and a second encoding (marker or line style).
"""
const SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                "#e87ba4", "#008300", "#6250d6", "#e34948"]

"Markers paired with `SERIES` slots (second encoding for > 3 series)."
const MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle, :star5, :hexagon, :cross]

"Line styles paired with `SERIES` slots."
const LINESTYLES = [:solid, :solid, :solid, :dash, :dash, :dot, :dot, :dashdot]

"""
Mesh families' fixed identity (colour, marker, label): A5 is always blue,
Cartesian always orange, in every figure.
"""
const MESH_STYLE = Dict("a5" => (label = "A5", color = SERIES[1], marker = :circle),
                        "cart" => (label = "Cartesian", color = SERIES[2], marker = :rect))

"""
    series_color(i) -> hex

Categorical colour of series `i` (1 … 8). An error past 8: fold the extra
series into small multiples instead of inventing a hue.
"""
function series_color(i::Integer)
    1 <= i <= length(SERIES) ||
        error("series_color: $i series — at most $(length(SERIES)) (use small multiples)")
    return SERIES[i]
end

"Line style per flow model (`run_styles(...; style = :model)`)."
const MODEL_LINESTYLE = Dict("standard" => :solid, "swe" => :dash, "sgs" => :dot)
"Display names of the flow models."
const MODEL_LABEL = Dict("standard" => "local-inertial", "swe" => "SWE", "sgs" => "SGS")

"""
    mesh_family(run) -> "a5" | "cart" | ""

The run's mesh family from its id scheme ("a5-v0", "a5-v1" → "a5";
"cartesian" → "cart").
"""
function mesh_family(r::FA.Run)
    s = r.mesh.id_scheme
    startswith(s, "a5") && return "a5"
    (s == "cartesian" || haskey(r.mesh.static, "cell_size_m")) && return "cart"
    return ""
end

"""
    run_styles(runs; style = :auto) -> Vector{(color, linestyle, marker)}

Line styles for several runs on one axis.

- `:series`: run i gets `SERIES[i]`, `LINESTYLES[i]`, `MARKERS[i]` (at most 8).
- `:model`: colour by mesh family (`MESH_STYLE`: A5 blue, Cartesian orange),
  line style by flow model (`MODEL_LINESTYLE`: local-inertial solid, SWE
  dashed, SGS dotted) — the same run type looks the same in every figure.
- `:auto`: `:model` when every run's (mesh family, flow model) is known and
  distinct, else `:series`.
"""
function run_styles(runs::AbstractVector; style::Symbol = :auto)
    keys_ = [(mesh_family(r), r.flow_model) for r in runs]
    model_ok = all(k -> haskey(MESH_STYLE, k[1]) && haskey(MODEL_LINESTYLE, k[2]), keys_) &&
               allunique(keys_)
    style === :auto && (style = model_ok ? :model : :series)
    if style === :model
        model_ok || error("run_styles: style = :model needs distinct runs of known mesh family and flow model")
        return [(color = MESH_STYLE[f].color, linestyle = MODEL_LINESTYLE[m], marker = MESH_STYLE[f].marker)
                for (f, m) in keys_]
    end
    style === :series || error("style must be :auto, :model or :series")
    return [(color = series_color(i), linestyle = LINESTYLES[i], marker = MARKERS[i]) for i in eachindex(runs)]
end

"Sequential blue ramp (light → dark): water depth, speed, any magnitude."
const SEQUENTIAL = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
                    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

"""
Diverging ramp for differences (run − reference): red (negative: shallower
than the reference) ← grey midpoint → blue (positive: deeper).
"""
const DIVERGING = ["#7d1a1c", "#b02e2e", "#e34948", "#ee8c88", "#f6c6c2", "#f0efec",
                   "#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]

"Basemap ramp for dry-cell elevation (low → high: mid grey → near-white)."
const TERRAIN = ["#a4a39b", "#c3c2b9", "#dcdbd4", "#eeede8"]

"Uniform basemap colour (`basemap = :grey`)."
const BASEMAP_GREY = "#e4e3dc"

"""
Extent classes (`extent_classes`): 1 hit (wet in both), 2 miss (reference
only), 3 false alarm (run only) — CVD-safe slots 1–3.
"""
const EXTENT_COLORS = [SERIES[1], SERIES[2], SERIES[3]]
const EXTENT_LABELS = ["Wet in both (hit)", "Reference only (miss)", "Run only (false alarm)"]

# ───────────────────────────── variables ───────────────────────────────────

"""
    var_spec(var) -> (label, unit, masked)

Axis/colourbar label and unit of a map variable, and whether dry cells
(depth ≤ `dry_threshold`) are masked. Variables: `:depth`, `:wse`, `:speed`,
`:volume`, `:elevation`, or any per-frame output variable name (as a Symbol
or String, e.g. `"saturation"`).
"""
function var_spec(var)
    v = Symbol(var)
    v === :depth     && return (label = "Depth", unit = "m", masked = true)
    v === :wse       && return (label = "Water surface", unit = "m", masked = true)
    v === :speed     && return (label = "Speed", unit = "m/s", masked = true)
    v === :volume    && return (label = "Volume", unit = "m³", masked = true)
    v === :elevation && return (label = "Elevation", unit = "m", masked = false)
    return (label = String(v), unit = "", masked = false)
end

"Colourbar/axis label: \"Depth (m)\"."
var_label(var) = (s = var_spec(var); isempty(s.unit) ? s.label : "$(s.label) ($(s.unit))")

"""
    frame_values(run, var, k) -> Vector{Float64}

Per-cell values of `var` (see `var_spec`) at frame `k`: `FA.depth`, `FA.wse`
(SGS-aware), `FA.speed`, `FA.volume`, the mesh elevation, or `FA.field`.
"""
function frame_values(r::FA.Run, var, k::Integer)
    v = Symbol(var)
    v === :depth     && return FA.depth(r, k)
    v === :wse       && return FA.wse(r, k)
    v === :speed     && return FA.speed(r, k)
    v === :volume    && return FA.volume(r, k)
    v === :elevation && return copy(r.mesh.elevation)
    return FA.field(r, String(v), k)
end

"Run label for titles: the label without a trailing `.h5`."
run_label(r::FA.Run) = endswith(r.label, ".h5") ? r.label[1:end-3] : r.label

# ───────────────────────────── colour ranges ───────────────────────────────

const _NICE = (1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)

# s·10^e with the division for negative e, so 3·10⁻¹ is exactly 0.3
_scale10(s, e) = e >= 0 ? s * 10.0^e : s / 10.0^(-e)

"""
    nice_ceil(x) -> Float64

The smallest "round" number ≥ x of the form {1, 1.2, 1.5, 2, 2.5, 3, 4, 5,
6, 8} × 10ⁿ (x > 0); 0 for x ≤ 0. `nice_ceil(3.7) == 4.0`.
"""
function nice_ceil(x::Real)
    (x > 0 && isfinite(x)) || return 0.0
    e = floor(log10(x)); m = _scale10(x, -e)
    for s in _NICE
        s >= m * (1 - 1e-12) && return _scale10(s, e)
    end
    return _scale10(1.0, e + 1)
end

"""
    robust_range(x; mask = nothing, lo = nothing, q = 0.98, symmetric = false,
                 nice = true) -> (lo, hi)

Colour range of `x` robust to outliers: over the finite values (where `mask`
is true), the `q` quantile as the top (rounded up with `nice_ceil` if `nice`).

- `lo = 0.0`: a magnitude (depth, speed) — (0, top).
- `symmetric = true`: a difference — (−m, m), m the `q` quantile of |x|.
- otherwise: the (1 − q, q) quantiles (elevation, WSE), not rounded.

No finite value → (0, 1); a zero-width range is widened so it is valid.

    robust_range(h; mask = h .> 0.01, lo = 0.0)       # depth
    robust_range(dh; symmetric = true)                # Δh
"""
function robust_range(x::AbstractVector; mask = nothing, lo = nothing, q::Real = 0.98,
                      symmetric::Bool = false, nice::Bool = true)
    v = [Float64(x[i]) for i in eachindex(x) if isfinite(x[i]) && (mask === nothing || mask[i])]
    isempty(v) && return (0.0, 1.0)
    if symmetric
        m = quantile(abs.(v), q)
        m = nice ? nice_ceil(m) : m
        m > 0 || (m = max(maximum(abs, v), 1e-3))
        return (-m, m)
    elseif lo !== nothing
        top = quantile(v, q)
        top = nice ? nice_ceil(top - lo) + lo : top
        top > lo || (top = lo + max(maximum(v) - lo, 1e-3))
        return (Float64(lo), top)
    end
    a, b = quantile(v, 1 - q), quantile(v, q)
    b > a || (a -= 0.5; b += 0.5)
    return (a, b)
end

"""
    default_range(var, values; mask = nothing) -> (lo, hi)

`robust_range` with the variable's convention: magnitudes (depth, speed,
volume, saturation) from 0; elevation and WSE between quantiles.
"""
default_range(var, values::AbstractVector; mask = nothing) =
    Symbol(var) in (:elevation, :wse) ? robust_range(values; mask, nice = false) :
                                        robust_range(values; mask, lo = 0.0)

"""
    color_index(x, (lo, hi), n) -> i

Index (1 … n) of value `x` on an n-colour ramp spanning `lo … hi`, clamped
at both ends (NaN → 0, i.e. not drawn).
"""
function color_index(x::Real, range::Tuple{<:Real,<:Real}, n::Integer)
    isfinite(x) || return 0
    lo, hi = range
    f = hi > lo ? (x - lo) / (hi - lo) : 0.5
    return clamp(1 + floor(Int, f * n), 1, n)
end

# ───────────────────────────── time and frames ─────────────────────────────

"""
    hhmm(t) -> "hh:mm"

Simulation time `t` (s) rounded to the minute: `hhmm(7210.0) == "02:00"`.
"""
function hhmm(t::Real)
    m = round(Int, t / 60)
    return string(lpad(div(m, 60), 2, '0'), ":", lpad(rem(m, 60), 2, '0'))
end

"`\"t = 02:00\"` (see `hhmm`)."
time_label(t::Real) = "t = " * hhmm(t)

"""
    time_axis(tmax) -> (seconds_per_unit, label)

Unit of a time axis running to `tmax` s: minutes up to 2 h, else hours.
"""
time_axis(tmax::Real) = tmax <= 7200 ? (60.0, "Time (min)") : (3600.0, "Time (h)")

"""
    frame_tol(run) -> s

Default tolerance for matching a requested time to a frame: 60 % of the
median frame interval (so the nearest frame is accepted, an off-range time
is not), at least 60 s.
"""
frame_tol(r::FA.Run) = FA.nframes(r) < 2 ? Inf : max(60.0, 0.6 * median(diff(r.t)))

"""
    frame_grid(runs, times; tol = nothing) -> Matrix{Int}

Frame indices, `length(runs) × length(times)`: for each run the frame
nearest each time (`FA.frame_at`), an error if it is more than `tol` s away
(default `frame_tol` of that run).
"""
function frame_grid(runs::AbstractVector, times::AbstractVector; tol = nothing)
    K = Matrix{Int}(undef, length(runs), length(times))
    for (i, r) in enumerate(runs), (j, t) in enumerate(times)
        K[i, j] = FA.frame_at(r, t; tol = tol === nothing ? frame_tol(r) : tol)
    end
    return K
end

"""
    matched_pair(ref, run; t = nothing, k = nothing, tol = nothing) -> (k_ref, k_run)

Frames of `ref` and `run` at one simulation time: `t` if given (each run's
nearest frame, within `tol`); else frame `k` of `run` and the reference
frame nearest its time; else the last time both runs reach (as `FA.compare`).
"""
function matched_pair(ref::FA.Run, run::FA.Run; t = nothing, k = nothing, tol = nothing)
    tl(r) = tol === nothing ? frame_tol(r) : tol
    if t !== nothing
        return FA.frame_at(ref, t; tol = tl(ref)), FA.frame_at(run, t; tol = tl(run))
    elseif k !== nothing
        return FA.frame_at(ref, run.t[k]; tol = tl(ref)), Int(k)
    end
    tc = min(ref.t[end], run.t[end])
    return FA.frame_at(ref, tc), FA.frame_at(run, tc)
end

"Check that two runs are on the same mesh in the same cell order."
function same_mesh(a::FA.Run, b::FA.Run)
    a.mesh.ids == b.mesh.ids ||
        error("$(run_label(a)) and $(run_label(b)) are not on the same mesh (cell ids differ)")
    return true
end

# ───────────────────────────── comparisons ─────────────────────────────────

"""
    extent_classes(h_ref, h; threshold = 0.1) -> Vector{Int8}

Per cell: 0 dry in both, 1 hit (h_ref > threshold and h > threshold),
2 miss (reference only), 3 false alarm (run only). NaN counts as dry, as in
`FA.compare_fields`, whose `hits`/`misses`/`false_alarms` these counts equal.
"""
function extent_classes(h_ref::AbstractVector, h::AbstractVector; threshold::Real = 0.1)
    length(h_ref) == length(h) || error("extent_classes: lengths differ")
    c = zeros(Int8, length(h))
    for i in eachindex(h, h_ref)
        a = h_ref[i] > threshold; b = h[i] > threshold
        c[i] = a & b ? 1 : a ? 2 : b ? 3 : 0
    end
    return c
end

"""
    extent_scores(classes) -> (hits, misses, false_alarms, csi)

Counts of each class and the critical success index
hits / (hits + misses + false alarms) (0 if all are dry), as `FA.compare_fields`.
"""
function extent_scores(c::AbstractVector{<:Integer})
    hits = count(==(1), c); miss = count(==(2), c); fa = count(==(3), c)
    return (hits = hits, misses = miss, false_alarms = fa, csi = hits / max(hits + miss + fa, 1))
end

"""
    difference_values(v_ref, v, h_ref, h; wet = 0.01) -> Vector

v − v_ref on the cells wet (depth > `wet`) in either run, NaN elsewhere and
where either value is not finite. For depth, `v = h`.
"""
function difference_values(v_ref::AbstractVector, v::AbstractVector,
                           h_ref::AbstractVector, h::AbstractVector; wet::Real = 0.01)
    d = fill(NaN, length(v))
    for i in eachindex(d)
        (h_ref[i] > wet || h[i] > wet) && isfinite(v[i]) && isfinite(v_ref[i]) &&
            (d[i] = v[i] - v_ref[i])
    end
    return d
end

# ───────────────────────────── map coordinates ─────────────────────────────

"""
    MapCoords(coords, lon0, lat0)

A map projection: `coords = :lonlat` (degrees, drawn with x stretched by
cos(lat0) so a degree of longitude has its true length) or `:local`
(equirectangular metres about (lon0, lat0), `FA.local_xy`).
"""
struct MapCoords
    coords::Symbol
    lon0::Float64
    lat0::Float64
    function MapCoords(coords::Symbol, lon0::Real, lat0::Real)
        coords in (:lonlat, :local) || error("coords must be :lonlat or :local, not :$coords")
        new(coords, lon0, lat0)
    end
end

"""
    map_coords(run_or_mesh; coords = :lonlat, origin = :centre) -> MapCoords

`origin`: `:centre` (the centre of the cell centres' bounding box),
`:source` (the run's first source point, `run_sources`; the centre if it has
none) or a `(lon, lat)` tuple. For maps of several runs, pass one
`MapCoords` to all of them so they share an origin.
"""
function map_coords(m::FA.MeshGeom; coords::Symbol = :lonlat, origin = :centre, sources = nothing)
    lon0, lat0 = if origin isa Tuple || origin isa AbstractVector
        Float64(origin[1]), Float64(origin[2])
    elseif origin === :source && sources !== nothing && !isempty(sources)
        sources[1].lon, sources[1].lat
    elseif origin === :centre || origin === :source
        (minimum(m.lon) + maximum(m.lon)) / 2, (minimum(m.lat) + maximum(m.lat)) / 2
    else
        error("origin must be :centre, :source or (lon, lat)")
    end
    return MapCoords(coords, lon0, lat0)
end
map_coords(r::FA.Run; coords::Symbol = :lonlat, origin = :centre) =
    map_coords(r.mesh; coords, origin, sources = origin === :source ? run_sources(r) : nothing)
map_coords(mc::MapCoords; kw...) = mc

"""
    project(mc, lon, lat) -> (x, y)

Map coordinates of points (scalars or vectors): degrees unchanged for
`:lonlat`, metres for `:local`.
"""
function project(mc::MapCoords, lon, lat)
    mc.coords === :lonlat && return (lon, lat)
    x = deg2rad.(lon .- mc.lon0) .* FA.R_EARTH .* cosd(mc.lat0)
    y = deg2rad.(lat .- mc.lat0) .* FA.R_EARTH
    return x, y
end

"""
    project_polygons(mc, mesh; cells = eachindex(mesh.polygons)) -> Vector{Matrix}

The cells' polygons (2 × nv, closed rings) in map coordinates.
"""
function project_polygons(mc::MapCoords, m::FA.MeshGeom; cells = eachindex(m.polygons))
    isempty(m.polygons) && error("the mesh has no polygons (format-1 run: load it with parquet = MESH)")
    out = Vector{Matrix{Float64}}(undef, length(cells))
    for (j, i) in enumerate(cells)
        P = m.polygons[i]
        x, y = project(mc, view(P, 1, :), view(P, 2, :))
        out[j] = permutedims(hcat(x, y))
    end
    return out
end

"""
    axis_labels(mc) -> (xlabel, ylabel)
"""
axis_labels(mc::MapCoords) = mc.coords === :lonlat ? ("Longitude (°)", "Latitude (°)") :
                                                     ("x (m)", "y (m)")

"""
    map_extent(mc, polys; pad = 0.01) -> (xlims, ylims, aspect)

Axis limits around projected polygons (padded by `pad` of the larger side)
and the axis-box aspect ratio (width / height) that draws them to scale:
for `:lonlat` the x extent counts cos(lat0) times.
"""
function map_extent(mc::MapCoords, polys::AbstractVector; pad::Real = 0.01)
    x0 = y0 = Inf; x1 = y1 = -Inf
    for P in polys, j in axes(P, 2)
        x0 = min(x0, P[1, j]); x1 = max(x1, P[1, j])
        y0 = min(y0, P[2, j]); y1 = max(y1, P[2, j])
    end
    sx = mc.coords === :lonlat ? cosd(mc.lat0) : 1.0
    w = (x1 - x0) * sx; h = y1 - y0
    p = pad * max(w, h)
    xl = (x0 - p / sx, x1 + p / sx); yl = (y0 - p, y1 + p)
    return xl, yl, (w + 2p) / (h + 2p)
end

# ───────────────────────────── points and sources ──────────────────────────

"""
    run_sources(run) -> Vector{(lon, lat, label)}

Where the run adds water, from its record (format 2): `--inflow-point`
(label: the hydrograph's label), `--injection-point`, `--rainpoint`, and
the `P` points of an `--inflow-bci` file if that file still exists at its
recorded path. Empty for format-1 files.
"""
function run_sources(r::FA.Run)
    out = NamedTuple{(:lon, :lat, :label),Tuple{Float64,Float64,String}}[]
    o = get(r.info, "options", nothing)
    o isa AbstractDict || return out
    num(x) = x isa Number ? Float64(x) : parse(Float64, strip(String(x)))
    for s in something(get(o, "inflow_points", nothing), [])          # lat, lon, file[, label]
        lab = length(s) >= 4 ? String(s[4]) : splitext(basename(String(s[3])))[1]
        push!(out, (lon = num(s[2]), lat = num(s[1]), label = lab))
    end
    for s in something(get(o, "injection_points", nothing), [])       # lon, lat, rate
        push!(out, (lon = num(s[1]), lat = num(s[2]), label = "injection"))
    end
    for s in something(get(o, "rainpoints", nothing), [])             # lon, lat, mm/h
        push!(out, (lon = num(s[1]), lat = num(s[2]), label = "rain point"))
    end
    bci = get(o, "inflow_bci", nothing)
    if bci isa AbstractString && isfile(bci)
        append!(out, parse_bci_points(read(bci, String)))
    end
    return out
end

"""
    parse_bci_points(text) -> Vector{(lon, lat, label)}

The point (`P x y TYPE [name]`) lines of a LISFLOOD-FP `.bci` file in WGS84.
"""
function parse_bci_points(text::AbstractString)
    out = NamedTuple{(:lon, :lat, :label),Tuple{Float64,Float64,String}}[]
    for l in eachline(IOBuffer(text))
        w = split(strip(l))
        (length(w) >= 3 && w[1] == "P") || continue
        x = tryparse(Float64, w[2]); y = tryparse(Float64, w[3])
        (x === nothing || y === nothing) && continue
        push!(out, (lon = x, lat = y, label = length(w) >= 5 ? String(w[5]) : ""))
    end
    return out
end

"""
    point_cells(mesh, points) -> Vector{Int}

Cell indices for `points`: integers are cell indices already; `(lon, lat)`
tuples (or named tuples with `lon`, `lat`) are the nearest cell
(`FA.source_cell`).
"""
point_cells(m::FA.MeshGeom, points) = [_cell(m, p) for p in points]
_cell(m::FA.MeshGeom, p::Integer) = (1 <= p <= length(m) || error("cell $p out of range"); Int(p))
_cell(m::FA.MeshGeom, p::NamedTuple) = FA.source_cell(m, p.lon, p.lat)
_cell(m::FA.MeshGeom, p) = FA.source_cell(m, p[1], p[2])

# ───────────────────────────── storage and mass balance ────────────────────

"""
    storage_series(run) -> Vector{Float64}

`FA.storage` at every frame (m³): Σ depth·area, or Σ volume for SGS runs.
"""
storage_series(r::FA.Run) = [FA.storage(r, k) for k in 1:FA.nframes(r)]

"""
    mass_balance_series(run) -> NamedTuple

From the per-frame `volume` variable and the end-of-run summary:
`t` (s), `stored` (Σ volume per frame, m³), `t_mid` and `rate` (the net
inflow dS/dt between frames, m³/s), and from the summary (NaN if absent)
`added`, `removed` (m³), `net = added − removed`, `stored_end` (the model's
final domain volume) and `error_rel` (its conservation error relative to
`added`).
"""
function mass_balance_series(r::FA.Run)
    S = [FA.stored_volume(r, k) for k in 1:FA.nframes(r)]
    t = r.t
    tm = length(t) > 1 ? (t[1:end-1] .+ t[2:end]) ./ 2 : Float64[]
    rate = length(t) > 1 ? diff(S) ./ diff(t) : Float64[]
    s = get(r.info, "summary", Dict{String,Any}())
    g(k) = Float64(get(s, k, NaN))
    added, removed = g("vol_added_m3"), g("vol_removed_m3")
    err = g("conservation_error_m3")
    return (t = t, stored = S, t_mid = tm, rate = rate, added = added, removed = removed,
            net = added - removed, stored_end = g("vol_domain_m3"),
            error_rel = added > 0 ? err / added : NaN)
end

"""
    volume_scale(vmax) -> (factor, unit)

Unit for volumes up to `vmax` m³: (1, "m³"), (1e3, "10³ m³") or (1e6, "10⁶ m³").
"""
volume_scale(v::Real) = abs(v) >= 1e6 ? (1e6, "10⁶ m³") : abs(v) >= 1e4 ? (1e3, "10³ m³") : (1.0, "m³")

"""
    format_volume(v) -> String

A volume for labels, three significant figures in `volume_scale`'s unit:
`format_volume(5.92e6) == "5.92 × 10⁶ m³"`.
"""
function format_volume(v::Real)
    isfinite(v) || return "–"
    f, u = volume_scale(v)
    s = @sprintf("%.3g", v / f)
    return f == 1 ? s * " " * u : s * " × " * u
end
