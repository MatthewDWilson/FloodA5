# visualisation/plots/makie.jl — the CairoMakie part of FloodA5Plots.
# Included by FloodA5Plots.jl only when CairoMakie loads; the Makie-free logic
# it uses is in core.jl.

# ───────────────────────────── theme ───────────────────────────────────────

"""
    house_theme(; fontsize = 13) -> Theme

The FloodA5 figure style: surface background, recessive grid and axes, left
bold titles, ink-coloured text. `set_theme!(FP.house_theme())` applies it to
every later figure; the FP plotting functions use it themselves.
"""
house_theme(; fontsize::Real = 13) = Theme(
    fontsize = fontsize, textcolor = INK, backgroundcolor = SURF,
    Axis = (backgroundcolor = SURF, xgridcolor = GRID, ygridcolor = GRID,
            leftspinecolor = MUTED, bottomspinecolor = MUTED, topspinevisible = false,
            rightspinevisible = false, xtickcolor = MUTED, ytickcolor = MUTED,
            xticklabelcolor = INK2, yticklabelcolor = INK2, titlealign = :left,
            titlefont = :bold, subtitlecolor = INK2),
    Colorbar = (ticklabelcolor = INK2, tickcolor = MUTED, spinewidth = 0.0,
                labelcolor = INK),
    Legend = (framevisible = false, labelcolor = INK2))

"Activate `house_theme()` for all later figures (notebooks)."
activate!(; kw...) = set_theme!(house_theme(; kw...))

_themed(f) = with_theme(f, house_theme())

"Map axes: no grid, all four spines (a map has no baseline)."
const _MAP_AXIS = (xgridvisible = false, ygridvisible = false, topspinevisible = true,
                   rightspinevisible = true, topspinecolor = MUTED, rightspinecolor = MUTED,
                   xticklabelsize = 11, yticklabelsize = 11, xticks = WilkinsonTicks(5; k_min = 3),
                   yticks = WilkinsonTicks(5; k_min = 3))

# ───────────────────────────── drawing cells ───────────────────────────────

const _NCOL = 256

"""
    _cell_mesh(polys) -> (points, faces, owner)

Fan triangulation of convex polygons (2 × nv closed rings): one vertex set
per cell, so a per-vertex colour equal to the cell's colour draws flat cells.
`owner[v]` is the cell of vertex v. One `mesh!` is faster than one path per
cell in CairoMakie (Carlisle, 30k cells: 0.8 vs 1.2 s per map).
"""
function _cell_mesh(polys::AbstractVector)
    nv = sum(P -> size(P, 2) - (P[:, 1] == P[:, end]), polys; init = 0)
    pts = Vector{Point2f}(undef, nv); owner = Vector{Int}(undef, nv)
    faces = Makie.GLTriangleFace[]; sizehint!(faces, nv)
    b = 0
    for (i, P) in enumerate(polys)
        n = size(P, 2) - (P[:, 1] == P[:, end])
        for j in 1:n
            pts[b + j] = Point2f(P[1, j], P[2, j]); owner[b + j] = i
        end
        for j in 2:n-1
            push!(faces, Makie.GLTriangleFace(b + 1, b + j, b + j + 1))
        end
        b += n
    end
    return pts, faces, owner
end

"Draw `polys` coloured by `values` on `ramp` over `range` (NaN: not drawn)."
function _draw_cells!(ax, polys, values, ramp, range)
    keep = [i for i in eachindex(values) if isfinite(values[i])]
    isempty(keep) && return nothing
    cols = resample_cmap(ramp, _NCOL)
    pts, faces, owner = _cell_mesh(polys[keep])
    c = [cols[color_index(values[keep[o]], range, _NCOL)] for o in owner]
    return mesh!(ax, pts, faces; color = c, shading = NoShading)
end

"Draw `polys` in one colour."
function _draw_flat!(ax, polys, color)
    isempty(polys) && return nothing
    pts, faces, _ = _cell_mesh(polys)
    return mesh!(ax, pts, faces; color = Makie.to_color(color), shading = NoShading)
end

function _basemap!(ax, polys, run::FA.Run, basemap)
    basemap === :none && return nothing
    basemap === :grey && return _draw_flat!(ax, polys, BASEMAP_GREY)
    basemap === :elevation || error("basemap must be :elevation, :grey or :none")
    z = run.mesh.elevation
    zf = filter(isfinite, z)
    # a flat bed (planar tests without a DEM) is drawn in the flat grey
    (isempty(zf) || maximum(zf) - minimum(zf) < 1e-6) && return _draw_flat!(ax, polys, BASEMAP_GREY)
    return _draw_cells!(ax, polys, z, TERRAIN, robust_range(z; nice = false))
end

"""
    MapPanel

What a map-drawing function returns: the `axis`, its `coords::MapCoords`
(pass it to `mark_points!` and to further maps that should share the
origin), the colour `range` and `colormap` drawn.
"""
struct MapPanel
    axis::Axis
    coords::MapCoords
    range::Tuple{Float64,Float64}
    colormap::Any
end

# `width` (px): a fixed-size axis drawn to scale (then `resize_to_layout!`);
# `nothing`: the axis fills its grid cell at the true aspect ratio.
function _map_axis(gp, mc::MapCoords, polys; title = "", subtitle = "", width = nothing, kw...)
    xl, yl, asp = map_extent(mc, polys)
    xlab, ylab = axis_labels(mc)
    size_kw = width === nothing ? (aspect = AxisAspect(asp),) :
              (width = Float64(width), height = Float64(width) / asp)
    return Axis(gp; _MAP_AXIS..., size_kw..., limits = (xl, yl),
                xlabel = xlab, ylabel = ylab, title, subtitle, kw...)
end

# ───────────────────────────── points ──────────────────────────────────────

"""
    mark_points!(panel_or_axis, [coords,] points; labels = nothing,
                 color = INK, marker = :circle, markersize = 9, fontsize = 11)

Markers (with optional labels) at `points` — `(lon, lat)` tuples or named
tuples with `lon`, `lat` (and `label`, used if `labels` is not given) — on a
map drawn with `coords`. A `MapPanel` carries its own coords.

    p = FP.fieldmap!(fig[1, 1], r, 10)
    FP.mark_points!(p, [(-2.93, 54.89)]; labels = ["gauge"])
"""
function mark_points!(ax::Axis, mc::MapCoords, points; labels = nothing, color = INK,
                      marker = :circle, markersize = 9, fontsize = 11)
    isempty(points) && return nothing
    lon = [p isa NamedTuple ? p.lon : p[1] for p in points]
    lat = [p isa NamedTuple ? p.lat : p[2] for p in points]
    x, y = project(mc, lon, lat)
    scatter!(ax, x, y; color = color, marker, markersize, strokecolor = SURF, strokewidth = 1.5)
    labs = labels !== nothing ? labels :
           [p isa NamedTuple && hasproperty(p, :label) ? p.label : "" for p in points]
    for (xi, yi, l) in zip(x, y, labs)
        isempty(l) && continue
        text!(ax, xi, yi; text = String(l), color = INK, fontsize, align = (:left, :bottom),
              offset = (5, 4))
    end
    return nothing
end
mark_points!(p::MapPanel, points; kw...) = mark_points!(p.axis, p.coords, points; kw...)

"""
    mark_sources!(panel, run; labels = true, kw...)

Mark where `run` adds water (`run_sources`: inflow, injection and rain
points, `.bci` points) on a map.
"""
function mark_sources!(p::MapPanel, r::FA.Run; labels::Bool = true, kw...)
    src = run_sources(r)
    # only points inside the map (a .bci may list points outside the mesh)
    (x0, x1), (y0, y1) = p.axis.limits[]
    src = filter(s -> (xy = project(p.coords, s.lon, s.lat);
                       x0 <= xy[1] <= x1 && y0 <= xy[2] <= y1), src)
    mark_points!(p, src; labels = labels ? nothing : fill("", length(src)),
                 marker = :star5, markersize = 13, kw...)
end

# ───────────────────────────── field maps ──────────────────────────────────

function _resolve_k(r::FA.Run, k, t)
    t !== nothing && return FA.frame_at(r, t; tol = frame_tol(r))
    k === nothing && return FA.nframes(r)
    FA.check_frame(r, k)
    return Int(k)
end

"Values to draw for `var` at frame k, NaN on dry cells (if the variable is masked)."
function _masked_values(r::FA.Run, var, k::Integer, dry_threshold::Real)
    v = frame_values(r, var, k)
    if var_spec(var).masked
        h = Symbol(var) === :depth ? v : FA.depth(r, k)
        v = [h[i] > dry_threshold ? v[i] : NaN for i in eachindex(v)]
    end
    return v
end

_default_cmap(var) = SEQUENTIAL

"""
    fieldmap!(gp, run, k = last; t = nothing, var = :depth, colorrange = nothing,
              colormap = nothing, dry_threshold = 0.01, basemap = :elevation,
              coords = :lonlat, origin = :centre, sources = true,
              colorbar = true, title = nothing, subtitle = nothing) -> MapPanel

Draw one frame of a run as a map of its cells (the mesh polygons) into the
grid position `gp` (an `Axis` goes in `gp`, its colourbar beside it unless
`colorbar = false`).

- frame: `k`, or the frame nearest time `t` (s).
- `var`: `:depth`, `:wse`, `:speed`, `:volume`, `:elevation` or an output
  variable name (`var_spec`). Depth-like variables are not drawn where the
  depth is ≤ `dry_threshold` (m); those cells show the `basemap`:
  `:elevation` (grey shades), `:grey` (flat) or `:none` (transparent).
- `colorrange`: default `default_range` — 0 to the rounded 98th percentile
  of the wet values; values outside are clamped (the colourbar shows a
  triangle).
- `coords = :lonlat` (degrees, true aspect) or `:local` (metres about
  `origin`: `:centre`, `:source` or `(lon, lat)`); or a `MapCoords` to share
  with other panels.
- `sources`: mark the run's water sources (`mark_sources!`).

    fig = Figure(); FP.fieldmap!(fig[1, 1], r, nothing; t = 7200); fig
"""
function fieldmap!(gp, r::FA.Run, k = nothing; t = nothing, var = :depth, colorrange = nothing,
                   colormap = nothing, dry_threshold::Real = 0.01, basemap = :elevation,
                   coords = :lonlat, origin = :centre, sources::Bool = true,
                   colorbar::Bool = true, title = nothing, subtitle = nothing,
                   colorbar_label = nothing, width = nothing, axis = (;))
    k = _resolve_k(r, k, t)
    mc = coords isa MapCoords ? coords : map_coords(r; coords, origin)
    polys = project_polygons(mc, r.mesh)
    v = _masked_values(r, var, k, dry_threshold)
    rng = colorrange === nothing ? default_range(var, v) : (Float64(colorrange[1]), Float64(colorrange[2]))
    cm = colormap === nothing ? _default_cmap(var) : colormap
    axgp = colorbar ? gp[1, 1] : gp
    ax = _map_axis(axgp, mc, polys;
                   title = title === nothing ? run_label(r) : title,
                   subtitle = subtitle === nothing ? time_label(r.t[k]) * " · " * lowercase(var_spec(var).label) : subtitle,
                   width, axis...)
    _basemap!(ax, polys, r, basemap)
    _draw_cells!(ax, polys, v, cm, rng)
    p = MapPanel(ax, mc, rng, cm)
    sources && mark_sources!(p, r)
    colorbar && _colorbar!(gp[1, 2], cm, rng, v; label = something(colorbar_label, var_label(var)))
    return p
end

"Colourbar for `ramp` over `range`, with clip triangles if `values` exceed it."
function _colorbar!(gp, ramp, range, values; label = "", vertical = true, kw...)
    fin = filter(isfinite, values)
    hi = !isempty(fin) && maximum(fin) > range[2] ? resample_cmap(ramp, 2)[end] : nothing
    lo = !isempty(fin) && minimum(fin) < range[1] ? resample_cmap(ramp, 2)[1] : nothing
    return Colorbar(gp; colormap = cgrad(ramp), limits = range, label, vertical,
                    highclip = hi, lowclip = lo, width = vertical ? 12 : Makie.automatic,
                    height = vertical ? Relative(0.8) : 12, ticks = WilkinsonTicks(5), kw...)
end

"""
    fieldmap(run, k = last; width = 560, kw...) -> Figure

One map, `width` pixels wide, the figure sized to fit (keywords: `fieldmap!`).

    FP.fieldmap(r; t = 4 * 3600)                       # depth at 4 h
    FP.fieldmap(r, 12; var = :speed, coords = :local)
"""
function fieldmap(r::FA.Run, k = nothing; width = 560, kw...)
    _themed() do
        fig = Figure()
        fieldmap!(fig[1, 1], r, k; width, kw...)
        resize_to_layout!(fig)
        fig
    end
end

"""
    fieldmaps(runs, times; var = :depth, colorrange = nothing, coords = :lonlat,
              origin = :centre, tol = nothing, panel = 300, labels = nothing,
              kw...) -> Figure

A grid of maps: one row per run (a single `Run` is one row), one column per
time (s; each run's nearest frame within `tol`, `frame_grid`), with one
colour range (over all panels unless given) and one colourbar. Rows are
labelled with the run labels (or `labels`), columns with the times. `panel`
is the width of one map in pixels. `colormap`, `dry_threshold`, `basemap`,
`sources`: as `fieldmap!`.

    FP.fieldmaps([li, swe], [2, 4, 6] .* 3600)
    FP.fieldmaps([a5, cart], [600, 1800, 3600]; coords = :local, origin = :source)
"""
function fieldmaps(runs::AbstractVector, times::AbstractVector; var = :depth, colorrange = nothing,
                   colormap = nothing, coords = :lonlat, origin = :centre, tol = nothing,
                   panel::Real = 300, labels = nothing, dry_threshold::Real = 0.01,
                   sources::Bool = true, basemap = :elevation)
    K = frame_grid(runs, times; tol)
    mc = coords isa MapCoords ? coords : map_coords(runs[1]; coords, origin)
    vals = [_masked_values(r, var, K[i, j], dry_threshold) for (i, r) in enumerate(runs), j in eachindex(times)]
    rng = colorrange !== nothing ? (Float64(colorrange[1]), Float64(colorrange[2])) :
          default_range(var, reduce(vcat, vec(vals)))
    cm = colormap === nothing ? _default_cmap(var) : colormap
    polys1 = project_polygons(mc, runs[1].mesh)
    nr, nc = length(runs), length(times)
    labs = labels === nothing ? run_label.(runs) : labels
    _themed() do
        fig = Figure()
        for (i, r) in enumerate(runs)
            polys = i == 1 ? polys1 : project_polygons(mc, r.mesh)
            for j in 1:nc
                k = K[i, j]
                ax = _map_axis(fig[i, j], mc, polys; title = labs[i],
                               subtitle = time_label(r.t[k]), width = panel)
                _basemap!(ax, polys, r, basemap)
                _draw_cells!(ax, polys, vals[i, j], cm, rng)
                sources && mark_sources!(MapPanel(ax, mc, rng, cm), r; labels = false)
                i < nr && hidexdecorations!(ax; grid = false, ticks = false)
                j > 1 && hideydecorations!(ax; grid = false, ticks = false)
            end
        end
        _colorbar!(fig[1:nr, nc + 1], cm, rng, reduce(vcat, vec(vals)); label = var_label(var))
        colgap!(fig.layout, 12); rowgap!(fig.layout, 10)
        resize_to_layout!(fig)
        fig
    end
end
fieldmaps(r::FA.Run, times::AbstractVector; kw...) = fieldmaps([r], times; kw...)

# ───────────────────────────── comparisons ─────────────────────────────────

"""
    difference_map!(gp, ref, run; t = nothing, k = nothing, var = :depth,
                    wet = 0.01, colorrange = nothing, coords = :lonlat,
                    origin = :centre, basemap = :elevation, colorbar = true,
                    title = nothing) -> MapPanel

`var` of `run` minus `var` of `ref` at one time (`matched_pair`: `t`, or
frame `k` of `run`, or the last common time), on the cells wet (depth >
`wet`) in either; other cells show the basemap. Diverging colours: blue =
run higher/deeper, red = lower/shallower; range symmetric about 0 at the
rounded 98th percentile of |Δ|. The runs must share a mesh.

    FP.difference_map(swe, li; t = 6 * 3600)
"""
function difference_map!(gp, ref::FA.Run, run::FA.Run; t = nothing, k = nothing, var = :depth,
                         wet::Real = 0.01, colorrange = nothing, coords = :lonlat,
                         origin = :centre, basemap = :elevation, colorbar::Bool = true,
                         title = nothing, sources::Bool = false, width = nothing)
    same_mesh(ref, run)
    kr, kk = matched_pair(ref, run; t, k)
    h_ref, h = FA.depth(ref, kr), FA.depth(run, kk)
    v_ref = Symbol(var) === :depth ? h_ref : frame_values(ref, var, kr)
    v = Symbol(var) === :depth ? h : frame_values(run, var, kk)
    d = difference_values(v_ref, v, h_ref, h; wet)
    rng = colorrange === nothing ? robust_range(d; symmetric = true) :
          (Float64(colorrange[1]), Float64(colorrange[2]))
    w = filter(isfinite, d)
    rms = isempty(w) ? NaN : sqrt(sum(abs2, w) / length(w))
    mc = coords isa MapCoords ? coords : map_coords(run; coords, origin)
    polys = project_polygons(mc, run.mesh)
    sym = Symbol(var) === :depth ? "Δh" : "Δ"
    u = var_spec(var).unit
    ax = _map_axis(colorbar ? gp[1, 1] : gp, mc, polys;
                   title = something(title, run_label(run) * " − " * run_label(ref)),
                   subtitle = time_label(run.t[kk]) * " · " * sym * @sprintf(" rms %.3g ", rms) * u *
                              " over cells wet in either", width)
    _basemap!(ax, polys, run, basemap)
    _draw_cells!(ax, polys, d, DIVERGING, rng)
    p = MapPanel(ax, mc, rng, DIVERGING)
    sources && mark_sources!(p, run; labels = false)
    colorbar && _colorbar!(gp[1, 2], DIVERGING, rng, d;
                           label = sym * " " * lowercase(var_spec(var).label) * (isempty(u) ? "" : " ($u)"))
    return p
end

"""
    difference_map(ref, run; width = 560, kw...) -> Figure

See `difference_map!`; the map is `width` pixels wide.
"""
function difference_map(ref::FA.Run, run::FA.Run; width = 560, kw...)
    _themed() do
        fig = Figure()
        difference_map!(fig[1, 1], ref, run; width, kw...)
        resize_to_layout!(fig)
        fig
    end
end

"""
    extent_map!(gp, ref, run; t = nothing, k = nothing, threshold = 0.1,
                coords = :lonlat, origin = :centre, basemap = :elevation,
                legend = true, title = nothing) -> MapPanel

Flood extent of `run` against `ref` at one time (`matched_pair`): cells wet
(depth > `threshold`, m) in both, in the reference only, in the run only
(`extent_classes`); the subtitle gives the CSI and counts (as
`FA.compare`). The runs must share a mesh.

    FP.extent_map(swe, li; t = 6 * 3600, threshold = 0.1)
"""
function extent_map!(gp, ref::FA.Run, run::FA.Run; t = nothing, k = nothing, threshold::Real = 0.1,
                     coords = :lonlat, origin = :centre, basemap = :elevation,
                     legend::Bool = true, title = nothing, width = nothing)
    same_mesh(ref, run)
    kr, kk = matched_pair(ref, run; t, k)
    c = extent_classes(FA.depth(ref, kr), FA.depth(run, kk); threshold)
    s = extent_scores(c)
    mc = coords isa MapCoords ? coords : map_coords(run; coords, origin)
    polys = project_polygons(mc, run.mesh)
    ax = _map_axis(legend ? gp[1, 1] : gp, mc, polys;
                   title = something(title, @sprintf("%s vs %s: CSI %.2f", run_label(run), run_label(ref), s.csi)),
                   subtitle = time_label(run.t[kk]) * @sprintf(" · wet = depth > %g m", threshold), width)
    _basemap!(ax, polys, run, basemap)
    for cls in 1:3
        sel = findall(==(cls), c)
        _draw_flat!(ax, polys[sel], EXTENT_COLORS[cls])
    end
    if legend
        n = (s.hits, s.misses, s.false_alarms)
        labs = [replace(EXTENT_LABELS[i], "Reference" => run_label(ref), "Run" => run_label(run)) *
                "  " * string(n[i]) for i in 1:3]
        Legend(gp[1, 2], [PolyElement(color = EXTENT_COLORS[i], strokecolor = :transparent) for i in 1:3],
               labs, "Cells"; framevisible = false, labelcolor = INK2, titlefont = :bold,
               titlehalign = :left, gridshalign = :left)
    end
    return MapPanel(ax, mc, (0.0, 3.0), EXTENT_COLORS)
end

"""
    extent_map(ref, run; width = 560, kw...) -> Figure

See `extent_map!`; the map is `width` pixels wide.
"""
function extent_map(ref::FA.Run, run::FA.Run; width = 560, kw...)
    _themed() do
        fig = Figure()
        extent_map!(fig[1, 1], ref, run; width, kw...)
        resize_to_layout!(fig)
        fig
    end
end

# ───────────────────────────── time series ─────────────────────────────────

_tsc(ts, t) = t ./ ts

"Legend entries for line styles from `run_styles`."
_line_elements(styles) = [[LineElement(color = s.color, linestyle = s.linestyle, linewidth = 2)]
                          for s in styles]

"Output variable read by `timeseries` for a map variable name."
_series_var(var) = Symbol(var) === :depth ? "water_depth" : Symbol(var) === :speed ? "velocity" :
                   String(Symbol(var))

"""
    timeseries!(gp, runs, points; var = :depth, labels = nothing,
                point_labels = nothing, style = :auto) -> Vector{Axis}

`var` (an output variable: `:depth`, `:speed` = the `velocity` variable,
`"volume"`, …) against time at `points` (cell indices, `(lon, lat)` or
`(lon = …, lat = …, label = …)`; the nearest cell, `point_cells`). One run:
one axis, a line per point. Several runs (each run's cell nearest each
point): one axis per point, a line per run (`run_styles(runs; style)`),
with a shared legend.

    FP.timeseries([li, swe], [(lon = -2.93, lat = 54.895, label = "floodplain")])
"""
function timeseries!(gp, runs::AbstractVector, points; var = :depth, labels = nothing,
                     point_labels = nothing, style::Symbol = :auto)
    v = _series_var(var)
    pl = point_labels !== nothing ? point_labels :
         [p isa NamedTuple && hasproperty(p, :label) ? p.label :
          p isa Integer ? "cell $p" : @sprintf("(%.4f, %.4f)", p[1], p[2]) for p in points]
    rl = labels === nothing ? run_label.(runs) : labels
    tmax = maximum(r.t[end] for r in runs)
    ts, tlab = time_axis(tmax)
    axes_ = Axis[]
    if length(runs) == 1
        r = runs[1]; cells = point_cells(r.mesh, points)
        ax = Axis(gp[1, 1]; xlabel = tlab, ylabel = var_label(var),
                  title = run_label(r))
        S = FA.series(r, v, cells)
        st = [(color = series_color(j), linestyle = LINESTYLES[j]) for j in eachindex(cells)]
        for j in eachindex(cells)
            lines!(ax, _tsc(ts, r.t), S[:, j]; color = st[j].color, linestyle = st[j].linestyle, linewidth = 2)
        end
        length(cells) > 1 && Legend(gp[2, 1], _line_elements(st), pl; orientation = :horizontal,
                                    framevisible = false, labelcolor = INK2)
        push!(axes_, ax)
        return axes_
    end
    st = run_styles(runs; style)
    for (j, p) in enumerate(points)
        ax = Axis(gp[1, j]; xlabel = tlab, ylabel = j == 1 ? var_label(var) : "",
                  title = pl[j])
        for (i, r) in enumerate(runs)
            c = point_cells(r.mesh, [p])
            lines!(ax, _tsc(ts, r.t), vec(FA.series(r, v, c)); color = st[i].color,
                   linestyle = st[i].linestyle, linewidth = 2)
        end
        push!(axes_, ax)
    end
    Legend(gp[2, 1:length(points)], _line_elements(st), rl;
           orientation = :horizontal, framevisible = false, labelcolor = INK2)
    return axes_
end
timeseries!(gp, r::FA.Run, points; kw...) = timeseries!(gp, [r], points; kw...)

"""
    timeseries(runs, points; size = nothing, kw...) -> Figure

See `timeseries!`.
"""
function timeseries(runs, points; size = nothing, kw...)
    rs = runs isa FA.Run ? [runs] : runs
    _themed() do
        n = length(rs) == 1 ? 1 : length(points)
        fig = Figure(size = something(size, (max(560, 380 * n), 400)))
        timeseries!(fig[1, 1], rs, points; kw...)
        fig
    end
end

"""
    storage_plot!(gp, runs; labels = nothing, style = :auto,
                  title = "Water stored in the domain") -> Axis

Storage (`FA.storage`: Σ depth·area, Σ volume for SGS) against time for
several runs on one axis, styled by `run_styles(runs; style)` (`:auto`:
A5 blue / Cartesian orange, local-inertial solid / SWE dashed / SGS dotted
when that identifies every run, else series colours), with a legend and —
for ≤ 4 runs — labels at the line ends.

    FP.storage_plot([a5_li, a5_sgs, a5_swe, cart_li, cart_swe])
"""
function storage_plot!(gp, runs::AbstractVector; labels = nothing, style::Symbol = :auto,
                       title = "Water stored in the domain", series = nothing)
    S = series === nothing ? [storage_series(r) for r in runs] : series
    vmax = maximum(maximum(filter(isfinite, s); init = 0.0) for s in S)
    f, unit = volume_scale(vmax)
    tmax = maximum(r.t[end] for r in runs)
    ts, tlab = time_axis(tmax)
    rl = labels === nothing ? run_label.(runs) : labels
    st = run_styles(runs; style)
    ax = Axis(gp[1, 1]; xlabel = tlab, ylabel = "Storage ($unit)", title)
    for (i, r) in enumerate(runs)
        lines!(ax, _tsc(ts, r.t), S[i] ./ f; color = st[i].color, linestyle = st[i].linestyle,
               linewidth = 2)
    end
    ylims!(ax, 0, nothing)
    if length(runs) <= 4        # direct labels at the line ends, nudged apart
        ends = [S[i][end] / f for i in eachindex(runs)]
        ys = _spread_labels(ends, 0.045 * vmax / f)
        for (i, r) in enumerate(runs)
            text!(ax, _tsc(ts, r.t[end]), ys[i]; text = rl[i], color = INK2, fontsize = 11,
                  align = (:left, :center), offset = (6, 0))
        end
        xlims!(ax, 0, _tsc(ts, tmax) * 1.18)
    end
    Legend(gp[2, 1], _line_elements(st), rl; orientation = :horizontal,
           framevisible = false, labelcolor = INK2, nbanks = length(runs) > 4 ? 2 : 1)
    return ax
end

"Positions for labels at `y` at least `gap` apart (order preserved)."
function _spread_labels(y::AbstractVector, gap::Real)
    o = sortperm(y); z = float.(y[o])
    for i in 2:length(z)
        z[i] = max(z[i], z[i-1] + gap)
    end
    out = similar(z); out[o] = z
    return out
end

"""
    storage_plot(runs; size = nothing, kw...) -> Figure

See `storage_plot!`.
"""
function storage_plot(runs::AbstractVector; size = nothing, kw...)
    _themed() do
        fig = Figure(size = something(size, (760, 460)))
        storage_plot!(fig[1, 1], runs; kw...)
        fig
    end
end

"""
    mass_balance_plot!(gp, run) -> (Axis, Axis)

Two panels from the per-frame `volume` and the run summary
(`mass_balance_series`): the stored volume Σ V against time, with the
summary's added − removed at the end (they agree to the conservation error,
quoted in the subtitle); and the net inflow rate dS/dt between frames.

    FP.mass_balance_plot(r)
"""
function mass_balance_plot!(gp, r::FA.Run)
    gl = gp isa GridLayout ? gp : GridLayout(gp)
    mb = mass_balance_series(r)
    smax = maximum(filter(isfinite, vcat(mb.stored, mb.net)); init = 0.0)
    f, unit = volume_scale(smax)
    ts, tlab = time_axis(r.t[end])
    ax1 = Axis(gl[1, 1]; ylabel = "Stored ($unit)",
               title = run_label(r) * ": mass balance",
               subtitle = isfinite(mb.added) ?
                   "added " * format_volume(mb.added) * " · removed " * format_volume(mb.removed) *
                   @sprintf(" · conservation error %.1e of added", mb.error_rel) :
                   "no run summary (unfinished run?)")
    lines!(ax1, _tsc(ts, mb.t), mb.stored ./ f; color = SERIES[1], linewidth = 2)
    if isfinite(mb.net)
        scatter!(ax1, [_tsc(ts, mb.t[end])], [mb.net / f]; color = SURF, strokecolor = INK,
                 strokewidth = 1.5, markersize = 10)
        text!(ax1, _tsc(ts, mb.t[end]), mb.net / f; text = "added − removed (summary)", color = INK2,
              fontsize = 11, align = (:right, :bottom), offset = (-8, 8))
    end
    ylims!(ax1, 0, 1.15 * smax / f)
    hidexdecorations!(ax1; grid = false, ticks = false)
    ax2 = Axis(gl[2, 1]; xlabel = tlab, ylabel = "dS/dt (m³/s)")
    hlines!(ax2, [0.0]; color = MUTED, linestyle = :dash, linewidth = 1)
    isempty(mb.rate) || stairs!(ax2, _tsc(ts, mb.t), vcat(mb.rate, mb.rate[end]); step = :post,
                                color = SERIES[1], linewidth = 2)
    linkxaxes!(ax1, ax2)
    rowsize!(gl, 2, Relative(0.4))
    return ax1, ax2
end

"""
    mass_balance_plot(run; size = nothing) -> Figure

See `mass_balance_plot!`.
"""
function mass_balance_plot(r::FA.Run; size = nothing)
    _themed() do
        fig = Figure(size = something(size, (720, 520)))
        mass_balance_plot!(fig[1, 1], r)
        fig
    end
end
