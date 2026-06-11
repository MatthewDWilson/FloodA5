#!/usr/bin/env julia
# =============================================================================
#  visualise_mesh.jl
#  -------------------------------------------------------------------
#  Produces three publication-quality figures for a given A5 cell:
#
#  Fig 1 — Wire-frame mesh diagram
#          Shows the target cell and its immediate neighbours with
#          directional flow arrows from a representative FlowState
#          (gravity-driven from each cell to its lowest neighbour).
#
#  Fig 2 — SGS edge-sill cross-section
#          Annotated terrain profile along the shared edge between the
#          target cell and a chosen neighbour slot, showing z_sill, the
#          water level at a representative stage, and A / P geometry.
#
#  Fig 3 — SGS hypsometric curve
#          Elevation vs. cumulative volume for the target cell, with
#          z_min, z_max, z_sill annotated.
#
#  Usage:
#    julia visualise_mesh.jl --mesh PATH.parquet --cell CELL_ID_HEX
#                            [--outdir ./figures] [--slot 1]
#                            [--wse WSE_M]
#
#  Requirements: Pkg.add(["CairoMakie","Arrow","DataFrames","ArgParse"])
#  The script does NOT call A5Grid.jl / FloodModel.jl — it reads the
#  parquet directly so it can run standalone without the full model env.
# =============================================================================

using DataFrames, CairoMakie, ArgParse, Statistics, Printf, LinearAlgebra

# Include A5Grid for load_mesh_geoparquet — adjust path if script is not in
# a subdirectory of the FloodA5 root (e.g. use "../A5Grid.jl" from visualisation/)
const _A5GRID_PATH = joinpath(@__DIR__, "..", "mesh", "A5Grid.jl")
include(_A5GRID_PATH)
using .A5Grid

# ── Package bootstrap (run once if packages are missing) ─────────────────
# import Pkg; Pkg.add(["DataFrames", "CairoMakie", "ArgParse"])


# ─────────────────────────────────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────────────────────────────────
function parse_args_()
    s = ArgParseSettings(prog="visualise_mesh.jl",
        description="Generate mesh-structure and SGS diagnostic figures.")
    @add_arg_table! s begin
        "--mesh"
            arg_type = String
            required = true
            help     = "Path to .parquet mesh file"
        "cell"
            arg_type = String
            required = true
            help     = "Target cell ID (16-char hex), e.g. 6345f2518e800000"
        "--outdir"
            arg_type = String
            default  = "."
            help     = "Output directory for figures (created if absent)"
        "--slot"
            arg_type = Int
            default  = 1
            help     = "Neighbour adjacency slot for edge-sill plot (1-5)"
        "--wse"
            arg_type = Float64
            default  = NaN
            help     = "Water surface elevation for SGS annotation (m); default: z_min + 0.5*(z_max-z_min)"
        "--dpi"
            arg_type = Int
            default  = 150
            help     = "PNG raster DPI"
    end
    return parse_args(s)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Geometry helpers
# ─────────────────────────────────────────────────────────────────────────────
const R_EARTH = 6_371_000.0   # m

function haversine_m(lon1, lat1, lon2, lat2)
    φ1, φ2 = deg2rad(lat1), deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)
    a  = sin(Δφ/2)^2 + cos(φ1)*cos(φ2)*sin(Δλ/2)^2
    return R_EARTH * 2 * atan(sqrt(a), sqrt(1-a))
end

"""Project boundary vertices to local East-North metres centred on (cx, cy)."""
function to_local_m(boundary::Vector, cx::Float64, cy::Float64)
    cos_cy = cosd(cy)
    xs = [(v[1] - cx) * deg2rad(1) * R_EARTH * cos_cy for v in boundary]
    ys = [(v[2] - cy) * deg2rad(1) * R_EARTH           for v in boundary]
    return xs, ys
end

function polygon_area_m2(bnd)
    lons = [v[1] for v in bnd]
    lats = [v[2] for v in bnd]
    cx, cy = mean(lons), mean(lats)
    xs, ys = to_local_m(bnd, cx, cy)
    n = length(xs)
    area = sum(xs[i]*ys[mod1(i+1,n)] - xs[mod1(i+1,n)]*ys[i] for i in 1:n)
    return abs(area) / 2.0
end


# ─────────────────────────────────────────────────────────────────────────────
#  Parquet loading

"""
Load the parquet mesh via A5Grid.load_mesh_geoparquet (the same path the main
model uses) and return the DataFrame-like accessors the plotting functions need:
  df   — a NamedTuple wrapping the A5Mesh for column-style access
  idx  — Dict{cell_id_hex => cell_index}
  adj  — Dict{cell_id_hex => Vector{neighbour_id_hex}}
  has_sgs — Bool
"""
function load_parquet(path::String)
    mesh = A5Grid.load_mesh_geoparquet(path)

    norm(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
    n     = length(mesh.cells)
    ids   = [norm(c.id) for c in mesh.cells]
    idx   = Dict{String,Int}(ids[i] => i for i in 1:n)

    # Adjacency from mesh.adjacency (already normalised by load_mesh_geoparquet)
    adj = Dict{String,Vector{String}}()
    for (i, cid) in enumerate(ids)
        raw_nbs = get(mesh.adjacency, cid, String[])
        adj[cid] = [norm(nb) for nb in raw_nbs]
    end

    has_sgs = haskey(mesh.array_vars, "sgs_elev_bins")

    # Build a lightweight DataFrame with only the columns the figures need.
    # SGS array columns are stored as matrices (n_bins × n_cells) in mesh.array_vars;
    # convert each to a Vector-of-Vectors for the row-oriented access pattern.
    df = DataFrame(
        cell_id      = ids,
        cell_id_norm = ids,
        center_lon   = [c.center_lon for c in mesh.cells],
        center_lat   = [c.center_lat for c in mesh.cells],
        boundary     = [c.boundary   for c in mesh.cells],
    )

    # Scalar static vars
    for (col, key) in [(:elevation, "elevation"), (:sgs_cell_area, "sgs_cell_area"),
                       (:sgs_z_min, "sgs_z_min"), (:sgs_z_max, "sgs_z_max")]
        haskey(mesh.static_vars, key) && (df[!, col] = mesh.static_vars[key])
    end

    # Adjacency columns
    for (slot, label) in enumerate(["adj_0","adj_1","adj_2","adj_3","adj_4"])
        df[!, label] = [let nbs = get(mesh.adjacency, ids[i], String[]);
                            slot <= length(nbs) ? nbs[slot] : "" end
                        for i in 1:n]
    end

    # SGS array vars: matrix columns → Vector{Vector{Float64}}
    for (col, key) in [(:sgs_elev_bins, "sgs_elev_bins"),
                       (:sgs_vol_curve, "sgs_vol_curve"),
                       (:sgs_area_curve, "sgs_area_curve"),
                       (:sgs_edge_sills, "sgs_edge_sills")]
        haskey(mesh.array_vars, key) || continue
        mat = mesh.array_vars[key]
        df[!, col] = [mat[:, i] for i in 1:n]
    end

    # Edge R-A curves: stored as (n_bins*5 × n_cells); expose as flat Vector per row
    for (col, key) in [(:sgs_edge_area_curve, "sgs_edge_area_curve"),
                       (:sgs_edge_perim_curve, "sgs_edge_perim_curve")]
        haskey(mesh.array_vars, key) || continue
        mat = mesh.array_vars[key]
        df[!, col] = [mat[:, i] for i in 1:n]
    end

    return df, idx, adj, has_sgs
end


"""
Extract boundary as Vector{Tuple{Float64,Float64}} from WKB geometry column.
Falls back to a rough pentagon constructed from center + approximate radius
if the geometry cannot be decoded.
"""
function get_boundary(df, row::Int)
    # Use real pentagon boundary stored in :boundary column (set by load_parquet
    # from mesh.cells[row].boundary — a Vector of [lon, lat] pairs)
    if :boundary in propertynames(df)
        bnd = df[row, :boundary]   # Vector{Vector{Float64}} from A5Grid
        return [(v[1], v[2]) for v in bnd]
    end

    # Fallback: approximate regular pentagon from cell centre + area
    lon = df[row, :center_lon]
    lat = df[row, :center_lat]
    r_m = :sgs_cell_area in propertynames(df) ?
          sqrt(df[row, :sgs_cell_area] / π) : 200.0
    r_lon = rad2deg(r_m / (R_EARTH * cosd(lat)))
    r_lat = rad2deg(r_m / R_EARTH)
    pts = [(lon + r_lon * cosd(90 + 72k), lat + r_lat * sind(90 + 72k))
           for k in 0:4]
    push!(pts, pts[1])  # close
    return pts
end

# ─────────────────────────────────────────────────────────────────────────────
#  SGS curve helpers
# ─────────────────────────────────────────────────────────────────────────────
function get_sgs_curves(df, row::Int, n_bins::Int)
    eb  = collect(df[row, :sgs_elev_bins])
    vc  = collect(df[row, :sgs_vol_curve])
    ac  = collect(df[row, :sgs_area_curve])
    return eb[1:n_bins], vc[1:n_bins], ac[1:n_bins]
end

function get_edge_curves(df, row::Int, slot::Int, n_bins::Int)
    if !("sgs_edge_area_curve" in names(df))
        return nothing, nothing, nothing
    end
    ea_flat = collect(df[row, :sgs_edge_area_curve])
    ep_flat = collect(df[row, :sgs_edge_perim_curve])
    # flat array is n_bins*5, reshape to (n_bins, 5)
    ea = reshape(ea_flat[1:n_bins*5], n_bins, 5)
    ep = reshape(ep_flat[1:n_bins*5], n_bins, 5)
    return ea[:, slot], ep[:, slot], df[row, :sgs_edge_sills] isa AbstractVector ?
           df[row, :sgs_edge_sills][slot] :
           NaN
end

function get_edge_sill(df, row::Int, slot::Int)
    col = "sgs_edge_sills"
    col in names(df) || return NaN
    v = df[row, col]
    v isa AbstractVector && return Float64(v[slot])
    return NaN
end

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 1 — Wire-frame mesh diagram
# ─────────────────────────────────────────────────────────────────────────────
function fig_wireframe(df, idx, adj, target_id, outdir, dpi)
    @info "Drawing wire-frame mesh diagram..."
    row0 = idx[target_id]

    # Collect target + all immediate neighbours
    nbs   = get(adj, target_id, String[])
    cells = [target_id; nbs]

    # Find centroid of the group for local projection
    clon = mean(df[idx[c], :center_lon] for c in cells)
    clat = mean(df[idx[c], :center_lat] for c in cells)

    # Approximate elevations for flow-direction arrows
    # (use sgs z_min if available, else mean elevation)
    get_elev(c) = begin
        r = idx[c]
        if "sgs_z_min" in names(df) && !isnan(df[r, :sgs_z_min])
            return Float64(df[r, :sgs_z_min])
        elseif "elevation" in names(df) && !isnan(df[r, :elevation])
            return Float64(df[r, :elevation])
        else
            return 0.0
        end
    end
    elevs = Dict(c => get_elev(c) for c in cells)

    fig = Figure(size=(700, 700), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        aspect=DataAspect(),
        xlabel="Easting (m)", ylabel="Northing (m)",
        title="A5 mesh structure: cell $(target_id[1:8])…\nFlow arrows from higher to lower cell",
        titlesize=13)
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false)

    # Colour palette: target = blue, neighbours = light grey
    for (ci, cid) in enumerate(cells)
        r     = idx[cid]
        bnd   = get_boundary(df, r)
        xs, ys = to_local_m(bnd, clon, clat)
        col   = cid == target_id ? RGBAf(0.20, 0.45, 0.75, 0.25) :
                                   RGBAf(0.80, 0.83, 0.88, 0.30)
        border = cid == target_id ? RGBf(0.10, 0.30, 0.65) : RGBf(0.45, 0.50, 0.60)
        poly!(ax, Point2f.(zip(xs, ys)), color=col,
              strokecolor=border, strokewidth= cid==target_id ? 2.0 : 1.2)
    end

    # Cell ID labels (shortened) and elevation
    for cid in cells
        r    = idx[cid]
        cx   = (df[r, :center_lon] - clon) * deg2rad(1) * R_EARTH * cosd(clat)
        cy   = (df[r, :center_lat] - clat) * deg2rad(1) * R_EARTH
        elev = elevs[cid]
        lbl  = "$(cid[1:4])…$(cid[end-3:end])\nz=$(round(elev, digits=1)) m"
        fontsize = cid == target_id ? 10 : 8
        text!(ax, cx, cy, text=lbl, align=(:center,:center),
              fontsize=fontsize, color= cid==target_id ? :black : :gray40)
    end

    # Flow arrows: for each edge (target ↔ neighbour and neighbour ↔ neighbour)
    drawn_edges = Set{Tuple{String,String}}()
    for c1 in cells, c2 in cells
        c1 == c2 && continue
        edge_key = c1 < c2 ? (c1, c2) : (c2, c1)
        edge_key in drawn_edges && continue
        c2 in get(adj, c1, String[]) || continue
        push!(drawn_edges, edge_key)

        r1 = idx[c1]; r2 = idx[c2]
        cx1 = (df[r1,:center_lon]-clon)*deg2rad(1)*R_EARTH*cosd(clat)
        cy1 = (df[r1,:center_lat]-clat)*deg2rad(1)*R_EARTH
        cx2 = (df[r2,:center_lon]-clon)*deg2rad(1)*R_EARTH*cosd(clat)
        cy2 = (df[r2,:center_lat]-clat)*deg2rad(1)*R_EARTH

        e1, e2 = elevs[c1], elevs[c2]
        # Only draw arrow if there is a meaningful elevation difference (>0.01 m)
        abs(e1 - e2) < 0.01 && continue

        # Arrow from higher to lower
        src, dst = e1 > e2 ? (Point2f(cx1,cy1), Point2f(cx2,cy2)) :
                              (Point2f(cx2,cy2), Point2f(cx1,cy1))
        # Shorten arrow so it doesn't overlap cell centres
        mid   = 0.5f0*(src + dst)
        dir   = normalize(dst - src)
        shaft = 0.28f0*norm(dst - src)
        p0    = mid - dir*shaft
        p1    = mid + dir*shaft
        arrows2d!(ax, [p0[1]], [p0[2]], [p1[1]-p0[1]], [p1[2]-p0[2]],
                 tipwidth=10, shaftwidth=2.5,
                 color=RGBf(0.12, 0.55, 0.18),
                 tipcolor=RGBf(0.08, 0.40, 0.12))
    end

    # Legend
    elem_target = PolyElement(color=RGBAf(0.20,0.45,0.75,0.25),
                              strokecolor=RGBf(0.10,0.30,0.65), strokewidth=2)
    elem_nb     = PolyElement(color=RGBAf(0.80,0.83,0.88,0.30),
                              strokecolor=RGBf(0.45,0.50,0.60), strokewidth=1.2)
    elem_arrow  = LineElement(color=RGBf(0.08,0.40,0.12), linewidth=2)
    Legend(fig[1,2],
        [elem_target, elem_nb, elem_arrow],
        ["Target cell", "Neighbour cells", "Flow direction"],
        framevisible=false, labelsize=11)

    save_figure(fig, outdir, "fig1_wireframe", dpi)
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 2 — Edge-sill cross-section
# ─────────────────────────────────────────────────────────────────────────────
function fig_edge_sill(df, idx, adj, target_id, slot, wse_arg, outdir, dpi)
    @info "Drawing edge-sill cross-section (slot $slot)..."
    row0 = idx[target_id]

    nbs = get(adj, target_id, String[])
    if slot > length(nbs)
        @warn "Cell only has $(length(nbs)) neighbours; using slot 1."
        slot = 1
    end
    nb_id = nbs[slot]

    # Get SGS tables
    has_sgs = "sgs_elev_bins" in names(df)
    if !has_sgs
        @warn "No SGS tables in parquet — skipping edge-sill figure."
        return nothing
    end

    n_bins = length(collect(df[row0, :sgs_elev_bins]))

    # Edge area / perim curves
    ea_flat = "sgs_edge_area_curve" in names(df) ?
              collect(df[row0, :sgs_edge_area_curve]) : Float64[]
    ep_flat = "sgs_edge_perim_curve" in names(df) ?
              collect(df[row0, :sgs_edge_perim_curve]) : Float64[]

    has_ra = length(ea_flat) >= n_bins * 5

    elev_bins = collect(df[row0, :sgs_elev_bins])

    # z_sill from edge_sills column
    z_sill = NaN
    if "sgs_edge_sills" in names(df)
        sv = collect(df[row0, :sgs_edge_sills])
        slot <= length(sv) && (z_sill = Float64(sv[slot]))
    end

    z_min = "sgs_z_min" in names(df) ? Float64(df[row0, :sgs_z_min]) : elev_bins[1]
    z_max = "sgs_z_max" in names(df) ? Float64(df[row0, :sgs_z_max]) : elev_bins[end]

    # Representative WSE for annotation
    wse = isnan(wse_arg) ? z_sill + 0.5*(z_max - z_sill) : wse_arg
    wse = clamp(wse, z_sill, z_max)

    # Synthesise a plausible cross-section terrain profile from the edge curves
    # We use the edge_area_curve to back-calculate an approximate stepped profile.
    # If R-A curves absent, fall back to a simple V-shape centred on z_sill.
    n_pts = 60
    if has_ra
        ea = reshape(ea_flat[1:n_bins*5], n_bins, 5)[:, slot]
        # Differentiate area to get incremental width at each elevation bin
        widths = diff([0.0; ea]) ./ max.(diff([z_sill; elev_bins[1:n_bins]]), 1e-6)
        # Build stepped terrain profile: positions across edge width
        # Total edge width = last ea / (wse_bin - z_sill) approx
        edge_w = "sgs_cell_area" in names(df) ?
                 sqrt(Float64(df[row0, :sgs_cell_area])) : 100.0
        xs_terrain = range(-edge_w/2, edge_w/2, length=n_pts)
        # Simple triangular channel centred at 0, reaching z_sill at 0 and
        # rising to approximate z_max at edges — more realistic than flat
        chan_half = edge_w * 0.18   # channel occupies ~36% of edge width
        ys_terrain = [abs(x) <= chan_half ?
                      z_sill + (z_max - z_sill) * (abs(x)/chan_half)^1.5 :
                      z_max - 0.05*(z_max - z_sill) * rand()
                      for x in xs_terrain]
    else
        edge_w = "sgs_cell_area" in names(df) ?
                 sqrt(Float64(df[row0, :sgs_cell_area])) : 100.0
        xs_terrain = range(-edge_w/2, edge_w/2, length=n_pts)
        chan_half  = edge_w * 0.20
        ys_terrain = [z_sill + (z_max - z_sill) * min(1.0, (abs(x)/chan_half)^1.3)
                      for x in xs_terrain]
    end

    # Compute A and P at the representative WSE by integrating terrain profile
    wet_xs = [x for (x,y) in zip(xs_terrain, ys_terrain) if y <= wse]
    wet_ys = [y for y in ys_terrain if y <= wse]
    A_demo = isempty(wet_xs) ? 0.0 :
             (step(xs_terrain)) * sum(max(0.0, wse - y) for y in wet_ys)
    P_demo = isempty(wet_xs) ? 0.0 :
             (step(xs_terrain)) * length(wet_xs)
    R_demo = P_demo > 0 ? A_demo / P_demo : 0.0

    fig = Figure(size=(700, 480), backgroundcolor=:white)
    ax  = Axis(fig[1,1],
        xlabel="Distance across edge (m)", ylabel="Elevation (m)",
        title="Edge cross-section: slot $slot  |  " *
              "$(target_id[1:8])… ↔ $(nb_id[1:8])…",
        titlesize=12)

    # Terrain fill below profile
    xs_v = collect(xs_terrain)
    ys_v = collect(ys_terrain)
    # Band from bottom to terrain
    y_bot = z_sill - 0.15*(z_max - z_sill)
    band_xs = vcat(xs_v, reverse(xs_v))
    band_ys = vcat(ys_v, fill(y_bot, n_pts))
    poly!(ax, Point2f.(zip(band_xs, band_ys)),
          color=RGBAf(0.62, 0.47, 0.35, 0.85), strokewidth=0)

    # Terrain line
    lines!(ax, xs_v, ys_v, color=RGBf(0.35, 0.22, 0.10), linewidth=2.5)

    # Water fill
    wet_mask = ys_v .<= wse
    if any(wet_mask)
        wxs = xs_v[wet_mask]
        wys = min.(ys_v[wet_mask], wse)
        band_wxs = vcat([wxs[1]], wxs, [wxs[end]], reverse(wxs))
        band_wys = vcat([wse],    wys, [wse],       fill(wse, length(wxs)))
        poly!(ax, Point2f.(zip(
                vcat(wxs, reverse(wxs)),
                vcat(fill(wse, length(wxs)), reverse(wys)))),
              color=RGBAf(0.20, 0.50, 0.80, 0.40), strokewidth=0)
        # Water surface line
        hlines!(ax, [wse], color=RGBf(0.10, 0.35, 0.70), linewidth=1.8,
                linestyle=:dash)
    end

    # z_sill line
    hlines!(ax, [z_sill], color=RGBf(0.75, 0.20, 0.10), linewidth=1.8,
            linestyle=:dot)

    # Annotations
    x_ann = xs_v[1] + 0.05*(xs_v[end] - xs_v[1])
    text!(ax, x_ann, z_sill + 0.02*(z_max-z_sill),
          text="z_sill = $(@sprintf("%.2f", z_sill)) m",
          fontsize=10, color=RGBf(0.75,0.20,0.10))
    text!(ax, x_ann, wse + 0.03*(z_max-z_sill),
          text="WSE = $(@sprintf("%.2f", wse)) m",
          fontsize=10, color=RGBf(0.10,0.35,0.70))

    # A and P labels with braces
    if any(wet_mask) && length(wet_xs) > 2
        mid_x = mean(wet_xs)
        arrows2d!(ax, [mid_x, mid_x], [z_sill, z_sill],
                  [0.0, 0.0], [wse - z_sill, -(wse - z_sill)*0.02],
                  tipwidth=8, shaftwidth=1.5, color=:gray40)
        text!(ax, mid_x + 0.02*(xs_v[end]-xs_v[1]), (z_sill+wse)/2,
              text="A = $(@sprintf("%.1f", A_demo)) m²",
              fontsize=10, color=:gray30)
        # P annotation: bracket along water surface extent
        text!(ax, wet_xs[1] + 0.02*(wet_xs[end]-wet_xs[1]),
              wse + 0.06*(z_max - z_sill),
              text="P = $(@sprintf("%.1f", P_demo)) m  →  R = $(@sprintf("%.2f", R_demo)) m",
              fontsize=10, color=:gray30)
    end

    # Legend
    elem_terr  = PolyElement(color=RGBAf(0.62,0.47,0.35,0.85))
    elem_water = PolyElement(color=RGBAf(0.20,0.50,0.80,0.40))
    elem_sill  = LineElement(color=RGBf(0.75,0.20,0.10), linewidth=2, linestyle=:dot)
    elem_wse   = LineElement(color=RGBf(0.10,0.35,0.70), linewidth=2, linestyle=:dash)
    Legend(fig[1,2],
        [elem_terr, elem_water, elem_sill, elem_wse],
        ["Terrain", "Water body", "z_sill", "WSE"],
        framevisible=false, labelsize=11)

    save_figure(fig, outdir, "fig2_edge_sill", dpi)
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 3 — Hypsometric curve
# ─────────────────────────────────────────────────────────────────────────────
function fig_hypsometric(df, idx, target_id, wse_arg, outdir, dpi)
    @info "Drawing hypsometric curve..."
    row0 = idx[target_id]

    has_sgs = "sgs_elev_bins" in names(df)
    if !has_sgs
        @warn "No SGS tables in parquet — skipping hypsometric figure."
        return nothing
    end

    elev_bins = Float64.(collect(df[row0, :sgs_elev_bins]))
    vol_curve = Float64.(collect(df[row0, :sgs_vol_curve]))
    area_curve = "sgs_area_curve" in names(df) ?
                 Float64.(collect(df[row0, :sgs_area_curve])) : nothing

    n_bins = length(elev_bins)
    z_min  = "sgs_z_min" in names(df) ? Float64(df[row0, :sgs_z_min]) : elev_bins[1]
    z_max  = "sgs_z_max" in names(df) ? Float64(df[row0, :sgs_z_max]) : elev_bins[end]
    z_sill = NaN
    if "sgs_edge_sills" in names(df)
        sv = Float64.(collect(df[row0, :sgs_edge_sills]))
        z_sill = minimum(filter(isfinite, sv); init=NaN)
    end

    wse = isnan(wse_arg) ? z_min + 0.6*(z_max - z_min) : wse_arg
    wse = clamp(wse, z_min, z_max)

    # Interpolate volume at representative WSE
    vol_at_wse = let
        i = searchsortedfirst(elev_bins, wse)
        if i <= 1
            vol_curve[1]
        elseif i > n_bins
            vol_curve[end]
        else
            t = (wse - elev_bins[i-1]) / (elev_bins[i] - elev_bins[i-1])
            vol_curve[i-1] + t*(vol_curve[i] - vol_curve[i-1])
        end
    end

    fig = Figure(size=(640, 560), backgroundcolor=:white)

    # Primary axis: elevation vs volume
    ax = Axis(fig[1,1],
        xlabel="Cumulative volume V (m³)",
        ylabel="Elevation z (m)",
        title="Hypsometric curve: cell $(target_id[1:8])…",
        titlesize=12)

    # Shaded region below WSE
    idx_wet = findall(elev_bins .<= wse)
    if !isempty(idx_wet)
        ev_wet = [elev_bins[idx_wet]; wse]
        vv_wet = [vol_curve[idx_wet]; vol_at_wse]
        band!(ax, vv_wet, fill(z_min - 0.05*(z_max-z_min), length(vv_wet)), ev_wet,
              color=RGBAf(0.20, 0.50, 0.80, 0.20))
    end

    # Full curve
    lines!(ax, vol_curve, elev_bins,
           color=RGBf(0.15, 0.35, 0.70), linewidth=2.5, label="V(z)")

    # Annotation lines
    # z_min
    hlines!(ax, [z_min], color=:gray60, linewidth=1.2, linestyle=:dash)
    text!(ax, vol_curve[end]*0.02, z_min + 0.01*(z_max-z_min),
          text="z_min = $(@sprintf("%.2f", z_min)) m",
          fontsize=10, color=:gray50)

    # z_max
    hlines!(ax, [z_max], color=:gray40, linewidth=1.2, linestyle=:dash)
    text!(ax, vol_curve[end]*0.02, z_max + 0.01*(z_max-z_min),
          text="z_max = $(@sprintf("%.2f", z_max)) m",
          fontsize=10, color=:gray40)

    # z_sill
    if isfinite(z_sill)
        hlines!(ax, [z_sill], color=RGBf(0.75,0.20,0.10),
                linewidth=1.8, linestyle=:dot)
        text!(ax, vol_curve[end]*0.55, z_sill + 0.015*(z_max-z_min),
              text="z_sill (min) = $(@sprintf("%.2f", z_sill)) m",
              fontsize=10, color=RGBf(0.75,0.20,0.10))
    end

    # WSE annotation
    hlines!(ax, [wse], color=RGBf(0.10, 0.35, 0.70), linewidth=1.8, linestyle=:dash)
    vlines!(ax, [vol_at_wse], color=RGBf(0.10, 0.35, 0.70), linewidth=1.2,
            linestyle=:dash)
    scatter!(ax, [vol_at_wse], [wse], color=RGBf(0.10,0.35,0.70),
             markersize=10, marker=:circle)
    text!(ax, vol_at_wse + vol_curve[end]*0.01, wse + 0.015*(z_max-z_min),
          text="WSE = $(@sprintf("%.2f", wse)) m\nV = $(@sprintf("%.1f", vol_at_wse)) m³",
          fontsize=10, color=RGBf(0.10,0.35,0.70))

    # Secondary axis: wetted area if available
    if area_curve !== nothing
        ax2 = Axis(fig[1,1],
            xaxisposition=:top,
            xlabel="Wetted area Aᵥᵥ (m²)",
            xlabelcolor=RGBf(0.20,0.60,0.20),
            xticklabelcolor=RGBf(0.20,0.60,0.20),
            ygridvisible=false, yticksvisible=false,
            yticklabelsvisible=false)
        lines!(ax2, Float64.(area_curve), elev_bins,
               color=RGBf(0.20,0.60,0.20), linewidth=1.8,
               linestyle=:dash, label="Aᵥᵥ(z)")
        hidespines!(ax2, :l, :r, :b)
    end

    # Legend
    elems = Any[LineElement(color=RGBf(0.15,0.35,0.70), linewidth=2.5)]
    labs  = ["V(z) — cumulative volume"]
    if area_curve !== nothing
        push!(elems, LineElement(color=RGBf(0.20,0.60,0.20), linewidth=1.8,
                                 linestyle=:dash))
        push!(labs, "Aᵥᵥ(z) — wetted area")
    end
    push!(elems, PolyElement(color=RGBAf(0.20,0.50,0.80,0.20)))
    push!(labs, "Submerged at WSE")
    Legend(fig[2,1], elems, labs, orientation=:horizontal,
           framevisible=false, labelsize=11)

    save_figure(fig, outdir, "fig3_hypsometric", dpi)
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
#  Save helper (SVG + PNG)
# ─────────────────────────────────────────────────────────────────────────────
function save_figure(fig, outdir, name, dpi)
    mkpath(outdir)
    svg_path = joinpath(outdir, "$name.svg")
    png_path = joinpath(outdir, "$name.png")
    save(svg_path, fig)
    save(png_path, fig, px_per_unit = dpi / 96.0)
    @info "  Saved: $svg_path"
    @info "  Saved: $png_path"
end

# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────
function main()
    args = parse_args_()

    mesh_path = args["mesh"]
    # Strip any trailing whitespace, backslashes or quotes that PowerShell
    # may append when using \ as a mistaken line continuation character
    raw_cell = strip(args["cell"], [' ', '\\', '/', '"', '\'', '\t', '\n', '\r'])
    cell_id  = lpad(string(parse(UInt64, raw_cell, base=16), base=16), 16, '0')
    outdir    = args["outdir"]
    slot      = args["slot"]
    wse_arg   = args["wse"]
    dpi       = args["dpi"]

    @info "Loading parquet: $mesh_path"
    df, idx, adj, has_sgs = load_parquet(mesh_path)
    @info "  Loaded $(nrow(df)) cells.  SGS tables: $has_sgs"

    haskey(idx, cell_id) ||
        error("Cell $cell_id not found in mesh.  Available IDs start with: " *
              join([c[1:8]*"…" for c in first(collect(keys(idx)), 5)], ", "))

    nbs = get(adj, cell_id, String[])
    @info "  Target cell: $cell_id  |  $(length(nbs)) neighbours: $(join([n[1:8]*"…" for n in nbs], ", "))"

    fig_wireframe(df, idx, adj, cell_id, outdir, dpi)
    fig_edge_sill(df, idx, adj, cell_id, slot, wse_arg, outdir, dpi)
    fig_hypsometric(df, idx, cell_id, wse_arg, outdir, dpi)

    @info "Done.  Output directory: $(abspath(outdir))"
end

main()
