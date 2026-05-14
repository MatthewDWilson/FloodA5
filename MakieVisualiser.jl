# MakieVisualiser.jl
# ------------------
# Native desktop flood visualiser using GLMakie.
#
# Layout (1200 × 960 window)
# --------------------------
#
#   ┌────────────────────────────────────┬───────────────────┐
#   │                                    │  DIAGNOSTICS      │
#   │   Map panel                        │  Sim time / step  │
#   │   (A5 pentagons, variable colour)  │  dt / frame       │
#   │   Source cells: bright stroke      │  Active variable  │
#   │                                    │  Max / wet stats  │
#   │                                    │  Mass balance     │
#   │                                    │  Mesh info        │
#   │                                    │  Lon/lat bounds   │
#   │                                    │  Threads / clock  │
#   ├────────────────────────────────────┤                   │
#   │  colorbar ───────────────────────  │                   │
#   ├────────────────────────────────────┤                   │
#   │  Display variable:  [Depth ▾]      │                   │
#   ├────────────┬──────────┬────────────┴──────────────────-┤
#   │  Volume    │  Mass    │  Wet cells │  Ring vol (bars)  │
#   │  budget    │  balance │            │  (test mode only) │
#   └────────────┴──────────┴────────────┴───────────────────┘
#
# Ring-index feature (test mode)
# --------------------------------
# Pass `source_indices` and `adjacency` to `start`.  A BFS from the source
# cells assigns each mesh cell a ring index (0 = source, 1 = first
# neighbours, etc., -1 = unreachable isolated cells).
#
# Two extra diagnostics become available:
#   • "Ring index" appears as an additional option in the Display Variable
#     menu.  Selecting it overlays the BFS ring number as the polygon colour
#     (fixed categorical colormap, fixed range).
#   • A fourth time-series panel shows a bar chart of the total water volume
#     per ring at the current frame, making the concentric flood cascade
#     immediately visible.  This panel is hidden when ring data is absent
#     (dense real-domain runs without source_indices).
#
# The ring-volume bar chart is useful on the 61-cell flat test mesh: you see
# ring 0 fill rapidly, then ring 1 (5 cells), ring 2 (11 cells), etc.
# On a real topographic mesh the rings are present but bar heights are
# uneven because topography modulates the cascade.  On very dense meshes
# (thousands of rings) you may want to pass `max_ring_display` to cap the
# number of bars shown.
#
# push_frame! keyword arguments
# --------------------------------
#   vol_added    — cumulative volume injected since t=0 (m³)
#   vol_domain   — current total volume in domain (m³)
#   vol_removed  — cumulative volume removed (m³); 0 until Phase 2 outflow BCs
#   n_wet        — count of cells with water_depth > 1e-4 m
#   sim_step     — integer simulation step counter
#   sim_dt       — last adaptive dt (s)
#
# Usage
# -----
#   include("MakieVisualiser.jl")
#   using .MakieVisualiser
#
#   vis = MakieVisualiser.start(mesh;
#             source_indices  = [idx1, idx2],
#             adjacency       = flow_state.adjacency)   # enables ring mode
#   MakieVisualiser.push_frame!(vis, cell_ids, depths, saturations, volumes,
#                                              velocities, t;
#                                              vol_added=…, vol_domain=…,
#                                              vol_removed=0.0,
#                                              n_wet=…, sim_step=…, sim_dt=…)
#   MakieVisualiser.stop(vis)

module MakieVisualiser

using GLMakie
using Printf
using Dates

export MakieVis, start, stop, push_frame!

# ---------------------------------------------------------------------------
# Colour palette (dark theme)
# ---------------------------------------------------------------------------

const BG_DARK   = RGBf(0.08, 0.10, 0.12)
const BG_MAP    = RGBf(0.05, 0.07, 0.10)
const BG_SIDE   = RGBf(0.07, 0.09, 0.12)
const BG_PLOT   = RGBf(0.06, 0.08, 0.11)
const COL_TEXT  = RGBf(0.75, 0.88, 0.92)
const COL_LABEL = RGBf(0.60, 0.75, 0.82)
const COL_GRID  = RGBAf(1.0, 1.0, 1.0, 0.06)
const COL_EDGE  = RGBAf(0.15, 0.25, 0.35, 0.55)

# Source-cell highlight stroke
const COL_SOURCE_STROKE = RGBAf(1.0, 0.92, 0.20, 0.90)   # bright amber-yellow
const SOURCE_STROKE_W   = 1.8f0

# Time-series line colours
const COL_ADDED   = RGBf(0.35, 0.85, 0.55)
const COL_DOMAIN  = RGBf(0.35, 0.65, 0.95)
const COL_REMOVED = RGBf(0.95, 0.45, 0.35)
const COL_MBERR   = RGBf(0.95, 0.80, 0.25)
const COL_WET     = RGBf(0.30, 0.90, 0.95)

# Ring index display
# Unreachable cells (ring == -1) mapped to a very dark colour so they
# visually recede — they're plotted at value 0 in the ring colormap.
const RING_COLORMAP = :viridis

# Display variable descriptors — extended with "Ring index" when ring data present
const BASE_DISPLAY_VARS = [
    ("Depth (m)",        :turbo,  "Water depth (m)"),
    ("Saturation (0–1)", :Blues,  "Wetted area fraction"),
    ("Volume (m³)",      :turbo,  "Stored volume (m³)"),
    ("Velocity (m/s)",   :plasma, "Velocity (m/s)"),
]
const RING_VAR = ("Ring index", RING_COLORMAP, "BFS ring from source")

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

mutable struct MakieVis
    # Fixed mesh geometry
    cell_ids   :: Vector{String}
    n_cells    :: Int
    resolution :: Int
    lon_range  :: Tuple{Float64, Float64}
    lat_range  :: Tuple{Float64, Float64}

    # Per-variable data Observables (aligned to mesh cell order)
    data_depth :: Observable{Vector{Float64}}
    data_sat   :: Observable{Vector{Float64}}
    data_vol   :: Observable{Vector{Float64}}
    data_vel   :: Observable{Vector{Float64}}

    # Ring index data (fixed after construction; Float64 for colormap use)
    # Empty when ring mode is inactive.
    ring_index    :: Vector{Int}      # BFS distance, -1 = unreachable
    data_ring     :: Vector{Float64}  # ring_index cast to Float64 (0-based for colormap)
    n_rings       :: Int              # max ring + 1; 0 when inactive
    ring_cell_map :: Vector{Vector{Int}}  # ring_cell_map[r+1] = cell indices at ring r

    # Display Observables — drive polygon colours and colorbar
    display_data :: Observable{Vector{Float64}}
    colorrange   :: Observable{Tuple{Float64,Float64}}
    colormap_obs :: Observable{Symbol}
    cbar_label   :: Observable{String}
    selected_var :: Observable{Int}

    # All display-variable labels (may include "Ring index" at end)
    var_labels   :: Vector{String}
    var_cmaps    :: Vector{Symbol}
    var_units    :: Vector{String}

    # Simulation metadata
    sim_time    :: Observable{Float64}
    frame_idx   :: Observable{Int}
    diag_text   :: Observable{String}

    # Time-series history — plain vectors (grown each frame)
    _ts_time        :: Vector{Float64}
    _ts_vol_added   :: Vector{Float64}
    _ts_vol_domain  :: Vector{Float64}
    _ts_vol_removed :: Vector{Float64}
    _ts_mb_err      :: Vector{Float64}
    _ts_wet         :: Vector{Float64}

    # Observables that lines! binds to — reassigned (not mutated) each frame
    ts_time        :: Observable{Vector{Float64}}
    ts_vol_added   :: Observable{Vector{Float64}}
    ts_vol_domain  :: Observable{Vector{Float64}}
    ts_vol_removed :: Observable{Vector{Float64}}
    ts_mb_err      :: Observable{Vector{Float64}}
    ts_wet         :: Observable{Vector{Float64}}

    # Ring bar chart Observable — heights = volume per ring at current frame
    # Empty vector when ring mode is inactive (bar plot simply shows nothing)
    ring_vol_obs :: Observable{Vector{Float64}}

    # Time-series + ring axes
    ax_vol  :: Axis
    ax_mb   :: Axis
    ax_wet  :: Axis
    ax_ring :: Union{Axis, Nothing}   # Nothing when ring mode inactive

    figure  :: Figure
    running :: Threads.Atomic{Bool}
end

# ---------------------------------------------------------------------------
# BFS ring-index computation
# ---------------------------------------------------------------------------

"""
    _bfs_ring_index(cell_ids, source_indices, adjacency) → (ring_index, ring_cell_map)

Compute BFS ring distances from source cells.

Returns:
  ring_index    :: Vector{Int}  — BFS distance for each cell (0=source, -1=unreachable)
  ring_cell_map :: Vector{Vector{Int}} — ring_cell_map[r+1] = cell indices at ring r
  n_rings       :: Int — number of distinct rings (max ring + 1)

The adjacency dict maps cell_id → [neighbour_id, ...].  Cells not present as
keys (boundary cells with no outgoing neighbours) are treated as unreachable
beyond themselves.
"""
function _bfs_ring_index(cell_ids     :: Vector{String},
                          source_idx  :: Vector{Int},
                          adjacency   :: Dict{String, Vector{String}})

    n         = length(cell_ids)
    id_to_idx = Dict{String, Int}(id => i for (i, id) in enumerate(cell_ids))

    ring = fill(-1, n)
    queue = Int[]

    # Initialise BFS from all source cells simultaneously (multi-source BFS).
    # This gives "distance to nearest source" when there are multiple sources,
    # which is the right thing to display on the map.
    for src in source_idx
        if 1 <= src <= n && ring[src] == -1
            ring[src] = 0
            push!(queue, src)
        end
    end

    head = 1
    while head <= length(queue)
        ci = queue[head]
        head += 1
        cid = cell_ids[ci]
        neighbours = get(adjacency, cid, String[])
        for nid in neighbours
            ni = get(id_to_idx, nid, 0)
            ni == 0 && continue
            if ring[ni] == -1
                ring[ni] = ring[ci] + 1
                push!(queue, ni)
            end
        end
    end

    n_rings = isempty(queue) ? 0 : maximum(ring[i] for i in 1:n if ring[i] >= 0) + 1

    ring_cell_map = [Int[] for _ in 1:n_rings]
    for i in 1:n
        r = ring[i]
        r >= 0 && push!(ring_cell_map[r + 1], i)
    end

    return ring, ring_cell_map, n_rings
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    start(mesh;
          source_indices  = Int[],
          adjacency       = Dict{String,Vector{String}}(),
          max_ring_display = 40,
          title           = "FloodA5 — Makie Viewer") → MakieVis

Open a GLMakie desktop window pre-loaded with the mesh geometry.

**Ring mode** is activated when both `source_indices` and `adjacency` are
provided.  A BFS from the source cells assigns each cell a ring index.  This
adds a "Ring index" option to the display variable menu and a per-ring volume
bar chart in the bottom strip.

`max_ring_display` caps the number of rings shown in the bar chart (default
40).  On dense real-domain meshes the BFS may produce hundreds of rings; cap
to avoid an unreadable chart.  The map overlay always shows all rings.
"""
function start(mesh;
               source_indices   :: Vector{Int}                   = Int[],
               adjacency        :: Dict{String, Vector{String}}  = Dict{String,Vector{String}}(),
               max_ring_display :: Int                           = 40,
               title            :: String                        = "FloodA5 — Makie Viewer")

    n        = length(mesh.cells)
    cell_ids = [c.id for c in mesh.cells]
    lons     = [c.center_lon for c in mesh.cells]
    lats     = [c.center_lat for c in mesh.cells]

    polys = Vector{Vector{Point2f}}(undef, n)
    for (i, c) in enumerate(mesh.cells)
        ring = [Point2f(Float32(v[1]), Float32(v[2])) for v in c.boundary]
        first(ring) ≈ last(ring) || push!(ring, first(ring))
        polys[i] = ring
    end

    lon_range = extrema(lons)
    lat_range = extrema(lats)
    zeros_n   = zeros(Float64, n)

    # --- Ring index ----------------------------------------------------------
    ring_mode = !isempty(source_indices) && !isempty(adjacency)
    ring_index, ring_cell_map, n_rings = if ring_mode
        _bfs_ring_index(cell_ids, source_indices, adjacency)
    else
        (fill(-1, n), Vector{Vector{Int}}(), 0)
    end

    # Float64 representation for colormap: unreachable (-1) → 0.0, others as-is
    data_ring_f64 = Float64[max(r, 0) for r in ring_index]

    if ring_mode
        @info "Ring mode active: $(n_rings) rings from $(length(source_indices)) source cell(s)"
        for r in 0:min(n_rings-1, 5)
            @info "  Ring $(r): $(length(ring_cell_map[r+1])) cells"
        end
        n_rings > 6 && @info "  ... ($(n_rings) rings total)"
    end

    # --- Display variable list (with optional Ring index) --------------------
    base_vars  = BASE_DISPLAY_VARS
    all_vars   = ring_mode ? vcat(base_vars, [RING_VAR]) : base_vars
    var_labels = [v[1] for v in all_vars]
    var_cmaps  = [v[2] for v in all_vars]
    var_units  = [v[3] for v in all_vars]

    # --- Per-variable data Observables ---------------------------------------
    data_depth = Observable(copy(zeros_n))
    data_sat   = Observable(copy(zeros_n))
    data_vol   = Observable(copy(zeros_n))
    data_vel   = Observable(copy(zeros_n))

    # --- Display Observables -------------------------------------------------
    selected_var = Observable(1)
    display_data = Observable(copy(zeros_n))
    colorrange   = Observable((0.0, 1.0))
    colormap_obs = Observable(var_cmaps[1])
    cbar_label   = Observable(var_units[1])

    # --- Simulation metadata -------------------------------------------------
    sim_time  = Observable(0.0)
    frame_idx = Observable(0)
    diag_text = Observable(
        _build_diag(0.0, 0, 0, 0.0, var_labels[1], 0.0, 0, 0.0, 0.0, 0.0,
                    n, mesh.resolution, lon_range, lat_range)
    )

    # --- Time-series plain vectors -------------------------------------------
    _ts_time        = Float64[]
    _ts_vol_added   = Float64[]
    _ts_vol_domain  = Float64[]
    _ts_vol_removed = Float64[]
    _ts_mb_err      = Float64[]
    _ts_wet         = Float64[]

    ts_time        = Observable(Float64[])
    ts_vol_added   = Observable(Float64[])
    ts_vol_domain  = Observable(Float64[])
    ts_vol_removed = Observable(Float64[])
    ts_mb_err      = Observable(Float64[])
    ts_wet         = Observable(Float64[])

    # Ring bar chart Observable — n_display_rings heights
    n_display = min(n_rings, max_ring_display)
    ring_vol_obs = Observable(zeros(Float64, max(n_display, 1)))

    # --- Figure & layout (1200 × 960) ----------------------------------------
    fig = Figure(size = (1200, 960), backgroundcolor = BG_DARK)

    # Row 1: map axis
    ax_map = Axis(fig[1, 1];
        title           = @lift("t = " * _fmt_time($(sim_time))),
        titlecolor      = COL_TEXT,
        titlesize       = 14,
        xlabel          = "Longitude",
        ylabel          = "Latitude",
        xlabelcolor     = COL_LABEL,
        ylabelcolor     = COL_LABEL,
        xticklabelcolor = COL_LABEL,
        yticklabelcolor = COL_LABEL,
        backgroundcolor = BG_MAP,
        xgridcolor      = COL_GRID,
        ygridcolor      = COL_GRID,
        aspect          = DataAspect(),
    )

    # Base polygon layer — uniform thin stroke, colour-mapped fill.
    # strokewidth MUST be a scalar: when passed as a per-cell Vector, GLMakie
    # expands it into a per-vertex buffer and errors because the buffer length
    # (n_cells) doesn't match the expanded line-primitive vertex count.
    pp = poly!(ax_map, polys;
        color       = display_data,
        colormap    = colormap_obs,
        colorrange  = colorrange,
        strokewidth = 0.4f0,
        strokecolor = COL_EDGE,
    )

    # Source-cell highlight overlay — second poly! on top, transparent fill,
    # bright thick stroke.  Only the source polygons are drawn so the rest of
    # the mesh is unaffected.  strokewidth stays scalar; no buffer-length error.
    if !isempty(source_indices)
        src_polys = [polys[i] for i in source_indices if 1 <= i <= n]
        poly!(ax_map, src_polys;
            color       = RGBAf(0, 0, 0, 0),   # transparent fill
            strokewidth = SOURCE_STROKE_W,
            strokecolor = COL_SOURCE_STROKE,
        )
    end

    # Row 2: colorbar
    Colorbar(fig[2, 1], pp;
        label          = cbar_label,
        vertical       = false,
        flipaxis       = false,
        labelcolor     = COL_LABEL,
        tickcolor      = COL_LABEL,
        ticklabelcolor = COL_LABEL,
        height         = 18,
    )

    # Row 3: variable selector Menu
    menu_box = fig[3, 1] = GridLayout()
    Label(menu_box[1, 1], "Display variable:";
        color     = COL_LABEL,
        fontsize  = 13,
        halign    = :right,
        tellwidth = false,
    )
    menu = Menu(menu_box[1, 2];
        options   = var_labels,
        default   = var_labels[1],
        width     = 180,
        tellwidth = false,
    )

    # Row 4: bottom time-series strip (3 panels always + 1 ring panel if active)
    ts_grid    = fig[4, 1:2] = GridLayout()
    n_ts_cols  = ring_mode ? 4 : 3

    # Left — volume budget
    ax_vol = Axis(ts_grid[1, 1];
        title           = "Volume budget",
        titlecolor      = COL_TEXT,
        titlesize       = 12,
        xlabel          = "Sim time (s)",
        ylabel          = "Volume (m³)",
        xlabelcolor     = COL_LABEL,
        ylabelcolor     = COL_LABEL,
        xticklabelcolor = COL_LABEL,
        yticklabelcolor = COL_LABEL,
        backgroundcolor = BG_PLOT,
        xgridcolor      = COL_GRID,
        ygridcolor      = COL_GRID,
    )
    lines!(ax_vol, ts_time, ts_vol_added;
           color = COL_ADDED,   linewidth = 1.8, label = "Added")
    lines!(ax_vol, ts_time, ts_vol_domain;
           color = COL_DOMAIN,  linewidth = 1.8, label = "Domain")
    lines!(ax_vol, ts_time, ts_vol_removed;
           color = COL_REMOVED, linewidth = 1.8, label = "Removed")
    axislegend(ax_vol;
        position        = :lt,
        backgroundcolor = RGBAf(0.07, 0.09, 0.12, 0.85),
        labelcolor      = COL_TEXT,
        framecolor      = COL_LABEL,
        fontsize        = 10,
    )

    # Centre — mass-balance error
    ax_mb = Axis(ts_grid[1, 2];
        title           = "Mass balance error",
        titlecolor      = COL_TEXT,
        titlesize       = 12,
        xlabel          = "Sim time (s)",
        ylabel          = "Error (m³)",
        xlabelcolor     = COL_LABEL,
        ylabelcolor     = COL_LABEL,
        xticklabelcolor = COL_LABEL,
        yticklabelcolor = COL_LABEL,
        backgroundcolor = BG_PLOT,
        xgridcolor      = COL_GRID,
        ygridcolor      = COL_GRID,
    )
    hlines!(ax_mb, [0.0]; color = COL_LABEL, linewidth = 0.8, linestyle = :dash)
    lines!(ax_mb, ts_time, ts_mb_err;
           color = COL_MBERR, linewidth = 1.8)

    # Right of centre — wet-cell count
    ax_wet = Axis(ts_grid[1, 3];
        title           = "Wet cells  (depth > 0.1 mm)",
        titlecolor      = COL_TEXT,
        titlesize       = 12,
        xlabel          = "Sim time (s)",
        ylabel          = "Cell count",
        xlabelcolor     = COL_LABEL,
        ylabelcolor     = COL_LABEL,
        xticklabelcolor = COL_LABEL,
        yticklabelcolor = COL_LABEL,
        backgroundcolor = BG_PLOT,
        xgridcolor      = COL_GRID,
        ygridcolor      = COL_GRID,
    )
    lines!(ax_wet, ts_time, ts_wet;
           color = COL_WET, linewidth = 1.8)

    # Far right — ring volume bar chart (only when ring mode active)
    ax_ring = if ring_mode
        ring_xs = Float64.(0:n_display-1)   # ring labels 0, 1, 2, …

        ax = Axis(ts_grid[1, 4];
            title           = "Volume by ring  (rings 0–$(n_display-1))",
            titlecolor      = COL_TEXT,
            titlesize       = 12,
            xlabel          = "Ring index",
            ylabel          = "Volume (m³)",
            xlabelcolor     = COL_LABEL,
            ylabelcolor     = COL_LABEL,
            xticklabelcolor = COL_LABEL,
            yticklabelcolor = COL_LABEL,
            backgroundcolor = BG_PLOT,
            xgridcolor      = COL_GRID,
            ygridcolor      = COL_GRID,
            # Show a tick every 5 rings; always include ring 0
            xticks          = 0:min(5, max(n_display-1, 1)):n_display-1,
        )

        # barplot! with Observable heights.
        # color must be a numeric vector (mapped through colormap), not a
        # colormap symbol — passing color=:viridis is interpreted as a literal
        # colour name and errors.  Passing color=ring_xs lets GLMakie map each
        # bar's x-position through the colormap, giving a gradient across rings.
        barplot!(ax, ring_xs, ring_vol_obs;
            color      = ring_xs,
            colormap   = :viridis,
            colorrange = (0.0, max(Float64(n_display - 1), 1.0)),
            gap        = 0.1,
        )
        ax
    else
        nothing
    end

    # --- Sidebar (col 2, rows 1–3) -------------------------------------------
    ax_side = Axis(fig[1:3, 2];
        backgroundcolor    = BG_SIDE,
        leftspinevisible   = false,
        rightspinevisible  = false,
        topspinevisible    = false,
        bottomspinevisible = false,
        xgridvisible       = false,
        ygridvisible       = false,
    )
    hidedecorations!(ax_side)
    # "DejaVu Sans Mono" is bundled with Julia/Makie and rendered as a vector
    # outline through GLMakie's own text pipeline (same as axis tick labels),
    # giving sharp, crisp glyphs at all sizes.  "Courier New" is a GDI bitmap
    # font that blurs when composited into an OpenGL framebuffer.
    text!(ax_side, 0.06, 0.97;
        text      = diag_text,
        align     = (:left, :top),
        space     = :relative,
        font      = "DejaVu Sans Mono",
        fontsize  = 14.5,
        color     = COL_TEXT,
    )

    # --- Column / row sizing -------------------------------------------------
    colsize!(fig.layout, 1, Relative(0.74))
    colsize!(fig.layout, 2, Relative(0.26))
    rowsize!(fig.layout, 1, Relative(0.46))
    rowsize!(fig.layout, 2, Fixed(42))
    rowsize!(fig.layout, 3, Fixed(36))
    rowsize!(fig.layout, 4, Relative(0.28))

    # Equal-width time-series panels
    for c in 1:n_ts_cols
        colsize!(ts_grid, c, Relative(1.0 / n_ts_cols))
    end

    # --- Menu callback -------------------------------------------------------
    # We capture all the local variables needed; `ring_index` etc. are in scope.
    on(menu.selection) do label
        idx = findfirst(==(label), var_labels)
        idx === nothing && return
        selected_var[] = idx

        # Ring index is the last entry (only present in ring mode)
        is_ring = ring_mode && idx == length(var_labels)
        new_data = if is_ring
            data_ring_f64
        else
            _pick_data(idx, data_depth, data_sat, data_vol, data_vel)
        end

        display_data[] = new_data
        colorrange[]   = if is_ring
            (0.0, Float64(max(n_rings - 1, 1)))
        else
            _auto_range(new_data, idx)
        end
        colormap_obs[] = var_cmaps[idx]
        cbar_label[]   = var_units[idx]

        cur_added  = isempty(_ts_vol_added)  ? 0.0 : last(_ts_vol_added)
        cur_domain = isempty(_ts_vol_domain) ? 0.0 : last(_ts_vol_domain)
        cur_mb     = isempty(_ts_mb_err)     ? 0.0 : last(_ts_mb_err)
        cur_wet    = isempty(_ts_wet)        ? 0   : round(Int, last(_ts_wet))
        diag_text[] = _build_diag(
            sim_time[], frame_idx[], 0, 0.0, label,
            maximum(new_data; init=0.0), cur_wet,
            cur_added, cur_domain, cur_mb,
            n, mesh.resolution, lon_range, lat_range,
        )
    end

    display(fig)
    @info "MakieVisualiser: window open — $(n) cells, res $(mesh.resolution)" *
          (ring_mode ? "  [ring mode: $(n_rings) rings]" : "")

    return MakieVis(
        cell_ids, n, mesh.resolution, lon_range, lat_range,
        data_depth, data_sat, data_vol, data_vel,
        ring_index, data_ring_f64, n_rings, ring_cell_map,
        display_data, colorrange, colormap_obs, cbar_label, selected_var,
        var_labels, var_cmaps, var_units,
        sim_time, frame_idx, diag_text,
        _ts_time, _ts_vol_added, _ts_vol_domain, _ts_vol_removed, _ts_mb_err, _ts_wet,
        ts_time, ts_vol_added, ts_vol_domain, ts_vol_removed, ts_mb_err, ts_wet,
        ring_vol_obs,
        ax_vol, ax_mb, ax_wet, ax_ring,
        fig,
        Threads.Atomic{Bool}(true),
    )
end

"""
    push_frame!(vis, cell_ids, depths, saturations, volumes, velocities, t;
                vol_added=0.0, vol_domain=0.0, vol_removed=0.0,
                n_wet=0, sim_step=0, sim_dt=0.0)

Update the viewer with a new simulation timestep.

**Observable reassignment:** time-series Observables are reassigned with a
fresh copy each frame (not mutated in-place + notified).  GLMakie's `lines!`
compares by object identity; `push!`+`notify` on the same vector silently
no-ops.  `obs[] = copy(vec)` guarantees a new object and reliably redraws.
"""
function push_frame!(vis         :: MakieVis,
                     cell_ids    :: Vector{String},
                     depths      :: Vector{Float64},
                     saturations :: Vector{Float64},
                     volumes     :: Vector{Float64},
                     velocities  :: Vector{Float64},
                     t           :: Float64;
                     vol_added   :: Float64 = 0.0,
                     vol_domain  :: Float64 = 0.0,
                     vol_removed :: Float64 = 0.0,
                     n_wet       :: Int     = 0,
                     sim_step    :: Int     = 0,
                     sim_dt      :: Float64 = 0.0)
    vis.running[] || return

    # Re-align arrays to the fixed mesh cell order
    id_to_pos = Dict{String,Int}(id => i for (i, id) in enumerate(cell_ids))
    function _align(vals)
        out = zeros(Float64, vis.n_cells)
        for (mesh_i, mesh_id) in enumerate(vis.cell_ids)
            p = get(id_to_pos, mesh_id, 0)
            p > 0 && (out[mesh_i] = vals[p])
        end
        out
    end

    vis.data_depth[] = _align(depths)
    vis.data_sat[]   = _align(saturations)
    vis.data_vol[]   = _align(volumes)
    vis.data_vel[]   = _align(velocities)

    # Update the displayed variable
    idx      = vis.selected_var[]
    is_ring  = vis.n_rings > 0 && idx == length(vis.var_labels)
    new_data = if is_ring
        vis.data_ring   # ring index is static — only needs to be set once
    else
        _pick_data(idx, vis.data_depth, vis.data_sat, vis.data_vol, vis.data_vel)
    end
    vis.display_data[] = new_data
    is_ring || (vis.colorrange[] = _auto_range(new_data, idx))

    frame            = vis.frame_idx[] + 1
    vis.sim_time[]   = t
    vis.frame_idx[]  = frame
    mb_err           = vol_added - vol_domain - vol_removed

    vis.diag_text[] = _build_diag(
        t, frame, sim_step, sim_dt, vis.var_labels[idx],
        maximum(new_data; init=0.0), n_wet,
        vol_added, vol_domain, mb_err,
        vis.n_cells, vis.resolution,
        vis.lon_range, vis.lat_range,
    )

    # --- Grow plain history vectors then reassign Observables ----------------
    push!(vis._ts_time,        t)
    push!(vis._ts_vol_added,   vol_added)
    push!(vis._ts_vol_domain,  vol_domain)
    push!(vis._ts_vol_removed, vol_removed)
    push!(vis._ts_mb_err,      mb_err)
    push!(vis._ts_wet,         Float64(n_wet))

    # Fresh copies so GLMakie sees a new object and fires lines! listeners
    vis.ts_time[]        = copy(vis._ts_time)
    vis.ts_vol_added[]   = copy(vis._ts_vol_added)
    vis.ts_vol_domain[]  = copy(vis._ts_vol_domain)
    vis.ts_vol_removed[] = copy(vis._ts_vol_removed)
    vis.ts_mb_err[]      = copy(vis._ts_mb_err)
    vis.ts_wet[]         = copy(vis._ts_wet)

    # --- Ring bar chart update -----------------------------------------------
    if vis.n_rings > 0
        n_disp   = length(vis.ring_vol_obs[])
        vol_data = vis.data_vol[]
        new_ring_vol = zeros(Float64, n_disp)
        for r in 0:n_disp-1
            cells_in_ring = vis.ring_cell_map[r + 1]
            new_ring_vol[r + 1] = sum(vol_data[i] for i in cells_in_ring; init=0.0)
        end
        vis.ring_vol_obs[] = new_ring_vol   # fresh vector → barplot! redraws

        # Scale ring bar chart y-axis to the maximum bar height with headroom
        if vis.ax_ring !== nothing
            hi = max(maximum(new_ring_vol; init=0.0) * 1.15, 1.0)
            ylims!(vis.ax_ring, 0.0, hi)
        end
    end

    # --- Auto-scale time-series axes -----------------------------------------
    t_hi = max(t, 1.0)
    xlims!(vis.ax_vol, 0.0, t_hi)
    xlims!(vis.ax_mb,  0.0, t_hi)
    xlims!(vis.ax_wet, 0.0, t_hi)
    _autolimits_vol!(vis.ax_vol, vol_added, vol_domain, vol_removed)
    _autolimits_mb!(vis.ax_mb,  vis._ts_mb_err)
    _autolimits_wet!(vis.ax_wet, vis._ts_wet, vis.n_cells)
end

"""
    stop(vis)

Close the GLMakie window and mark the visualiser as stopped.
"""
function stop(vis::MakieVis)
    vis.running[] = false
    try; GLMakie.closeall(); catch; end
    @info "MakieVisualiser: window closed"
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _pick_data(idx, d_depth, d_sat, d_vol, d_vel)
    idx == 1 && return d_depth[]
    idx == 2 && return d_sat[]
    idx == 3 && return d_vol[]
    return d_vel[]
end

function _auto_range(data::Vector{Float64}, var_idx::Int)
    var_idx == 2 && return (0.0, 1.0)
    hi = maximum(data; init=0.0)
    hi = hi > 1e-6 ? hi : 1.0
    return (0.0, hi)
end

function _fmt_time(t::Float64)::String
    t < 60   && return @sprintf("%d s",     round(Int, t))
    t < 3600 && return @sprintf("%.1f min", t / 60)
    return              @sprintf("%.2f h",  t / 3600)
end

function _autolimits_vol!(ax, vol_added, vol_domain, vol_removed)
    hi = max(vol_added, vol_domain, vol_removed, 0.0)
    ylims!(ax, 0.0, hi > 0.0 ? hi * 1.10 : 1.0)
end

function _autolimits_mb!(ax, mb_err_vec::Vector{Float64})
    isempty(mb_err_vec) && return
    peak = maximum(abs, mb_err_vec)
    # Scale tightly to the actual error magnitude with 20% headroom.
    # No hard floor: if the error is genuinely 3e-10 the axis shows ±~3.6e-10
    # so any deviation from zero is visible.  Only fall back to a non-zero
    # range when the error is exactly zero (all values identical), to prevent
    # a degenerate flat axis at the very first frame.
    half = peak > 0.0 ? peak * 1.20 : 1e-10
    ylims!(ax, -half, half)
end

function _autolimits_wet!(ax, wet_vec::Vector{Float64}, n_cells::Int)
    hi = isempty(wet_vec) ? Float64(n_cells) : max(maximum(wet_vec), 1.0)
    ylims!(ax, 0.0, hi * 1.10)
end

function _build_diag(
    t          :: Float64,
    frame      :: Int,
    sim_step   :: Int,
    sim_dt     :: Float64,
    var_label  :: String,
    max_val    :: Float64,
    n_wet,
    vol_added  :: Float64,
    vol_domain :: Float64,
    mb_err     :: Float64,
    n_cells    :: Int,
    res        :: Int,
    lon_range  :: Tuple{Float64, Float64},
    lat_range  :: Tuple{Float64, Float64},
)::String
    s  = "─────────────────────"
    vl = length(var_label) > 13 ? var_label[1:13] : rpad(var_label, 13)
    lines = [
        " FloodA5  Diagnostics",
        " $s",
        @sprintf(" Sim time   %10s",     _fmt_time(t)),
        @sprintf(" Step       %10d",     sim_step),
        @sprintf(" dt         %10.2f s", sim_dt),
        @sprintf(" Frame      %10d",     frame),
        " $s",
        " Variable",
        "  $vl",
        @sprintf(" Max value  %10.4g", max_val),
        @sprintf(" Wet cells  %10d",   Int(n_wet)),
        " $s",
        " Mass Balance",
        @sprintf(" Vol added  %10.4g", vol_added),
        @sprintf(" Vol domain %10.4g", vol_domain),
        @sprintf(" MB error   %10.4g", mb_err),
        " $s",
        @sprintf(" Mesh cells %10d",   n_cells),
        @sprintf(" Resolution %10d",   res),
        " $s",
        " Lon range",
        @sprintf("  %.4f → %.4f", lon_range[1], lon_range[2]),
        " Lat range",
        @sprintf("  %.4f → %.4f", lat_range[1], lat_range[2]),
        " $s",
        @sprintf(" Threads    %10d",  Threads.nthreads()),
        @sprintf(" Updated    %s",    Dates.format(now(), "HH:MM:SS")),
    ]
    return join(lines, "\n")
end

end # module MakieVisualiser
