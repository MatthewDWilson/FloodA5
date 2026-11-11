# MakieVisualiser.jl
# ------------------
# Native desktop flood visualiser using GLMakie.
#
# Layout (1200 × 960 window)
# --------------------------
#
#   ┌────────────────────────────────────┬────────────────────┐
#   │  fig[1,1]                          │  fig[1,2]          │
#   │  Map panel  (ax_map)               │  DIAGNOSTICS       │
#   │  A5 pentagons, variable colour     │  (ax_side)         │
#   │  Source cells: bright stroke       │  Sim time/step/dt  │
#   │  Quiver overlay (linesegments! +   │  Mass balance      │
#   │    scatter!)                       │  Mesh info         │
#   │                                    ├────────────────────┤
#   ├────────────────────────────────────┤  fig[2,2]          │
#   │  fig[2,1]  colorbar               │  CONTROLS          │
#   ├────────────────────────────────────┤  (ctrl_grid)       │
#   │  fig[3,1]  Pause [||]  (toolbar)  │  Display var [▾]   │
#   ├────────────┬──────────┬────────────┤  ──────────────    │
#   │  fig[4,1:2] time-series strip      │  Flow arrows ○     │
#   │  Volume    │  Mass    │  Wet cells │  Scale   [────]    │
#   │  budget    │  balance │  (+Ring)   │  Stride  [────]    │
#   └────────────┴──────────┴────────────┴  Vel min [────]    │
#                                        └────────────────────┘
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
#
# Replay of an output file (same window)
# --------------------------------------
#   rd  = FloodA5Replay.open_replay("run.h5")        # visualisation/replay_data.jl
#   vis = MakieVisualiser.replay(rd; fps = 4)         # blocks until closed
# adds a frame slider and a speed slider (row 3; the strip moves to row 4);
# the Pause/Play button drives playback and the strip shows the whole run
# with a marker at the shown frame. Command line: visualisation/replay.jl.

module MakieVisualiser

using GLMakie
using Printf
using Dates
using Statistics

export MakieVis, start, stop, push_frame!, replay, replay_goto!

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

    # Velocity vector components (geographic, m/s) — used for quiver overlay
    data_vel_u :: Observable{Vector{Float64}}
    data_vel_v :: Observable{Vector{Float64}}

    # Quiver overlay state
    quiver_visible  :: Observable{Bool}
    quiver_segs_obs  :: Observable{Vector{Point2f}}  # linesegments shaft data
    quiver_tips_obs  :: Observable{Vector{Point2f}}  # scatter tip positions
    quiver_rots_obs  :: Observable{Vector{Float32}}  # scatter tip rotations
    quiver_sizes_obs :: Observable{Vector{Float32}}  # scatter markersize per arrow
    quiver_lw_obs    :: Observable{Float32}          # linesegments linewidth (scalar, peak-based)
    quiver_scale    :: Observable{Float64}  # arrow length multiplier (user slider)
    quiver_stride   :: Observable{Int}      # plot every N-th cell (density control)
    cell_lons       :: Vector{Float64}    # cell centre longitudes (mesh order)
    cell_lats       :: Vector{Float64}    # cell centre latitudes  (mesh order)
    paused          :: Threads.Atomic{Bool}  # true = sim loop is paused

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

    # The window's screen (`nothing` if not displayed) and, in replay mode,
    # the replay controls and data (see `replay`); `nothing` when live.
    screen    :: Any
    replay_ui :: Any
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
               title            :: String                        = "FloodA5 — Makie Viewer",
               show_window      :: Bool                          = true,
               _replay                                           = nothing)

    n        = length(mesh.cells)
    cell_ids = [c.id for c in mesh.cells]
    lons     = [c.center_lon for c in mesh.cells]
    lats     = [c.center_lat for c in mesh.cells]
    cell_lons_fixed = copy(lons)   # fixed mesh-order copy for quiver
    cell_lats_fixed = copy(lats)

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
    data_vel_u = Observable(copy(zeros_n))
    data_vel_v = Observable(copy(zeros_n))
    quiver_visible = Observable(false)

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

    # --- Figure & layout ----------------------------------------
    # (replay mode: 40 px taller for the frame-slider row)
    fig = Figure(size = (1200, _replay === nothing ? 800 : 840), backgroundcolor = BG_DARK)

    # Geographic limits: expand lon half-span by 1/cos(lat) so that one degree
    # of longitude occupies the same data-unit length as one degree of latitude.
    # DataAspect() then maps equal data units to equal pixels, giving a
    # Cartesian-correct display that honours window resize without distortion.
    # Compute extent from polygon boundaries, not cell centres,
    # so boundary cells are not clipped at the axis edge.
    all_boundary_lons = [Float64(v[1]) for c in mesh.cells for v in c.boundary]
    all_boundary_lats = [Float64(v[2]) for c in mesh.cells for v in c.boundary]
    bnd_lon_range = extrema(all_boundary_lons)
    bnd_lat_range = extrema(all_boundary_lats)

    _pad      = 0.04          # small cosmetic padding only — extent already covers boundaries
    _mean_lat = mean(lats)
    _cos_lat  = cos(deg2rad(_mean_lat))
    _lon_c    = (bnd_lon_range[1] + bnd_lon_range[2]) / 2
    _lat_c    = (bnd_lat_range[1] + bnd_lat_range[2]) / 2

    # Each half-span must (a) contain the data with padding, and (b) maintain
    # the correct geographic aspect so DataAspect() renders without distortion.
    # Take the max of both constraints in each axis.
    _lat_half_data = (bnd_lat_range[2] - bnd_lat_range[1]) / 2 * (1 + _pad)
    _lon_half_data = (bnd_lon_range[2] - bnd_lon_range[1]) / 2 * (1 + _pad)
    _lon_half = max(_lon_half_data, _lat_half_data / _cos_lat)
    _lat_half = max(_lat_half_data, _lon_half_data * _cos_lat)

    ax_map = Axis(fig[1, 1];
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
        limits          = (_lon_c - _lon_half, _lon_c + _lon_half,
                        _lat_c - _lat_half, _lat_c + _lat_half),
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

    # Quiver control Observables — defined before ctrl_grid widgets that use them
    quiver_scale   = Observable(1.0)
    quiver_stride  = Observable(1)
    _quiver_updating = Ref(false)

    # Helper: recompute quiver arrays with current settings
    function _refresh_quiver()
        _quiver_updating[] && return
        _quiver_updating[] = true
        try
            _update_quiver!(quiver_segs_obs, quiver_tips_obs, quiver_rots_obs, quiver_sizes_obs, quiver_lw_obs,
                            cell_lons_fixed, cell_lats_fixed,
                            data_vel_u[], data_vel_v[],
                            ax_map, quiver_scale[], quiver_stride[],
                            vel_thresh_obs[] * 1e-3)
        finally
            _quiver_updating[] = false
        end
    end

    # ── Row 4: bottom time-series strip ──────────────────────────────────────
    # Replay mode (`replay`): a frame-slider row goes above the strip
    is_replay  = _replay !== nothing
    ts_grid    = fig[is_replay ? 4 : 3, 1:2] = GridLayout()
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
           color = COL_REMOVED, linewidth = 1.8,
           label = is_replay ? _replay.removed_label : "Removed")
    axislegend(ax_vol;
        position        = :lt,
        backgroundcolor = RGBAf(0.07, 0.09, 0.12, 0.85),
        labelcolor      = COL_TEXT,
        framecolor      = COL_LABEL,
        fontsize        = 10,
    )

    # Centre — mass-balance error
    ax_mb = Axis(ts_grid[1, 2];
        title           = is_replay ? _replay.mb_title : "Mass balance error",
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

    # Time label: drawn inside the map axis at top-centre so it creates no
    # protrusion that would push the sidebar grid down.
    text!(ax_map, @lift("t = " * _fmt_time($(sim_time)));
        position  = (0.5, 0.99),
        align     = (:center, :top),
        space     = :relative,
        fontsize  = 14,
        color     = COL_TEXT,
    )

    # --- Sidebar col 2 --------------------------------------------------------
    # Right sidebar: a nested GridLayout spanning fig rows 1:2 (map+colorbar
    # height) so it is fully independent of the left column's row structure.
    # Internally: row 1 = controls (shrink-wraps), row 2 = diagnostics (fills).
    sidebar_grid = fig[1:2, 2] = GridLayout()

    # Controls panel — top of sidebar
    ctrl_grid = sidebar_grid[1, 1] = GridLayout()

    # Diagnostics panel — fills remaining sidebar height below controls
    ax_side = Axis(sidebar_grid[2, 1];
        backgroundcolor    = BG_SIDE,
        leftspinevisible   = false,
        rightspinevisible  = false,
        topspinevisible    = false,
        bottomspinevisible = false,
        xgridvisible       = false,
        ygridvisible       = false,
    )
    hidedecorations!(ax_side)
    text!(ax_side, 0.06, 0.97;
        text      = diag_text,
        align     = (:left, :top),
        space     = :relative,
        font      = "DejaVu Sans Mono",
        fontsize  = 11.0,
        color     = COL_TEXT,
    )

    # Controls row shrinks to content; diagnostics row fills the rest
    rowsize!(sidebar_grid, 1, Auto())       # controls: shrinks to content
    rowsize!(sidebar_grid, 2, Auto())       # diagnostics: fills remaining height
    rowgap!(sidebar_grid, 4)

    # Background box (Axis with no decorations gives a coloured background)
    ax_ctrl = Axis(ctrl_grid[1:8, 1:2];
        backgroundcolor    = BG_SIDE,
        leftspinevisible   = false,
        rightspinevisible  = false,
        topspinevisible    = false,
        bottomspinevisible = false,
        xgridvisible       = false,
        ygridvisible       = false,
    )
    hidedecorations!(ax_ctrl)

    # Title
    Label(ctrl_grid[1, 1:2], "Controls";
        color    = COL_TEXT,
        fontsize = 11,
        font     = "DejaVu Sans Mono",
        halign   = :left,
        tellwidth = false,
    )

    # Separator line via a thin coloured Label
    Label(ctrl_grid[2, 1:2], "─────────────────────────────";
        color    = COL_LABEL,
        fontsize = 10,
        halign   = :left,
        tellwidth = false,
    )

    # Row 3: Display variable selector
    Label(ctrl_grid[3, 1], "Display:";
        color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = false)
    menu = Menu(ctrl_grid[3, 2];
        options   = var_labels,
        default   = var_labels[1],
        width     = 130,
        tellwidth = false,
    )

    # Row 4: Flow arrows toggle
    Label(ctrl_grid[4, 1], "Flow arrows:";
        color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = false)
    quiver_toggle = Toggle(ctrl_grid[4, 2]; active = false, halign = :left, tellwidth = false)

    # Row 5: arrow scale — full widget-column slider
    Label(ctrl_grid[5, 1], "Arrow scale:";
        color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = false)
    quiver_scale_slider = Slider(ctrl_grid[5, 2];
        range = 0.1:0.1:5.0, startvalue = 1.0, width = 130, tellwidth = false)
    on(quiver_scale_slider.value) do v; quiver_scale[] = v; end

    # Row 6: arrow stride
    Label(ctrl_grid[6, 1], "Arrow stride:";
        color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = false)
    quiver_stride_slider = Slider(ctrl_grid[6, 2];
        range = 1:1:10, startvalue = 1, width = 130, tellwidth = false)
    on(quiver_stride_slider.value) do v; quiver_stride[] = v; end

    # Row 7: velocity threshold
    Label(ctrl_grid[7, 1], "Vel min (mm/s):";
        color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = false)
    vel_thresh_obs = Observable(0.01)
    vel_thresh_slider = Slider(ctrl_grid[7, 2];
        range = vcat(0.001f0, 0.01f0:0.01f0:0.1f0, 0.2f0:0.1f0:2.0f0),
        startvalue = 0.01f0, width = 130, tellwidth = false)
    on(vel_thresh_slider.value) do v; vel_thresh_obs[] = Float64(v); end

    # Row 8: pause / play
    Label(ctrl_grid[8, 1], is_replay ? "Playback:" : "Simulation:";
        color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = false)
    _sim_paused = Ref(false)
    _paused_obs = Observable(false)
    pause_btn = Button(ctrl_grid[8, 2];
        label     = @lift($(_paused_obs) ? "Play  [>]" : "Pause [||]"),
        width = 90, tellwidth = false)
    on(pause_btn.clicks) do _
        _sim_paused[] = !_sim_paused[]
        _paused_obs[] = _sim_paused[]
    end

    colsize!(ctrl_grid, 1, Auto())
    colsize!(ctrl_grid, 2, Auto())
    rowgap!(ctrl_grid, 3)

    # Wire quiver widgets → update function (all defined above)
    on(quiver_toggle.active) do val
        quiver_visible[] = val
        val && !isempty(data_vel_u[]) && _refresh_quiver()
    end
    on(quiver_scale)   do _; quiver_visible[] && _refresh_quiver(); end
    on(quiver_stride)  do _; quiver_visible[] && _refresh_quiver(); end
    on(vel_thresh_obs) do _; quiver_visible[] && _refresh_quiver(); end

    # ── Quiver overlay — built from linesegments! + scatter! ────────────────
    # We avoid arrows!/arrows2d! due to version-specific dispatch issues with
    # empty Observable arrays at construction time.  Instead each arrow is drawn
    # as a shaft (linesegments!) plus a filled triangle head (scatter! with
    # marker=:utriangle, rotated per-arrow via rotations observable).
    #
    # _update_quiver! populates three Observables every frame:
    #   quiver_segs_obs  — Vec2f pairs [tail, tip] for linesegments!
    #   quiver_tips_obs  — tip positions for scatter! arrowheads
    #   quiver_rots_obs  — rotation angle (radians) for each arrowhead
    #   quiver_sizes_obs — markersize (pixels) per arrowhead, scales with speed
    _domain_lon_span = lon_range[2] - lon_range[1]

    quiver_segs_obs  = Observable(Point2f[])   # interleaved [tail, tip, tail, tip …]
    quiver_tips_obs  = Observable(Point2f[])   # one per arrow
    quiver_rots_obs  = Observable(Float32[])   # one per arrow (radians)
    quiver_sizes_obs = Observable(Float32[])   # markersize per arrowhead
    quiver_lw_obs    = Observable(1.5f0)       # scalar linewidth, updated each frame

    # Capture plot objects so we can set .visible directly in the on() listener.
    # Passing an Observable as `visible=` at construction is not guaranteed to
    # create a live subscription in Makie 0.24 — we set it explicitly instead.
    quiver_lines = linesegments!(ax_map, quiver_segs_obs;
        color       = :white,
        linewidth   = quiver_lw_obs,
        visible     = false,
        inspectable = false,
    )
    quiver_heads = scatter!(ax_map, quiver_tips_obs;
        color        = :white,
        marker       = :utriangle,
        markersize   = quiver_sizes_obs,
        rotation     = quiver_rots_obs,
        visible      = false,
        inspectable  = false,
    )

    # Explicit visibility listener — directly sets the plot attribute.
    on(quiver_visible) do val
        quiver_lines.visible[] = val
        quiver_heads.visible[] = val
    end

    # Recompute filtered quiver arrays whenever vel_u changes.
    # vel_v.val is updated silently before vel_u[] is notified (see push_frame!).
    on(data_vel_u) do _
        _quiver_updating[] && return
        _quiver_updating[] = true
        try
            _update_quiver!(quiver_segs_obs, quiver_tips_obs, quiver_rots_obs, quiver_sizes_obs, quiver_lw_obs,
                            cell_lons_fixed, cell_lats_fixed,
                            data_vel_u[], data_vel_v[],
                            ax_map, quiver_scale[], quiver_stride[])
        finally
            _quiver_updating[] = false
        end
    end

    # --- Column / row sizing -------------------------------------------------
    # Col 1 (map+ts): fills all width not used by the fixed sidebar.
    # Col 2 (sidebar): Fixed 300px — controls and diagnostics never reflow.
    # Row 1 (map+sidebar): fills all height not used by fixed rows below.
    # Row 2 (colorbar): Fixed 50px — always the same thin strip.
    # Row 3 (time-series): Fixed 240px — keeps plots readable at any size.
    #
    # Result: resizing the window stretches the map and colorbar; the sidebar,
    # time-series panels, and colorbar strip stay a constant size.
    # Auto() = fill all space not claimed by Fixed columns/rows.
    # This is correct — Relative(1.0) means 100% of figure, overriding Fixed.
    # Follow Makie tutorial pattern: Auto() for scalable content.
    # Col 1 (map): Auto — expands to fill all width not taken by sidebar.
    # Col 2 (sidebar): Auto — sizes to its content (ctrl_grid + ax_side).
    # Row 1 (map+sidebar): Auto — fills all height not taken by fixed rows.
    # Row 2 (colorbar): Fixed 44px thin strip.
    # Row 3 (time-series): Fixed 200px — keeps plots readable at any window size.
    colsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 2, Fixed(280))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Fixed(44))
    rowsize!(fig.layout, 3, Fixed(200))
    rowgap!(fig.layout, 1, 4)
    rowgap!(fig.layout, 2, 4)

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

    # Replay mode: frame slider, frame/time readout and playback speed in a row
    # above the time-series strip, and a current-time marker on the strip.
    replay_ui = nothing
    if is_replay
        nfr = _replay.nframes
        rp_grid = fig[3, 1:2] = GridLayout()
        Label(rp_grid[1, 1], "Frame:";
            color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = true)
        frame_slider = Slider(rp_grid[1, 2]; range = 1:max(nfr, 1), startvalue = 1)
        frame_text = Observable("")
        Label(rp_grid[1, 3], frame_text;
            color = COL_TEXT, fontsize = 11, font = "DejaVu Sans Mono",
            halign = :left, width = 250, tellwidth = true)
        Label(rp_grid[1, 4], "Speed (frames/s):";
            color = COL_LABEL, fontsize = 11, halign = :right, tellwidth = true)
        fps_slider = Slider(rp_grid[1, 5];
            range = [0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 30.0],
            startvalue = 4.0, width = 140, tellwidth = true)
        Label(rp_grid[1, 6], lift(v -> @sprintf("%g", v), fps_slider.value);
            color = COL_TEXT, fontsize = 11, halign = :left, width = 30, tellwidth = true)
        colgap!(rp_grid, 8)
        rowsize!(fig.layout, 3, Fixed(30))
        rowsize!(fig.layout, 4, Fixed(200))
        rowgap!(fig.layout, 3, 4)
        t_marker = lift(t -> [t], sim_time)
        for ax in (ax_vol, ax_mb, ax_wet)
            vlines!(ax, t_marker; color = COL_LABEL, linewidth = 1.0, linestyle = :dot)
        end
        replay_ui = (slider = frame_slider, fps_slider = fps_slider, frame_text = frame_text,
                     menu = menu, pause_btn = pause_btn, paused_obs = _paused_obs,
                     nframes = nfr)
    end

    screen = show_window ? display(fig) : nothing
    @info "MakieVisualiser: window open — $(n) cells, res $(mesh.resolution)" *
          (ring_mode ? "  [ring mode: $(n_rings) rings]" : "")

    _vis = MakieVis(
        cell_ids, n, mesh.resolution, lon_range, lat_range,
        data_depth, data_sat, data_vol, data_vel,
        data_vel_u, data_vel_v,
        quiver_visible,
        quiver_segs_obs, quiver_tips_obs, quiver_rots_obs, quiver_sizes_obs, quiver_lw_obs,
        quiver_scale, quiver_stride,
        cell_lons_fixed, cell_lats_fixed,
        Threads.Atomic{Bool}(false),   # paused — starts unpaused
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
        screen, replay_ui,
    )
    # Wire _sim_paused Ref to vis.paused after construction.
    # The on(pause_btn.clicks) handler writes to _sim_paused[]; here we add
    # a second listener that mirrors it to vis.paused so the sim loop sees it.
    on(pause_btn.clicks) do _
        _vis.paused[] = _sim_paused[]
    end
    return _vis
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
                     vel_u       :: Vector{Float64},
                     vel_v       :: Vector{Float64},
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

    idx, new_data = _set_fields!(vis, _align(depths), _align(saturations), _align(volumes),
                                 _align(velocities), _align(vel_u), _align(vel_v))

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
    _update_ring_bars!(vis)

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
# Replay of an output file
# ---------------------------------------------------------------------------

"""
    replay(rd; fps = 4.0, var = 1, from = 1, playing = true,
           show_window = true, block = show_window, max_ring_display = 40) → MakieVis

Replay a run from its output file in the live viewer's window: same map,
colorbar, controls, diagnostics and time-series strip, plus a frame slider
(scrub to any output time), a playback speed slider (output frames per
second) and the Pause/Play button driving playback. `rd` is a
`FloodA5Replay.ReplayData` (`visualisation/replay_data.jl`; duck-typed:
`mesh`, `t`, `series`, `source_indices`, `adjacency`, `label`, `notes`,
`frame(k)`). `visualisation/replay.jl` is the command line.

- `var`: initial display variable (index or label, e.g. "Velocity (m/s)").
- `from`: first frame shown (index).
- The strip shows the whole run at once with a dotted marker at the shown
  frame; the volume series come from `rd.series` (see replay_data.jl for
  what is read from the file and what is rebuilt from the run record).
- `show_window = false` builds the figure without a window (for `save`);
  `block = true` returns only when the window is closed.
"""
function replay(rd; fps::Real = 4.0, var = 1, from::Integer = 1, playing::Bool = true,
                show_window::Bool = true, block::Bool = show_window,
                max_ring_display::Int = 40)
    t  = Float64.(rd.t)
    nf = length(t)
    nf > 0 || error("replay: no frames")
    s  = rd.series
    vis = start(rd.mesh;
                source_indices   = Vector{Int}(rd.source_indices),
                adjacency        = rd.adjacency,
                max_ring_display = max_ring_display,
                title            = "FloodA5 replay — " * rd.label,
                show_window      = show_window,
                _replay          = (nframes = nf, removed_label = s.removed_label,
                                    mb_title = s.mb_title))
    ui = merge(vis.replay_ui, (rd = rd, current = Ref(0)))
    vis.replay_ui = ui

    # Whole-run time series (fixed during replay)
    for (buf, obs, v) in ((vis._ts_time,        vis.ts_time,        t),
                          (vis._ts_vol_added,   vis.ts_vol_added,   s.vol_added),
                          (vis._ts_vol_domain,  vis.ts_vol_domain,  s.vol_domain),
                          (vis._ts_vol_removed, vis.ts_vol_removed, s.vol_removed),
                          (vis._ts_mb_err,      vis.ts_mb_err,      s.mb_err),
                          (vis._ts_wet,         vis.ts_wet,         s.n_wet))
        append!(buf, Float64.(v))
        obs[] = copy(buf)
    end
    fmax(v) = maximum((x for x in v if isfinite(x)); init = 0.0)
    # (2% headroom on the right when end-of-run markers sit at the last time)
    t_hi = max(t[end], 1.0) * (s.end_mark === nothing ? 1.0 : 1.02)
    for ax in (vis.ax_vol, vis.ax_mb, vis.ax_wet)
        xlims!(ax, 0.0, t_hi)
    end
    _autolimits_vol!(vis.ax_vol, fmax(s.vol_added), fmax(s.vol_domain), fmax(s.vol_removed))
    mbv = Float64[x for x in s.mb_err if isfinite(x)]
    if s.end_mark !== nothing
        te, rem, mb = s.end_mark
        scatter!(vis.ax_vol, [te], [rem]; color = COL_REMOVED, marker = :diamond,
                 markersize = 11, label = "Removed (end)")
        if isfinite(mb)
            scatter!(vis.ax_mb, [te], [mb]; color = COL_MBERR, marker = :diamond, markersize = 11)
            push!(mbv, mb)
        end
        _autolimits_vol!(vis.ax_vol, fmax(s.vol_added), fmax(s.vol_domain),
                         max(fmax(s.vol_removed), rem))
    end
    isempty(mbv) ? ylims!(vis.ax_mb, -1.0, 1.0) : _autolimits_mb!(vis.ax_mb, mbv)
    _autolimits_wet!(vis.ax_wet, Float64.(s.n_wet), vis.n_cells)

    # Frame display: the slider drives everything
    on(k -> _replay_show!(vis, k), ui.slider.value)
    on(ui.menu.selection) do _
        ui.current[] > 0 && _replay_diag!(vis, ui.current[])
    end
    vi = var isa Integer ? Int(var) : findfirst(==(String(var)), vis.var_labels)
    vi === nothing && error("replay: unknown variable $(repr(var)); one of $(vis.var_labels)")
    vi == 1 || (ui.menu.i_selected[] = vi)
    replay_goto!(vis, from)

    # Play/Pause: Play at the last frame starts again from the first
    on(ui.pause_btn.clicks) do _
        !vis.paused[] && ui.slider.value[] >= nf && replay_goto!(vis, 1)
    end
    playing || (ui.pause_btn.clicks[] = ui.pause_btn.clicks[] + 1)
    set_close_to!(ui.fps_slider, Float64(fps))

    foreach(n -> @info("Replay: " * n), rd.notes)
    if show_window
        if block
            try
                _replay_loop(vis)
            catch e
                e isa InterruptException || rethrow(e)
            finally
                stop(vis)
            end
        else
            @async _replay_loop(vis)
        end
    end
    return vis
end

"""
    replay_goto!(vis, k) → vis

Show frame `k` of a replay (moves the slider).
"""
function replay_goto!(vis::MakieVis, k::Integer)
    ui = vis.replay_ui
    ui === nothing && error("replay_goto!: not a replay window")
    k = clamp(Int(k), 1, ui.nframes)
    if ui.slider.value[] == k
        _replay_show!(vis, k)             # the slider would not move: show directly
    else
        set_close_to!(ui.slider, k)       # → on(slider.value) → _replay_show!
    end
    ui.current[] == k || _replay_show!(vis, k)
    return vis
end

# Playback: advance one frame per 1/fps s while not paused; pause at the end.
function _replay_loop(vis::MakieVis)
    ui = vis.replay_ui
    while vis.running[] && (vis.screen === nothing || isopen(vis.screen))
        if !vis.paused[]
            k = ui.slider.value[]
            if k >= ui.nframes
                ui.pause_btn.clicks[] = ui.pause_btn.clicks[] + 1
            else
                set_close_to!(ui.slider, k + 1)
            end
        end
        sleep(1.0 / ui.fps_slider.value[])
    end
end

function _replay_show!(vis::MakieVis, k::Integer)
    ui = vis.replay_ui
    rd = ui.rd
    fr = rd.frame(k)
    _set_fields!(vis, fr.depth, fr.saturation, fr.volume, fr.velocity, fr.vel_u, fr.vel_v)
    ui.current[]    = k
    vis.sim_time[]  = Float64(rd.t[k])
    vis.frame_idx[] = k
    _update_ring_bars!(vis)
    ui.frame_text[] = @sprintf("%d / %d   t = %s", k, ui.nframes, _fmt_time(Float64(rd.t[k])))
    _replay_diag!(vis, k)
    return nothing
end

function _replay_diag!(vis::MakieVis, k::Integer)
    ui  = vis.replay_ui
    rd  = ui.rd
    s   = rd.series
    idx = vis.selected_var[]
    shown = vis.display_data[]
    mx  = maximum((x for x in shown if isfinite(x)); init = 0.0)
    vis.diag_text[] = _build_replay_diag(
        Float64(rd.t[k]), k, ui.nframes, rd.label, vis.var_labels[idx], mx,
        Int(round(s.n_wet[k])), s.vol_added[k], s.vol_domain[k], s.vol_removed[k],
        s.mb_err[k], s.removed_label, s.inputs_known, s.removed,
        vis.n_cells, vis.resolution)
end

function _fmt_val(x::Float64)
    isfinite(x) || return @sprintf("%10s", "n/a")
    return @sprintf("%10.4g", x)
end

function _build_replay_diag(t, k, nf, label, var_label, max_val, n_wet,
                            vol_added, vol_domain, vol_removed, mb_err,
                            removed_label, inputs_known, removed, n_cells, res)::String
    s   = "─────────────────────"
    vl  = length(var_label) > 13 ? var_label[1:13] : rpad(var_label, 13)
    fl  = length(label) > 19 ? label[1:16] * "…" : label
    rl  = removed === :closed ? " Removed   " : " Rem+MB err"
    out = removed === :closed ? "closed" : removed === :end_only ? "end only" : "unknown"
    lines = [
        " FloodA5  Replay",
        " $s",
        @sprintf(" Sim time   %10s", _fmt_time(t)),
        @sprintf(" Frame      %10s", "$k / $nf"),
        " File  $fl",
        " $s",
        " Variable",
        "  $vl",
        @sprintf(" Max value  %10.4g", max_val),
        @sprintf(" Wet cells  %10d",   n_wet),
        " $s",
        " Mass Balance",
        " Vol added  " * _fmt_val(vol_added),
        " Vol domain " * _fmt_val(vol_domain),
        rl * _fmt_val(vol_removed),
        " MB error   " * _fmt_val(mb_err),
        " $s",
        @sprintf(" Mesh cells %10d", n_cells),
        @sprintf(" Resolution %10d", res),
        " $s",
        @sprintf(" Inputs     %10s", inputs_known ? "rebuilt" : "n/a"),
        @sprintf(" Outflow    %10s", out),
    ]
    return join(lines, "\n")
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
    _update_quiver!(segs_obs, tips_obs, rots_obs, sizes_obs, lw_obs, lons, lats, vel_u, vel_v, ax, usr_scale=1.0, stride=1, vel_thresh=QUIVER_MIN_SPEED)

Filter cells to those with non-trivial speed and update the three quiver Observables:
- `segs_obs`  — interleaved [tail, tip, tail, tip …] Point2f pairs for linesegments!
- `tips_obs`  — tip Point2f positions for scatter! arrowheads
- `rots_obs`  — arrowhead rotation angle (radians) for scatter! rotations

Arrow length is scaled so the peak-speed cell draws an arrow ≈2% of the
domain longitude span.  Cells with speed below `QUIVER_MIN_SPEED` (1 mm/s)
are omitted to avoid rendering zero-length arrows over dry cells.

The arrowhead rotation follows Makie's scatter convention: 0 rad = pointing
upward (+y), so we use `atan(u, v)` (note argument order) to get the angle
from north, matching the arrow direction in lon/lat space.
"""
const QUIVER_MIN_SPEED = 1e-5   # m/s — below this, cell is skipped in quiver (~0.01 mm/s)

function _update_quiver!(segs_obs   :: Observable{Vector{Point2f}},
                          tips_obs   :: Observable{Vector{Point2f}},
                          rots_obs   :: Observable{Vector{Float32}},
                          sizes_obs  :: Observable{Vector{Float32}},
                          lw_obs     :: Observable{Float32},
                          lons       :: Vector{Float64},
                          lats       :: Vector{Float64},
                          vel_u      :: Vector{Float64},
                          vel_v      :: Vector{Float64},
                          ax         :: Axis,
                          usr_scale  :: Float64 = 1.0,
                          stride     :: Int      = 1,
                          vel_thresh :: Float64  = QUIVER_MIN_SPEED)
    n = length(lons)

    peak = maximum(sqrt(vel_u[i]^2 + vel_v[i]^2) for i in 1:n; init=0.0)
    if peak <= vel_thresh
        segs_obs[]  = Point2f[]
        tips_obs[]  = Point2f[]
        rots_obs[]  = Float32[]
        sizes_obs[] = Float32[]
        lw_obs[]    = 1.5f0
        return
    end

    # Convert shaft lengths from screen pixels → data (lon/lat) degrees.
    # We read the axis pixel area and data limits at the moment of the call,
    # so the conversion is always current (works after zoom/pan too).
    # deg_per_px_lon and deg_per_px_lat give the data-space size of one pixel.
    lims      = ax.finallimits[]          # Rect2{Float32} in data coords
    px_area   = ax.scene.viewport[]       # Rect2 in device-independent units
    deg_per_px_lon = width(lims)  / width(px_area)
    deg_per_px_lat = height(lims) / height(px_area)
    # Use the average so arrows are symmetric regardless of aspect ratio
    deg_per_px = (deg_per_px_lon + deg_per_px_lat) / 2.0

    # Shaft: peak-speed → SHAFT_PEAK_PX pixels; floor → SHAFT_MIN_PX pixels.
    # Both scale linearly with usr_scale (slider 0.1 – 5.0).
    SHAFT_PEAK_PX = 28.0 * usr_scale   # pixels for the fastest cell
    SHAFT_MIN_PX  =  6.0 * usr_scale   # minimum visible shaft for any wet cell

    shaft_at_peak = SHAFT_PEAK_PX * deg_per_px          # in data degrees
    shaft_min     = SHAFT_MIN_PX  * deg_per_px

    # Head: screen-pixel size, 4px (dot) → 16px, scaled by usr_scale
    HEAD_MIN = 3.0f0
    HEAD_MAX = 16.0f0

    segs  = Point2f[]
    tips  = Point2f[]
    rots  = Float32[]
    sizes = Float32[]
    peak_lw = 1.5f0   # will be updated to scaled value for peak-speed arrow

    for i in 1:stride:n
        u = vel_u[i]; v = vel_v[i]
        speed = sqrt(u^2 + v^2)
        speed < vel_thresh && continue

        # Convert velocity direction to screen-pixel space so the shaft and
        # arrowhead always point in exactly the same visual direction.
        # ux_px/uy_px is the unit vector in screen pixels.
        ux_px = u / deg_per_px_lon
        uy_px = v / deg_per_px_lat
        px_len = sqrt(ux_px^2 + uy_px^2)
        px_len < 1e-12 && continue
        ux_px /= px_len;  uy_px /= px_len   # normalise to unit screen vector

        # Shaft length in data coordinates, using the screen-space direction
        sc     = max(shaft_at_peak * (speed / peak), shaft_min)
        tail   = Point2f(lons[i], lats[i])
        tip    = Point2f(lons[i] + Float32(ux_px * sc * deg_per_px_lon),
                         lats[i] + Float32(uy_px * sc * deg_per_px_lat))

        # Head size scales with speed and usr_scale
        frac      = Float32(speed / peak)
        head_size = HEAD_MIN + frac * (HEAD_MAX - HEAD_MIN)
        head_size = clamp(head_size * Float32(usr_scale), HEAD_MIN,
                          HEAD_MAX * Float32(usr_scale))

        # Shaft linewidth scales with head size so thick-head arrows have
        # proportionally thicker shafts (min 0.8px, max ~3px)
        lw = clamp(head_size * 0.12f0, 0.8f0, 3.0f0)

        # Rotation: atan(screen_x, screen_y) because Makie's :utriangle
        # points up (+y screen) at rotation=0, so we rotate from +y.
        # Note ux_px is the screen x-component, uy_px is screen y-component.
        rot = Float32(atan(ux_px, uy_px))

        push!(segs, tail);  push!(segs, tip)
        push!(tips, tip)
        push!(rots, rot)
        push!(sizes, head_size)
        peak_lw = max(peak_lw, lw)   # track max linewidth for this frame
    end

    segs_obs[]  = segs
    tips_obs[]  = tips
    rots_obs[]  = rots
    sizes_obs[] = sizes
    lw_obs[]    = peak_lw
end

"""
    _set_fields!(vis, depth, sat, vol, vel, vel_u, vel_v) → (idx, shown)

Set the per-cell data Observables (mesh order) and the displayed variable
with its colour range. Shared by `push_frame!` and replay.
"""
function _set_fields!(vis::MakieVis, depth, sat, vol, vel, vel_u, vel_v)
    vis.data_depth[]  = depth
    vis.data_sat[]    = sat
    vis.data_vol[]    = vol
    vis.data_vel[]    = vel
    # Update vel_v silently first (no notification), then notify vel_u once.
    # This ensures _update_quiver! sees the current vel_v when it fires.
    vis.data_vel_v.val = vel_v
    vis.data_vel_u[]   = vel_u   # triggers on(data_vel_u) with both current

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
    return idx, new_data
end

"Per-ring volume bar chart (ring mode only) from the current `data_vol`."
function _update_ring_bars!(vis::MakieVis)
    vis.n_rings > 0 || return
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
    return
end

function _pick_data(idx, d_depth, d_sat, d_vol, d_vel)
    idx == 1 && return d_depth[]
    idx == 2 && return d_sat[]
    idx == 3 && return d_vol[]
    return d_vel[]
end

function _auto_range(data::Vector{Float64}, var_idx::Int)
    var_idx == 2 && return (0.0, 1.0)   # saturation: fixed 0–1

    # For velocity (var_idx == 4) use the 99th percentile of positive values
    # to prevent a handful of near-dry cells with stale flux from collapsing
    # the colorscale. For all other variables use the true maximum.
    hi = if var_idx == 4
        pos = filter(v -> isfinite(v) && v > 1e-4, data)
        if length(pos) >= 10
            sort!(pos)
            pos[clamp(round(Int, 0.99 * length(pos)), 1, length(pos))]
        else
            maximum(data; init=0.0)
        end
    else
        maximum(data; init=0.0)
    end
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
        @sprintf(" Threads    %10d",  Threads.nthreads()),
        @sprintf(" Updated    %s",    Dates.format(now(), "HH:MM:SS")),
    ]
    return join(lines, "\n")
end

end # module MakieVisualiser
