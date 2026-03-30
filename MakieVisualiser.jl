# MakieVisualiser.jl
# ------------------
# Native desktop flood visualiser using GLMakie.
#
# Layout (1200 × 760 window)
# --------------------------
#
#   ┌────────────────────────────────────┬───────────────────┐
#   │                                    │  DIAGNOSTICS      │
#   │   Map panel                        │                   │
#   │   (A5 pentagons, variable colour)  │  Sim time / frame │
#   │                                    │  Active variable  │
#   │                                    │  Max / wet stats  │
#   │                                    │  Mesh info        │
#   │                                    │  Lon/lat bounds   │
#   │                                    │  Threads / clock  │
#   ├────────────────────────────────────┤                   │
#   │  [Depth ▾] colorbar label──────    │                   │
#   └────────────────────────────────────┴───────────────────┘
#
#  The dropdown (Menu widget) in row 3 selects which variable is mapped
#  to the polygon colour and colorbar.  Options:
#    • Depth (m)         — water_depth
#    • Saturation (0–1)  — fractional wetted area (SGS only; 1.0 for standard)
#    • Volume (m³)       — stored water volume per cell
#    • Velocity (m/s)    — scalar velocity magnitude
#
# Usage
# -----
#   include("MakieVisualiser.jl")
#   using .MakieVisualiser
#
#   vis = MakieVisualiser.start(mesh)
#   MakieVisualiser.push_frame!(vis, cell_ids, depths, saturations, volumes,
#                                              velocities, t)
#   MakieVisualiser.stop(vis)
#
# Install dependency:  using Pkg; Pkg.add("GLMakie")

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
const COL_TEXT  = RGBf(0.75, 0.88, 0.92)
const COL_LABEL = RGBf(0.60, 0.75, 0.82)
const COL_GRID  = RGBAf(1.0, 1.0, 1.0, 0.06)
const COL_EDGE  = RGBAf(0.15, 0.25, 0.35, 0.55)

# Display variable descriptors -----------------------------------------------
# Each entry: (menu label, field symbol, colormap, unit string)
const DISPLAY_VARS = [
    ("Depth (m)",        :depth,      :turbo,   "Water depth (m)"),
    ("Saturation (0–1)", :saturation, :Blues,   "Wetted area fraction"),
    ("Volume (m³)",      :volume,     :turbo,   "Stored volume (m³)"),
    ("Velocity (m/s)",   :velocity,   :plasma,  "Velocity (m/s)"),
]
const VAR_LABELS   = [v[1] for v in DISPLAY_VARS]
const VAR_SYMBOLS  = [v[2] for v in DISPLAY_VARS]
const VAR_CMAPS    = [v[3] for v in DISPLAY_VARS]
const VAR_UNITS    = [v[4] for v in DISPLAY_VARS]

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
State container for the Makie visualiser.

All mutable display state is held in `Observable`s so that GLMakie
re-renders automatically whenever `push_frame!` writes new values.
The `figure` field keeps the window alive for the lifetime of this struct.
"""
mutable struct MakieVis
    # Fixed mesh geometry — set once at construction, never mutated
    cell_ids   :: Vector{String}
    n_cells    :: Int
    resolution :: Int
    lon_range  :: Tuple{Float64, Float64}
    lat_range  :: Tuple{Float64, Float64}

    # Per-variable data Observables — one full vector per variable,
    # always aligned to mesh cell order.
    data_depth :: Observable{Vector{Float64}}
    data_sat   :: Observable{Vector{Float64}}
    data_vol   :: Observable{Vector{Float64}}
    data_vel   :: Observable{Vector{Float64}}

    # Display Observables — drive the polygon colours and colorbar
    display_data  :: Observable{Vector{Float64}}   # whichever var is selected
    colorrange    :: Observable{Tuple{Float64,Float64}}
    colormap_obs  :: Observable{Symbol}
    cbar_label    :: Observable{String}

    # Which variable the Menu is showing (index into DISPLAY_VARS)
    selected_var  :: Observable{Int}

    # Simulation metadata Observables
    sim_time   :: Observable{Float64}
    frame_idx  :: Observable{Int}
    diag_text  :: Observable{String}

    # GLMakie figure — must stay in scope to keep the window open
    figure     :: Figure

    # Shutdown flag — set false by stop(), polled by run_flood_model
    running    :: Threads.Atomic{Bool}
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    start(mesh; title="FloodA5 — Makie Viewer") → MakieVis

Open a GLMakie desktop window pre-loaded with the mesh geometry.
Returns immediately — the simulation loop can begin pushing frames at once.

The colorbar range is auto-scaled to the current frame's data range each
time `push_frame!` is called, so no fixed ceiling needs to be configured.
"""
function start(mesh; title::String = "FloodA5 — Makie Viewer")

    n        = length(mesh.cells)
    cell_ids = [c.id for c in mesh.cells]
    lons     = [c.center_lon for c in mesh.cells]
    lats     = [c.center_lat for c in mesh.cells]

    # Build closed polygon rings as Vector{Point2f} — one ring per cell
    polys = Vector{Vector{Point2f}}(undef, n)
    for (i, c) in enumerate(mesh.cells)
        ring = [Point2f(Float32(v[1]), Float32(v[2])) for v in c.boundary]
        first(ring) ≈ last(ring) || push!(ring, first(ring))
        polys[i] = ring
    end

    lon_range = extrema(lons)
    lat_range = extrema(lats)

    zeros_n = zeros(Float64, n)

    # --- Per-variable data Observables --------------------------------------
    data_depth = Observable(copy(zeros_n))
    data_sat   = Observable(copy(zeros_n))
    data_vol   = Observable(copy(zeros_n))
    data_vel   = Observable(copy(zeros_n))

    # --- Display Observables ------------------------------------------------
    selected_var = Observable(1)   # starts on Depth
    display_data = Observable(copy(zeros_n))
    colorrange   = Observable((0.0, 1.0))
    colormap_obs = Observable(VAR_CMAPS[1])
    cbar_label   = Observable(VAR_UNITS[1])

    # --- Simulation metadata ------------------------------------------------
    sim_time  = Observable(0.0)
    frame_idx = Observable(0)
    diag_text = Observable(
        _build_diag(0.0, 0, VAR_LABELS[1], 0.0, 0.0,
                    n, mesh.resolution, lon_range, lat_range)
    )

    # --- Figure & layout ----------------------------------------------------
    fig = Figure(size = (1200, 760), backgroundcolor = BG_DARK)

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

    pp = poly!(ax_map, polys;
        color       = display_data,
        colormap    = colormap_obs,
        colorrange  = colorrange,
        strokewidth = 0.4,
        strokecolor = COL_EDGE,
    )

    # Row 2: colorbar + dropdown in an HBox
    # Colorbar takes most of the width; the Menu sits to its right.
    cb = Colorbar(fig[2, 1], pp;
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
        color    = COL_LABEL,
        fontsize = 13,
        halign   = :right,
        tellwidth = false,
    )
    menu = Menu(menu_box[1, 2];
        options     = VAR_LABELS,
        default     = VAR_LABELS[1],
        width       = 180,
        tellwidth   = false,
    )

    # --- Sidebar (column 2, rows 1–3) ---------------------------------------
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

    text!(ax_side, 0.06, 0.97;
        text     = diag_text,
        align    = (:left, :top),
        space    = :relative,
        font     = "Courier New",
        fontsize = 12,
        color    = COL_TEXT,
    )

    # --- Column / row sizing ------------------------------------------------
    colsize!(fig.layout, 1, Relative(0.74))
    colsize!(fig.layout, 2, Relative(0.26))
    rowsize!(fig.layout, 2, Fixed(42))   # colorbar
    rowsize!(fig.layout, 3, Fixed(36))   # dropdown

    # --- Menu callback: switch display variable ------------------------------
    on(menu.selection) do label
        idx = findfirst(==(label), VAR_LABELS)
        idx === nothing && return
        selected_var[] = idx
        # Swap the data being displayed
        new_data = _pick_data(idx, data_depth, data_sat, data_vol, data_vel)
        display_data[] = new_data
        colorrange[]   = _auto_range(new_data, idx)
        colormap_obs[] = VAR_CMAPS[idx]
        cbar_label[]   = VAR_UNITS[idx]
        # Update sidebar to reflect new variable
        diag_text[] = _build_diag(
            sim_time[], frame_idx[], label,
            maximum(new_data; init=0.0),
            count(>(1e-4), new_data),
            n, mesh.resolution, lon_range, lat_range,
        )
    end

    # --- Open window --------------------------------------------------------
    display(fig)
    @info "MakieVisualiser: window open — $(n) cells at resolution $(mesh.resolution)"

    return MakieVis(
        cell_ids, n, mesh.resolution, lon_range, lat_range,
        data_depth, data_sat, data_vol, data_vel,
        display_data, colorrange, colormap_obs, cbar_label,
        selected_var,
        sim_time, frame_idx, diag_text,
        fig,
        Threads.Atomic{Bool}(true),
    )
end

"""
    push_frame!(vis, cell_ids, depths, saturations, volumes, velocities, t)

Update the viewer with a new simulation timestep.

All five data arrays are parallel to `cell_ids` and cover whatever cells
have non-zero state this step.  Each is re-aligned to the full mesh cell
order (absent cells default to 0.0).  Only the currently selected variable
is re-rendered; the others are stored and become visible when the user
switches the dropdown.
"""
function push_frame!(vis         :: MakieVis,
                     cell_ids    :: Vector{String},
                     depths      :: Vector{Float64},
                     saturations :: Vector{Float64},
                     volumes     :: Vector{Float64},
                     velocities  :: Vector{Float64},
                     t           :: Float64)
    vis.running[] || return

    # Re-align all four arrays to the fixed mesh cell order
    function _align(vals)
        d = Dict{String,Float64}(zip(cell_ids, vals))
        [get(d, id, 0.0) for id in vis.cell_ids]
    end

    new_depth = _align(depths)
    new_sat   = _align(saturations)
    new_vol   = _align(volumes)
    new_vel   = _align(velocities)

    # Write all four data stores
    vis.data_depth[] = new_depth
    vis.data_sat[]   = new_sat
    vis.data_vol[]   = new_vol
    vis.data_vel[]   = new_vel

    # Update the displayed variable
    idx      = vis.selected_var[]
    new_data = _pick_data(idx, vis.data_depth, vis.data_sat,
                          vis.data_vol, vis.data_vel)
    vis.display_data[] = new_data
    vis.colorrange[]   = _auto_range(new_data, idx)

    # Simulation metadata
    frame      = vis.frame_idx[] + 1
    vis.sim_time[]  = t
    vis.frame_idx[] = frame

    vis.diag_text[] = _build_diag(
        t, frame, VAR_LABELS[idx],
        maximum(new_data; init=0.0),
        count(>(1e-4), new_data),
        vis.n_cells, vis.resolution,
        vis.lon_range, vis.lat_range,
    )
end

"""
    stop(vis)

Close the GLMakie window and mark the visualiser as stopped.
Safe to call from any thread.
"""
function stop(vis::MakieVis)
    vis.running[] = false
    try
        GLMakie.closeall()
    catch
    end
    @info "MakieVisualiser: window closed"
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""Return the data Observable for the given variable index."""
function _pick_data(idx, d_depth, d_sat, d_vol, d_vel)
    idx == 1 && return d_depth[]
    idx == 2 && return d_sat[]
    idx == 3 && return d_vol[]
    return d_vel[]
end

"""
Auto-scale colorrange for the current data.
Saturation is always [0, 1].  Others use [0, max] with a small floor
so the colorbar doesn't collapse to a single colour on a dry mesh.
"""
function _auto_range(data::Vector{Float64}, var_idx::Int)
    var_idx == 2 && return (0.0, 1.0)   # saturation: fixed [0,1]
    hi = maximum(data; init=0.0)
    hi = hi > 1e-6 ? hi : 1.0           # floor: avoid zero-width range
    return (0.0, hi)
end

"""Format elapsed simulation seconds as a human-readable string."""
function _fmt_time(t::Float64)::String
    t < 60   && return @sprintf("%d s",     round(Int, t))
    t < 3600 && return @sprintf("%.1f min", t / 60)
    return              @sprintf("%.2f h",  t / 3600)
end

"""
Build the pre-formatted monospace diagnostics string for the sidebar.
`max_val` and `active_n` refer to whichever variable is currently selected.
"""
function _build_diag(
    t         :: Float64,
    frame     :: Int,
    var_label :: String,
    max_val   :: Float64,
    active_n,           # count of cells above threshold — Int or Float
    n_cells   :: Int,
    res       :: Int,
    lon_range :: Tuple{Float64, Float64},
    lat_range :: Tuple{Float64, Float64},
)::String
    s = "─────────────────────"
    # Trim var_label to fit the fixed-width column
    vl = length(var_label) > 13 ? var_label[1:13] : rpad(var_label, 13)
    lines = [
        " FloodA5  Diagnostics",
        " $s",
        @sprintf(" Sim time   %10s",  _fmt_time(t)),
        @sprintf(" Frame      %10d",  frame),
        " $s",
        " Variable",
        "  $vl",
        @sprintf(" Max value  %10.4g", max_val),
        @sprintf(" Active     %10d",  Int(active_n)),
        " $s",
        @sprintf(" Mesh cells %10d",  n_cells),
        @sprintf(" Resolution %10d",  res),
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
