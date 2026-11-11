"""
FloodViewer.jl
--------------
Standalone post-processing viewer for FloodA5 simulation output.

Loads mesh geometry (from the output file itself for format-2 files,
FloodA5 0.4.0+; from the GeoParquet mesh if one is given — needed for
format-1 files) and time-series data from a sim.h5 HDF5 file, then
displays an interactive GLMakie window with:
  - Pentagon map coloured by the selected variable
  - Frame slider + Play/Pause button
  - Variable selector (water_depth, saturation, volume)
  - Diagnostics panel (time, wet cells, max value, mass conservation)
  - Export frame as PNG

Usage
-----
  julia --project=. --threads auto visualisation/FloodViewer.jl sim.h5            # format 2
  julia --project=. --threads auto visualisation/FloodViewer.jl mesh.parquet sim.h5
  julia --project=. --threads auto visualisation/FloodViewer.jl sim.h5 --var saturation
  julia --project=. --threads auto visualisation/FloodViewer.jl sim.h5 --fps 10

The main viewer for output files is the live viewer's replay mode
(visualisation/replay.jl); this one is a simpler alternative.

Dependencies
------------
  The root project (GLMakie, HDF5); the mesh is read in pure Julia
  (analysis/FloodA5Analysis.jl), no Python.
"""

using GLMakie
using HDF5
isdefined(Main, :FloodA5Output) || include(joinpath(dirname(@__DIR__), "io", "FloodA5Output.jl"))
const FO = FloodA5Output
isdefined(Main, :FloodA5Analysis) ||
    include(joinpath(dirname(@__DIR__), "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
using Printf
using Dates
using Statistics: quantile

# ---------------------------------------------------------------------------
# Load mesh geometry
# ---------------------------------------------------------------------------

"""
    load_mesh_geometry(path) → (polys, lons, lats, elevations, n_cells)

Load the mesh geometry — from a GeoParquet mesh (`FA.load_mesh_parquet`,
row order) or from a format-2 output file's `/mesh` (`FA.mesh_geometry`) —
and return:
  polys      — Vector of polygon vertex arrays, each (2 × n_verts), for poly!
  lons/lats  — cell centre coordinates
  elevations — bed elevation per cell (NaN if not sampled)
  n_cells    — number of cells
"""
function load_mesh_geometry(path::String)
    @info "Loading mesh geometry from $(basename(path))..."
    mesh = if lowercase(splitext(path)[2]) in (".h5", ".hdf5")
        FO.format_version(path) >= 2 ||
            error("$path is a format-1 output file (FloodA5 ≤ 0.3.x) without mesh " *
                  "geometry: give the mesh's .parquet as well")
        FA.mesh_geometry(path)
    else
        FA.load_mesh_parquet(path)
    end
    n = length(mesh)
    isempty(mesh.polygons) && error("no cell polygons in $path")

    # GLMakie poly! takes (2 × n_verts) matrices (x=lon, y=lat columns)
    polys = mesh.polygons
    lons  = mesh.lon
    lats  = mesh.lat
    elevs = mesh.elevation

    @info "  $(n) cells loaded"
    return polys, lons, lats, elevs, n
end

# ---------------------------------------------------------------------------
# Load HDF5 frames
# ---------------------------------------------------------------------------

"""
    SimData

All simulation frames loaded into memory.

Fields
------
  t           — simulation times (s), length n_frames
  water_depth — (n_cells × n_frames)
  volume      — (n_cells × n_frames)
  saturation  — (n_cells × n_frames)
  velocity    — (n_cells × n_frames)
  n_frames    — number of frames
  n_cells     — number of cells
"""
struct SimData
    t           :: Vector{Float64}
    water_depth :: Matrix{Float64}
    volume      :: Matrix{Float64}
    saturation  :: Matrix{Float64}
    velocity    :: Matrix{Float64}
    n_frames    :: Int
    n_cells     :: Int
end

function load_sim_data(h5_path::String)::SimData
    @info "Loading simulation data from $(basename(h5_path))..."
    isfile(h5_path) || error("HDF5 file not found: $h5_path")

    HDF5.h5open(h5_path, "r") do fid
        # Output format 1 or 2 (io/FloodA5Output.jl)
        n_frames = FO.n_frames(fid)
        n_frames > 0 || error("No frames found in $h5_path")

        t           = FO.frame_times(fid)
        water_depth = FO.read_var(fid, "water_depth")
        n_cells     = size(water_depth, 1)
        volume      = FO.read_var(fid, "volume")
        saturation  = FO.has_var(fid, "saturation") ?
                      FO.read_var(fid, "saturation") : Float64.(water_depth .> 1e-4)
        velocity    = FO.read_var(fid, "velocity")

        @info "  $n_frames frames × $n_cells cells  " *
              "(t = $(round(t[1],digits=0))s … $(round(t[end],digits=0))s)"

        return SimData(t, water_depth, volume, saturation, velocity,
                       n_frames, n_cells)
    end
end

# ---------------------------------------------------------------------------
# Variable helpers
# ---------------------------------------------------------------------------

const VAR_LABELS = Dict(
    "water_depth" => "Water Depth (m)",
    "saturation"  => "Wetted Fraction (0–1)",
    "volume"      => "Stored Volume (m³)",
    "velocity"    => "Velocity (m/s)",
)

function get_var_data(sim::SimData, var::String, frame::Int)::Vector{Float64}
    frame = clamp(frame, 1, sim.n_frames)
    if     var == "water_depth";  return sim.water_depth[:, frame]
    elseif var == "saturation";   return sim.saturation[:, frame]
    elseif var == "volume";       return sim.volume[:, frame]
    elseif var == "velocity";     return sim.velocity[:, frame]
    else;  error("Unknown variable: $var")
    end
end

function var_colormap(var::String)
    if     var == "water_depth";  return :blues
    elseif var == "saturation";   return :YlGnBu
    elseif var == "volume";       return :blues
    elseif var == "velocity";     return :plasma
    else;  return :viridis
    end
end

# ---------------------------------------------------------------------------
# Main viewer
# ---------------------------------------------------------------------------

"""
    launch_viewer(mesh_path, h5_path; initial_var="water_depth", fps=5)

Open the interactive GLMakie flood viewer. `mesh_path`: the GeoParquet mesh,
or the output file itself (format 2: it carries the mesh).
"""
function launch_viewer(mesh_path::String, h5_path::String;
                       initial_var::String = "water_depth",
                       fps::Float64        = 5.0)

    polys, lons, lats, elevs, n_cells = load_mesh_geometry(mesh_path)
    sim = load_sim_data(h5_path)

    sim.n_cells == n_cells ||
        @warn "Cell count mismatch: mesh=$(n_cells), HDF5=$(sim.n_cells). " *
              "Are you loading the right parquet file?"
    n = min(n_cells, sim.n_cells)

    # ── Observable state ────────────────────────────────────────────────────
    frame_obs   = Observable(1)
    var_obs     = Observable(initial_var)
    playing_obs = Observable(false)
    clim_obs    = Observable((0.0, 1.0))

    # Current data vector (colour values)
    data_obs = @lift begin
        d = get_var_data(sim, $var_obs, $frame_obs)
        replace(d, NaN => 0.0)
    end

    # Colour limits: 0 to 99th percentile of max-over-time per cell
    # (recomputed when variable changes, stable during playback)
    function compute_clim(var::String)
        all_data = get_var_data.(Ref(sim), Ref(var), 1:sim.n_frames)
        flat = filter(isfinite, vcat(all_data...))
        isempty(flat) && return (0.0, 1.0)
        hi = quantile(flat, 0.99)
        hi <= 0.0 && return (0.0, 1.0)
        return (0.0, hi)
    end

    # ── Window layout ────────────────────────────────────────────────────────
    fig = Figure(size=(1300, 750), backgroundcolor=:white)

    # Map panel
    ax_map = Axis(fig[1, 1],
        xlabel = "Longitude",
        ylabel = "Latitude",
        title  = @lift(@sprintf("Frame %d / %d   t = %.0f s",
                                $frame_obs, sim.n_frames, sim.t[$frame_obs])),
        aspect = DataAspect(),
    )

    # Draw polygons — split into dry (grey) and wet (coloured)
    # We use one poly! call per cell to allow per-cell colouring.
    # For performance, we draw all cells as a mesh with per-face colours
    # using the poly! approach with a colour observable.
    cell_colors = @lift begin
        vals = $data_obs
        lo, hi = $clim_obs
        span  = hi - lo
        cmap  = to_colormap(var_colormap($var_obs))
        n_cm  = length(cmap)
        cols  = Vector{RGBAf}(undef, n)
        for i in 1:n
            v = isfinite(vals[i]) ? vals[i] : 0.0
            if v < 1e-4 || span <= 0.0
                cols[i] = RGBAf(0.88, 0.88, 0.88, 1.0)
            else
                t_norm  = clamp((v - lo) / span, 0.0, 1.0)
                idx     = clamp(round(Int, t_norm * (n_cm - 1)) + 1, 1, n_cm)
                cols[i] = RGBAf(cmap[idx])
            end
        end
        cols
    end

    # Draw all pentagons in a single poly! call with a vector colour observable.
    # GLMakie accepts a Vector{<:AbstractVector{Point2f}} for the polygons
    # and a matching color vector — much faster than one poly! per cell.
    all_pts = [[Point2f(polys[i][1, j], polys[i][2, j])
                for j in axes(polys[i], 2)]
               for i in 1:n]
    poly!(ax_map, all_pts;
          color       = cell_colors,
          strokecolor = (:black, 0.20),
          strokewidth = 0.4)

    # Colorbar
    cb_cmap = @lift to_colormap(var_colormap($var_obs))
    cb_lims = lift(identity, clim_obs)
    Colorbar(fig[2, 1];
        colormap  = cb_cmap,
        limits    = cb_lims,
        label     = @lift(get(VAR_LABELS, $var_obs, $var_obs)),
        vertical  = false,
        width     = Relative(0.6),
        tellwidth = false,
    )

    # ── Controls panel (right column) ────────────────────────────────────────
    panel = fig[1:2, 2]
    panel_grid = GridLayout(panel)

    # Title
    Label(panel_grid[1, 1],
        "FloodA5 Viewer";
        fontsize=16, tellwidth=false)

    # Variable selector
    Label(panel_grid[2, 1], "Variable:"; halign=:left, tellwidth=false)
    var_menu = Menu(panel_grid[3, 1];
        options  = ["water_depth", "saturation", "volume", "velocity"],
        default  = initial_var,
        tellwidth = true)

    on(var_menu.selection) do sel
        var_obs[] = sel
        clim_obs[] = compute_clim(sel)
    end

    # Frame slider
    Label(panel_grid[4, 1], "Frame:"; halign=:left, tellwidth=false)
    slider = Slider(panel_grid[5, 1];
        range    = 1:sim.n_frames,
        startvalue = 1,
        tellwidth = true)
    connect!(frame_obs, slider.value)

    # Play / Pause
    play_btn = Button(panel_grid[6, 1];
        label     = @lift($playing_obs ? "Pause [||]" : "Play  [>]"),
        tellwidth = true)

    on(play_btn.clicks) do _
        playing_obs[] = !playing_obs[]
    end

    # FPS display
    Label(panel_grid[7, 1],
        "Speed: $(fps) fps  ($(round(1/fps, digits=2))s/frame)";
        halign=:left, tellwidth=false)

    # Diagnostics
    Label(panel_grid[8, 1], "─── Diagnostics ───";
        tellwidth=false)

    diag_text = @lift begin
        f    = $frame_obs
        var  = $var_obs
        vals = get_var_data(sim, var, f)
        t_s  = sim.t[f]
        n_wet = count(>(1e-4), sim.water_depth[:, f])
        max_d = maximum(filter(isfinite, sim.water_depth[:, f]); init=0.0)
        max_v = maximum(filter(isfinite, vals); init=0.0)
        tot_v = sum(filter(isfinite, sim.volume[:, f]))

        @sprintf("""
Time         : %.0f s (%.2f hr)
Frame        : %d / %d
Wet cells    : %d / %d  (%.0f%%)
Max depth    : %.4f m
Max %s : %.4f
Total volume : %.0f m³
        """,
        t_s, t_s/3600,
        f, sim.n_frames,
        n_wet, sim.n_cells, 100*n_wet/sim.n_cells,
        max_d,
        rpad(var, 8), max_v,
        tot_v)
    end

    Label(panel_grid[9, 1], diag_text;
        halign=:left, tellwidth=false, fontsize=11)

    # Export button
    export_btn = Button(panel_grid[10, 1];  # (no emoji: Makie 0.24 has no font for them)
        label="Export PNG", tellwidth=true)

    on(export_btn.clicks) do _
        fname = @sprintf("flooda5_frame_%04d_%s.png",
                         frame_obs[], var_obs[])
        save(fname, fig)
        @info "Saved: $fname"
    end

    # Mesh info
    Label(panel_grid[11, 1],
        @sprintf("Mesh: %d cells  |  %d frames", n, sim.n_frames);
        halign=:left, tellwidth=false, fontsize=11, color=:grey50)

    lon_range = extrema(lons)
    lat_range = extrema(lats)
    Label(panel_grid[12, 1],
        @sprintf("Lon %.4f–%.4f  Lat %.4f–%.4f",
                 lon_range..., lat_range...);
        halign=:left, tellwidth=false, fontsize=10, color=:grey50)

    colsize!(fig.layout, 1, Relative(0.78))
    colsize!(fig.layout, 2, Relative(0.22))
    rowsize!(fig.layout, 2, Auto())

    # Initialise colour limits
    clim_obs[] = compute_clim(initial_var)

    # ── Playback loop ─────────────────────────────────────────────────────────
    @async begin
        while true
            sleep(1.0 / fps)
            isopen(fig.scene) || break
            playing_obs[] || continue
            next = frame_obs[] + 1
            if next > sim.n_frames
                playing_obs[] = false
            else
                set_close_to!(slider, next)
            end
        end
    end

    display(fig)
    @info "Viewer ready — close the window to exit."
    return fig
end

# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

function main()
    args = copy(ARGS)

    # positional: [mesh.parquet] sim.h5 (the values of --var/--fps excluded)
    pos = [a for (i, a) in enumerate(args)
           if !startswith(a, "--") && !(i > 1 && args[i-1] in ("--var", "--fps"))]
    if "--help" in args || "-h" in args || !(1 <= length(pos) <= 2)
        println("""
FloodViewer.jl — interactive post-processing viewer for FloodA5 output

Usage:
    julia --project=. --threads auto visualisation/FloodViewer.jl [mesh.parquet] <sim.h5> [options]

The mesh is read from sim.h5 for format-2 files (FloodA5 0.4.0+); format-1
files (≤ 0.3.x) need the GeoParquet mesh as well.

Options:
    --var VAR     Initial variable to display (default: water_depth)
                  Choices: water_depth, saturation, volume, velocity
    --fps N       Playback speed in frames per second (default: 5)
    --help        Show this help

Examples:
    julia --project=. --threads auto visualisation/FloodViewer.jl sim.h5
    julia --project=. --threads auto visualisation/FloodViewer.jl meshdem2.parquet sim.h5
    julia --project=. --threads auto visualisation/FloodViewer.jl sim.h5 --var saturation --fps 10
""")
        exit(0)
    end

    h5_path   = pos[end]
    mesh_path = length(pos) == 2 ? pos[1] : h5_path

    isfile(mesh_path) || error("Mesh file not found: $mesh_path")
    isfile(h5_path)   || error("HDF5 file not found: $h5_path")

    # --var
    var_idx = findfirst(==("--var"), args)
    initial_var = var_idx !== nothing ? args[var_idx + 1] : "water_depth"
    initial_var in keys(VAR_LABELS) ||
        error("Unknown variable '$initial_var'. Choose from: $(join(keys(VAR_LABELS), ", "))")

    # --fps
    fps_idx = findfirst(==("--fps"), args)
    fps = fps_idx !== nothing ? parse(Float64, args[fps_idx + 1]) : 5.0

    fig = launch_viewer(mesh_path, h5_path;
                        initial_var = initial_var,
                        fps         = fps)

    # Keep alive until window is closed
    try
        while isopen(fig.scene)
            sleep(0.5)
        end
    catch e
        e isa InterruptException || rethrow(e)
    end
end

# Run from the command line; `include` only defines the functions
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
