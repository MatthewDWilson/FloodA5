#!/usr/bin/env julia
# =============================================================================
# view_h5_output.jl
#
# Interactive viewer and animation tool for FloodA5 HDF5 simulation output.
#
# Usage:
#
#   Interactive map at final frame:
#     julia --project=. visualisation/view_h5_output.jl sim.h5
#
#   Interactive map + cell time series:
#     julia --project=. visualisation/view_h5_output.jl sim.h5 \
#         --cell 6345ddac19800000
#
#   Interactive map at a specified time:
#     julia --project=. visualisation/view_h5_output.jl sim.h5 \
#         --time 7200 --var depth
#
#   Interactive map at a specified figure resolution:
#     julia --project=. visualisation/view_h5_output.jl sim.h5 \
#         --resolution 1600x900
#
#   Export animation frames:
#     julia --project=. visualisation/view_h5_output.jl sim.h5 \
#         --animate \
#         --mesh test/carlisle/carlisle_mesh18_sgs.parquet \
#         --frame-step 100 \
#         --var depth \
#         --output-dir testanimation
#
#   Export animation using local East-North metre coordinates:
#     julia --project=. visualisation/view_h5_output.jl sim.h5 \
#         --animate \
#         --mesh test/carlisle/carlisle_mesh18_sgs.parquet \
#         --frame-step 100 \
#         --var depth \
#         --output-dir testanimation \
#         --local-units
#
# Animation mode renders the actual A5 pentagonal cell boundaries from the
# supplied GeoParquet mesh. The HDF5 file deliberately does not contain mesh
# geometry, so --mesh is required for --animate.
#
# By default animation geometry is plotted in geographic longitude/latitude
# coordinates. Use --local-units to transform the geometry to a local
# East-North coordinate system in metres centred on the simulation domain.
#
# Frames are written as:
#     frame_000001.png
#     frame_000002.png
#     ...
#
# They can subsequently be assembled into an MP4 with ffmpeg.
#
# CLI options:
#
#   --cell ID
#       Cell to highlight and plot as a time series.
#
#   --var NAME
#       Variable to display:
#         depth, volume, saturation, velocity
#
#   --time SECONDS
#       Requested time for the interactive map.
#
#   --animate
#       Export a sequence of PNG animation frames instead of opening the
#       interactive viewer.
#
#   --mesh PATH
#       GeoParquet mesh containing the actual A5 cell boundaries.
#       Required with --animate.
#
#   --frame-step N
#       Output every Nth simulation frame. The final frame is always included.
#
#   --output-dir PATH
#       Directory for animation frames. Default: animation
#
#   --resolution WIDTHxHEIGHT
#       Figure/output resolution in pixels. Default: 1920x1080.
#       Examples: 1920x1080, 1600x900, 1200x800.
#
#   --local-units
#       For animation, transform mesh longitude/latitude coordinates to a
#       local East-North coordinate system in metres. Without this flag,
#       animation uses geographic longitude/latitude coordinates directly.
#
# =============================================================================

using HDF5
using GLMakie
import CairoMakie
using Printf
using Statistics
using DataFrames

# A5Grid provides load_mesh_geoparquet(), _to_hex(), and cell boundaries.
const _A5GRID_PATH = joinpath(@__DIR__, "..", "mesh", "A5Grid.jl")
include(_A5GRID_PATH)
using .A5Grid

const R_EARTH = 6_371_000.0


# =============================================================================
# CLI parsing
# =============================================================================

"""
Parse a resolution specification of the form WIDTHxHEIGHT.

Examples:
    1400x820
    1920x1080
    1600X900
    1600,900

Returns:
    (width, height)
"""
function parse_resolution(spec::String)
    s = strip(spec)

    # Accept x, X, or comma as separator.
    parts = split(s, r"[xX,]")

    length(parts) == 2 ||
        error(
            "--resolution must have the form WIDTHxHEIGHT, " *
            "for example 1400x820; got: $spec"
        )

    width = tryparse(Int, strip(parts[1]))
    height = tryparse(Int, strip(parts[2]))

    isnothing(width) &&
        error("Invalid resolution width in '$spec'")

    isnothing(height) &&
        error("Invalid resolution height in '$spec'")

    width > 0 ||
        error("Resolution width must be greater than zero")

    height > 0 ||
        error("Resolution height must be greater than zero")

    return (width, height)
end


function parse_args(args)
    length(args) < 1 && error("""
        Usage:
          julia --project=. view_h5_output.jl <h5file>
              [--cell ID] [--var NAME] [--time SECONDS]
              [--resolution WIDTHxHEIGHT]

        Animation:
          julia --project=. view_h5_output.jl <h5file>
              --animate --mesh MESH.parquet
              [--frame-step N] [--var NAME] [--output-dir DIR]
              [--resolution WIDTHxHEIGHT] [--local-units]
    """)

    h5file       = args[1]
    cell_id      = nothing
    varname      = "depth"
    time_req     = nothing

    animate      = false
    mesh_path    = nothing
    frame_step   = 1
    output_dir   = "animation"

    resolution   = (1920, 1080)
    local_units  = false

    i = 2
    while i <= length(args)
        arg = args[i]

        if arg == "--cell" && i < length(args)
            cell_id = lowercase(lstrip(args[i + 1], '0'))
            i += 2

        elseif arg == "--var" && i < length(args)
            varname = args[i + 1]
            i += 2

        elseif arg == "--time" && i < length(args)
            time_req = tryparse(Float64, args[i + 1])
            isnothing(time_req) &&
                error("--time must be a number (seconds), got: $(args[i + 1])")
            i += 2

        elseif arg == "--animate"
            animate = true
            i += 1

        elseif arg == "--mesh" && i < length(args)
            mesh_path = args[i + 1]
            i += 2

        elseif arg == "--frame-step" && i < length(args)
            frame_step = tryparse(Int, args[i + 1])
            isnothing(frame_step) &&
                error("--frame-step must be an integer, got: $(args[i + 1])")
            frame_step > 0 ||
                error("--frame-step must be greater than zero")
            i += 2

        elseif arg == "--output-dir" && i < length(args)
            output_dir = args[i + 1]
            i += 2

        elseif arg == "--resolution" && i < length(args)
            resolution = parse_resolution(args[i + 1])
            i += 2

        elseif arg == "--local-units"
            local_units = true
            i += 1

        else
            error("Unknown or incomplete argument: $arg")
        end
    end

    return (
        h5file      = h5file,
        cell_id     = cell_id,
        varname     = varname,
        time_req    = time_req,
        animate     = animate,
        mesh_path   = mesh_path,
        frame_step  = frame_step,
        output_dir  = output_dir,
        resolution  = resolution,
        local_units = local_units,
    )
end


# =============================================================================
# Frame selection
# =============================================================================

function select_frame(times::Vector{Float64}, time_req)
    isnothing(time_req) && return length(times), nothing

    t_end = times[end]

    if time_req > t_end
        @warn "Requested time $(time_req)s is beyond the simulation end " *
              "($(t_end)s). Showing final frame."
        return length(times), t_end
    end

    diffs = abs.(times .- time_req)
    fi = argmin(diffs)
    nearest = times[fi]

    if diffs[fi] > 0.5
        @warn "Requested time $(time_req)s not in output. Nearest frame: " *
              "$(round(nearest, digits=1))s " *
              "(Δ=$(round(diffs[fi], digits=1))s)."
    end

    return fi, nearest
end


# =============================================================================
# HDF5 reading
# =============================================================================

function load_h5(h5file::String)
    isfile(h5file) || error("File not found: $h5file")
    println("Reading: $h5file")

    h5open(h5file, "r") do f
        cell_ids = read(f["mesh/cell_ids"])
        lons     = Float64.(read(f["mesh/center_lons"]))
        lats     = Float64.(read(f["mesh/center_lats"]))
        elevs    = Float64.(read(f["mesh/elevations"]))

        frames   = sort(keys(f["frames"]))
        n_frames = length(frames)
        n_cells  = length(cell_ids)

        println("  Cells   : $n_cells")
        println("  Frames  : $n_frames")

        times    = Vector{Float64}(undef, n_frames)
        depth_ts = Matrix{Float64}(undef, n_frames, n_cells)
        vol_ts   = Matrix{Float64}(undef, n_frames, n_cells)
        sat_ts   = Matrix{Float64}(undef, n_frames, n_cells)
        vel_ts   = Matrix{Float64}(undef, n_frames, n_cells)

        for (fi, fr) in enumerate(frames)
            g = f["frames/$fr"]
            times[fi]       = Float64(read(g["t"]))
            depth_ts[fi, :] = read(g["water_depth"])
            vol_ts[fi, :]   = read(g["volume"])
            sat_ts[fi, :]   = read(g["saturation"])
            vel_ts[fi, :]   = read(g["velocity"])
        end

        t_hr = times ./ 3600.0

        final_depth = depth_ts[end, :]
        final_vol   = vol_ts[end, :]

        println("  Duration: $(round(times[end] / 3600, digits=2)) h")
        println("  Final domain volume: " *
                "$(round(sum(final_vol) / 1e6, digits=3)) × 10⁶ m³")
        println("  Final max depth    : " *
                "$(round(maximum(final_depth), digits=3)) m")

        wet_final = final_depth[final_depth .> 0.001]

        if !isempty(wet_final)
            println("  Final mean depth   : " *
                    "$(round(mean(wet_final), digits=3)) m (wet cells)")
        else
            println("  Final mean depth   : no wet cells")
        end

        return (
            cell_ids   = cell_ids,
            lons       = lons,
            lats       = lats,
            elevs      = elevs,
            times      = times,
            t_hr       = t_hr,
            depth      = depth_ts,
            volume     = vol_ts,
            saturation = sat_ts,
            velocity   = vel_ts,
            n_cells    = n_cells,
            n_frames   = n_frames,
        )
    end
end


# =============================================================================
# Cell ID helpers
# =============================================================================

"""
Normalise an HDF5 or mesh cell ID to FloodA5's 16-character lower-case
hexadecimal representation.
"""
function normalise_cell_id(id)
    if id isa Integer
        return A5Grid._to_hex(UInt64(id))
    else
        s = lowercase(strip(string(id)))
        s = startswith(s, "0x") ? s[3:end] : s
        return A5Grid._to_hex(parse(UInt64, s, base=16))
    end
end


function find_cell(cell_ids::Vector, target_id::String)
    t = normalise_cell_id(target_id)

    for (i, cid) in enumerate(cell_ids)
        normalise_cell_id(cid) == t && return i
    end

    return nothing
end


# =============================================================================
# Variable selector
# =============================================================================

function get_var(data, varname::String)
    varname == "depth" &&
        return data.depth, "Water depth (m)", "turbo"

    varname == "volume" &&
        return data.volume, "Volume (m³)", "turbo"

    varname == "saturation" &&
        return data.saturation, "Saturation (0–1)", "blues"

    varname == "velocity" &&
        return data.velocity, "Velocity (m/s)", "plasma"

    error(
        "Unknown variable '$varname'. " *
        "Choose: depth, volume, saturation, velocity"
    )
end


# =============================================================================
# Mesh loading
# =============================================================================

"""
Load the GeoParquet mesh through the same A5Grid pathway used by
visualise_mesh.jl.

Returns:
    df   DataFrame containing cell geometry and relevant mesh information
    idx  Dict{String,Int} mapping normalised cell ID -> row
"""
function load_parquet(path::String)
    isfile(path) || error("Mesh file not found: $path")

    mesh = A5Grid.load_mesh_geoparquet(path)

    n = length(mesh.cells)
    ids = [normalise_cell_id(c.id) for c in mesh.cells]

    idx = Dict{String,Int}(ids[i] => i for i in 1:n)

    df = DataFrame(
        cell_id    = ids,
        center_lon = [Float64(c.center_lon) for c in mesh.cells],
        center_lat = [Float64(c.center_lat) for c in mesh.cells],
        boundary   = [c.boundary for c in mesh.cells],
    )

    # Retain the same static information as visualise_mesh.jl where present.
    for (col, key) in [
        (:elevation, :elevation),
        (:sgs_cell_area, :sgs_cell_area),
        (:sgs_z_min, :sgs_z_min),
        (:sgs_z_max, :sgs_z_max),
    ]
        key_string = String(key)

        if haskey(mesh.static_vars, key_string)
            df[!, col] = mesh.static_vars[key_string]
        end
    end

    return df, idx
end


"""
Return the actual pentagonal boundary stored in the GeoParquet mesh.

This is deliberately the same geometry used by visualise_mesh.jl. No
approximation from cell centres is used for animation.
"""
function get_boundary(df, row::Int)
    bnd = df[row, :boundary]
    return [(Float64(v[1]), Float64(v[2])) for v in bnd]
end


"""
Project geographic boundary vertices to a local East-North coordinate system
centred on (cx, cy). This is the same local projection used by visualise_mesh.jl.
"""
function to_local_m(
    boundary::Vector,
    cx::Float64,
    cy::Float64,
)
    cos_cy = cosd(cy)

    xs = [
        (v[1] - cx) * deg2rad(1) * R_EARTH * cos_cy
        for v in boundary
    ]

    ys = [
        (v[2] - cy) * deg2rad(1) * R_EARTH
        for v in boundary
    ]

    return xs, ys
end


"""
Return either geographic longitude/latitude coordinates or local East/North
coordinates in metres, depending on the --local-units option.
"""
function project_boundary(
    boundary::Vector,
    cx::Float64,
    cy::Float64,
    local_units::Bool,
)
    if local_units
        return to_local_m(boundary, cx, cy)
    else
        xs = [Float64(v[1]) for v in boundary]
        ys = [Float64(v[2]) for v in boundary]
        return xs, ys
    end
end


"""
Match every HDF5 cell to the supplied GeoParquet mesh.

The HDF5 deliberately contains no geometry. This validation prevents an
incorrect mesh (for example a different resolution or simulation domain)
from silently producing a misleading animation.
"""
function match_mesh_cells(data, mesh_idx)
    mesh_rows = Vector{Int}(undef, data.n_cells)
    missing = String[]

    for i in 1:data.n_cells
        cid = normalise_cell_id(data.cell_ids[i])

        if haskey(mesh_idx, cid)
            mesh_rows[i] = mesh_idx[cid]
        else
            push!(missing, cid)
        end
    end

    if !isempty(missing)
        preview = join(missing[1:min(5, length(missing))], ", ")

        error(
            "$(length(missing)) HDF5 cells were not found in the supplied " *
            "GeoParquet mesh. First missing cell(s): $preview"
        )
    end

    return mesh_rows
end


# =============================================================================
# Animation
# =============================================================================

"""
Construct the static polygon geometry for the animation.

The polygons are created once. Their colour is then updated for each frame,
rather than destroying and recreating thousands of Makie polygon plots.
"""
function build_animation_polygons!(
    ax,
    data,
    mesh_df,
    mesh_rows,
    color_values,
    colormap,
    colorrange,
    local_units,
)
    cx = mean(data.lons)
    cy = mean(data.lats)

    plots = Vector{Any}(undef, data.n_cells)

    for i in 1:data.n_cells
        bnd = get_boundary(mesh_df, mesh_rows[i])

        xs, ys = project_boundary(
            bnd,
            cx,
            cy,
            local_units,
        )

        # No stroke: the A5 polygons share their actual boundary vertices,
        # so neighbouring cells abut without a distracting grid overlay.
        plots[i] = poly!(
            ax,
            Point2f.(zip(xs, ys));
            color = Observable(color_values[i]),
            colormap = colormap,
            colorrange = colorrange,
            strokewidth = 0,
        )
    end

    return plots
end


function update_animation_colours!(
    plots,
    values,
)
    for i in eachindex(plots)
        plots[i].color[] = values[i]
    end
end


function animate_h5(
    data,
    var_data,
    var_label,
    colormap,
    mesh_path,
    frame_step,
    output_dir,
    resolution,
    local_units,
)
    isnothing(mesh_path) &&
        error("--mesh is required when using --animate")

    isfile(mesh_path) ||
        error("Mesh file not found: $mesh_path")

    frame_step > 0 ||
        error("--frame-step must be greater than zero")

    println()
    println("Animation mode")
    println("  Mesh       : $mesh_path")
    println("  Frame step : $frame_step")
    println("  Output     : $output_dir")
    println("  Resolution : $(resolution[1]) × $(resolution[2]) px")
    println(
        "  Coordinates: " *
        (local_units ? "local East-North metres" : "longitude/latitude")
    )

    mesh_df, mesh_idx = load_parquet(mesh_path)

    println("  Mesh cells : $(nrow(mesh_df))")

    mesh_rows = match_mesh_cells(data, mesh_idx)

    # Use one fixed colour range for the entire animation. For physical
    # variables whose natural minimum is zero, retaining zero as the lower
    # bound makes the colour scale consistent and interpretable.
    finite_values = var_data[isfinite.(var_data)]

    isempty(finite_values) &&
        error("Selected variable '$var_label' contains no finite values")

    data_min = minimum(finite_values)
    data_max = maximum(finite_values)

    color_min = data_min

    if var_label in (
        "Water depth (m)",
        "Volume (m³)",
        "Saturation (0–1)",
        "Velocity (m/s)",
    )
        color_min = min(0.0, data_min)
    end

    color_max = data_max

    if color_max <= color_min
        color_max = color_min + 1e-6
    end

    colorrange = (color_min, color_max)

    println(
        "  Colour range: " *
        "$(round(color_min, digits=4)) – " *
        "$(round(color_max, digits=4))"
    )

    frames = collect(1:frame_step:data.n_frames)

    # Always include the final simulation frame, even when it is not an
    # exact multiple of frame_step.
    if frames[end] != data.n_frames
        push!(frames, data.n_frames)
    end

    mkpath(output_dir)

    fig = Figure(
        size = resolution,
        backgroundcolor = :white,
    )

    ax = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        xlabel = local_units ? "Easting (m)" : "Longitude",
        ylabel = local_units ? "Northing (m)" : "Latitude",
        title = "",
    )

    plots = build_animation_polygons!(
        ax,
        data,
        mesh_df,
        mesh_rows,
        var_data[frames[1], :],
        colormap,
        colorrange,
        local_units,
    )

    Colorbar(
        fig[1, 2],
        limits = colorrange,
        colormap = colormap,
        label = var_label,
    )

    # The map itself is the presentation focus; no cell outlines.
    hidedecorations!(
        ax,
        ticks = false,
        ticklabels = false,
        label = false,
    )
    hidespines!(ax)

    for (output_index, frame_idx) in enumerate(frames)
        values = var_data[frame_idx, :]

        if output_index > 1
            update_animation_colours!(plots, values)
        end

        t_s = data.times[frame_idx]
        t_hr = t_s / 3600.0

        ax.title =
            "$(var_label)  |  t = $(round(t_hr, digits=2)) h"

        output_path = joinpath(
            output_dir,
            @sprintf("frame_%06d.png", output_index),
        )

        CairoMakie.save(output_path, fig)

        if output_index == 1 ||
           output_index == length(frames) ||
           output_index % 10 == 0

            @printf(
                "  Frame %d/%d  |  simulation t = %.2f h  |  %s\n",
                output_index,
                length(frames),
                t_hr,
                output_path,
            )
        end
    end

    println()
    println(
        "Animation frames written to: " *
        "$(abspath(output_dir))"
    )
    println("  Frames: $(length(frames))")
    println("  First : frame_000001.png")
    println(
        "  Last  : " *
        @sprintf("frame_%06d.png", length(frames))
    )

    return nothing
end


# =============================================================================
# Interactive viewer
# =============================================================================

function interactive_view(
    data,
    cell_id,
    varname,
    time_req,
    resolution,
)
    var_data, var_label, colormap = get_var(data, varname)

    map_fi, map_t_actual = select_frame(
        data.times,
        time_req,
    )

    map_t_s = data.times[map_fi]
    map_t_hr = map_t_s / 3600.0

    if isnothing(time_req)
        map_title_suffix =
            "t = $(round(map_t_hr, digits=2)) h (final frame)"
    else
        map_title_suffix =
            "t = $(round(map_t_hr, digits=2)) h " *
            "($(round(map_t_s, digits=0))s)"
    end

    # -------------------------------------------------------------------------
    # Cell lookup
    # -------------------------------------------------------------------------

    cell_idx = nothing

    if !isnothing(cell_id)
        cell_idx = find_cell(
            data.cell_ids,
            cell_id,
        )

        if isnothing(cell_idx)
            @warn "Cell ID '$cell_id' not found in mesh. " *
                  "Available IDs (first 5): " *
                  "$(data.cell_ids[1:min(5, end)])"
        else
            println()
            println("Cell $cell_id (index $cell_idx):")

            println(
                "  Lon/Lat      : " *
                "$(round(data.lons[cell_idx], digits=5)), " *
                "$(round(data.lats[cell_idx], digits=5))",
            )

            println(
                "  Elevation    : " *
                "$(round(data.elevs[cell_idx], digits=2)) m",
            )

            println(
                "  Final depth  : " *
                "$(round(data.depth[end, cell_idx], digits=3)) m",
            )

            println(
                "  Final volume : " *
                "$(round(data.volume[end, cell_idx], digits=1)) m³",
            )

            println(
                "  Final sat    : " *
                "$(round(data.saturation[end, cell_idx], digits=3))",
            )

            peak_depth_idx =
                argmax(data.depth[:, cell_idx])

            println(
                "  Peak depth   : " *
                "$(round(maximum(data.depth[:, cell_idx]), digits=3)) m " *
                "at t=$(round(data.t_hr[peak_depth_idx], digits=2)) h",
            )

            peak_volume_idx =
                argmax(data.volume[:, cell_idx])

            println(
                "  Peak volume  : " *
                "$(round(maximum(data.volume[:, cell_idx]), digits=1)) m³ " *
                "at t=$(round(data.t_hr[peak_volume_idx], digits=2)) h",
            )
        end
    end

    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------

    fig_height = isnothing(cell_idx) ?
                 resolution[2] :
                 max(resolution[2], 950)

    fig = Figure(
        size = (
            resolution[1],
            fig_height,
        ),
    )

    ax_map = Axis(
        fig[1, 1],
        title = "$(var_label) — $map_title_suffix",
        xlabel = "Longitude",
        ylabel = "Latitude",
        aspect = DataAspect(),
    )

    map_vals = var_data[map_fi, :]

    finite_map = map_vals[isfinite.(map_vals)]

    isempty(finite_map) &&
        error("Selected variable '$varname' contains no finite values")

    vmin = minimum(finite_map)
    vmax = maximum(finite_map)

    if vmax <= vmin
        vmax = vmin + 1e-6
    end

    sc = scatter!(
        ax_map,
        data.lons,
        data.lats;
        color = map_vals,
        colormap = colormap,
        colorrange = (vmin, vmax),
        markersize = 8,
    )

    Colorbar(
        fig[1, 2],
        sc,
        label = var_label,
    )

    # Highlight queried cell.
    if !isnothing(cell_idx)
        scatter!(
            ax_map,
            [data.lons[cell_idx]],
            [data.lats[cell_idx]];
            color = :red,
            marker = :star5,
            markersize = 18,
        )

        text!(
            ax_map,
            data.lons[cell_idx],
            data.lats[cell_idx];
            text = "  $(cell_id[1:min(8, length(cell_id))])…",
            fontsize = 10,
            color = :red,
        )
    end

    # -------------------------------------------------------------------------
    # Cell time series
    # -------------------------------------------------------------------------

    if !isnothing(cell_idx)
        ax_ts = Axis(
            fig[2, 1:2],
            title = "Cell $cell_id — time series",
            xlabel = "Time (h)",
            ylabel = var_label,
        )

        cell_series = var_data[:, cell_idx]

        lines!(
            ax_ts,
            data.t_hr,
            cell_series;
            color = :steelblue,
            linewidth = 2,
        )

        scatter!(
            ax_ts,
            data.t_hr,
            cell_series;
            color = :steelblue,
            markersize = 5,
        )

        vlines!(
            ax_ts,
            [map_t_hr];
            color = :darkorange,
            linewidth = 1.5,
            linestyle = :dash,
        )

        peak_idx = argmax(cell_series)
        peak_y = cell_series[peak_idx]

        text_y = if isfinite(peak_y) && peak_y != 0
            peak_y * 0.05
        else
            0.0
        end

        text!(
            ax_ts,
            map_t_hr,
            text_y;
            text = " map\n $(round(map_t_s))s",
            fontsize = 9,
            color = :darkorange,
        )

        scatter!(
            ax_ts,
            [data.t_hr[peak_idx]],
            [cell_series[peak_idx]];
            color = :red,
            markersize = 10,
        )

        text!(
            ax_ts,
            data.t_hr[peak_idx],
            cell_series[peak_idx];
            text =
                " peak\n " *
                "$(round(cell_series[peak_idx], digits=2))",
            fontsize = 10,
            color = :red,
        )

        # Domain total on secondary y-axis.
        ax_dom = Axis(
            fig[2, 1:2],
            ylabel = "Domain total volume (×10⁶ m³)",
            yaxisposition = :right,
            ylabelcolor = :gray50,
            yticklabelcolor = :gray50,
        )

        hidespines!(ax_dom)
        hidexdecorations!(ax_dom)

        domain_vol = [
            sum(data.volume[fi, :])
            for fi in 1:data.n_frames
        ]

        lines!(
            ax_dom,
            data.t_hr,
            domain_vol ./ 1e6;
            color = :gray50,
            linewidth = 1.5,
            linestyle = :dash,
        )

        println()
        println(
            "Time series summary " *
            "($varname in cell $cell_id):"
        )
        println("(→ = map snapshot frame)")

        for (fi, t) in enumerate(data.t_hr)
            marker = fi == map_fi ? "→" : " "

            println(
                "$marker Frame $fi/$(data.n_frames) " *
                "— t = $(round(t, digits=3)) h " *
                "— $(round(cell_series[fi], digits=4))"
            )
        end
    end

    display(fig)

    println()
    println("Close the window to exit.")

    # Keep the interactive window alive until it is closed.
    for screen in copy(GLMakie.ALL_SCREENS)
        if screen.window_open[]
            wait(screen)
        end
    end
end


# =============================================================================
# Main
# =============================================================================

args = parse_args(ARGS)

data = load_h5(args.h5file)

if args.animate
    animate_h5(
        data,
        get_var(data, args.varname)...,
        args.mesh_path,
        args.frame_step,
        args.output_dir,
        args.resolution,
        args.local_units,
    )
else
    if !isnothing(args.mesh_path)
        @warn "--mesh was supplied without --animate; ignoring it."
    end

    if args.local_units
        @warn "--local-units currently applies to animation geometry only; " *
              "the interactive viewer uses HDF5 longitude/latitude coordinates."
    end

    interactive_view(
        data,
        args.cell_id,
        args.varname,
        args.time_req,
        args.resolution,
    )
end
