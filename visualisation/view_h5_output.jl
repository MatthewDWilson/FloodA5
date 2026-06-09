"""
view_h5_output.jl — Interactive viewer and query tool for FloodA5 HDF5 simulation output.

Usage:
    julia --project=. view_h5_output.jl <h5file> [--cell CELL_ID] [--var VARNAME] [--time SECONDS]

Arguments:
    <h5file>         Path to FloodA5 .h5 output file
    --cell ID        Hex cell ID to plot time series for (e.g. 6345ddac19800000)
    --var NAME       Variable to plot: depth (default), volume, saturation, velocity
    --time SECONDS   Simulation time in seconds for the map snapshot (default: final frame).
                     If the requested time is not exactly in the output, the nearest frame
                     is used and a warning is printed.  If beyond the simulation end, the
                     final frame is used.

Examples:
    # Map at final frame + time series for a cell
    julia --project=. view_h5_output.jl test/carlisle/carlisle_sgs_ra_res18.h5 \\
        --cell 6345ddac19800000

    # Map at t=7200s (2h) with volume variable
    julia --project=. view_h5_output.jl sim.h5 --time 7200 --var volume \\
        --cell 6345ddac19800000
"""

using HDF5, GLMakie, Printf, Statistics

# ── CLI parsing ────────────────────────────────────────────────────────────────
function parse_args(args)
    length(args) < 1 && error("Usage: view_h5_output.jl <h5file> [--cell ID] [--var NAME] [--time SECONDS]")
    h5file      = args[1]
    cell_id     = nothing
    varname     = "depth"
    time_req    = nothing   # requested simulation time in seconds (nothing = final frame)
    i = 2
    while i <= length(args)
        if args[i] == "--cell" && i < length(args)
            cell_id = lowercase(lstrip(args[i+1], '0'))
            i += 2
        elseif args[i] == "--var" && i < length(args)
            varname = args[i+1]
            i += 2
        elseif args[i] == "--time" && i < length(args)
            time_req = tryparse(Float64, args[i+1])
            isnothing(time_req) && error("--time must be a number (seconds), got: $(args[i+1])")
            i += 2
        else
            i += 1
        end
    end
    return h5file, cell_id, varname, time_req
end

# ── Frame selection with nearest-match and warnings ───────────────────────────
function select_frame(times::Vector{Float64}, time_req)
    # No time requested → final frame
    isnothing(time_req) && return length(times), nothing

    t_end = times[end]

    # Beyond simulation end
    if time_req > t_end
        @warn "Requested time $(time_req)s is beyond the simulation end ($(t_end)s). Showing final frame."
        return length(times), t_end
    end

    # Exact or nearest match
    diffs = abs.(times .- time_req)
    fi    = argmin(diffs)
    nearest = times[fi]

    if diffs[fi] > 0.5   # more than 0.5s away — warn
        @warn "Requested time $(time_req)s not in output. Nearest frame: $(round(nearest, digits=1))s (Δ=$(round(diffs[fi], digits=1))s)."
    end

    return fi, nearest
end

# ── HDF5 reading ───────────────────────────────────────────────────────────────
function load_h5(h5file::String)
    isfile(h5file) || error("File not found: $h5file")
    println("Reading: $h5file")

    h5open(h5file, "r") do f
        cell_ids = read(f["mesh/cell_ids"])
        lons     = read(f["mesh/center_lons"])
        lats     = read(f["mesh/center_lats"])
        elevs    = read(f["mesh/elevations"])

        frames   = sort(keys(f["frames"]))
        n_frames = length(frames)
        n_cells  = length(cell_ids)

        println("  Cells   : $n_cells")
        println("  Frames  : $n_frames")

        times    = Float64[]
        depth_ts = Matrix{Float64}(undef, n_frames, n_cells)
        vol_ts   = Matrix{Float64}(undef, n_frames, n_cells)
        sat_ts   = Matrix{Float64}(undef, n_frames, n_cells)
        vel_ts   = Matrix{Float64}(undef, n_frames, n_cells)

        for (fi, fr) in enumerate(frames)
            g = f["frames/$fr"]
            push!(times, read(g["t"]))
            depth_ts[fi, :] = read(g["water_depth"])
            vol_ts[fi, :]   = read(g["volume"])
            sat_ts[fi, :]   = read(g["saturation"])
            vel_ts[fi, :]   = read(g["velocity"])
        end

        t_hr = times ./ 3600.0

        # Summary stats at final frame
        final_depth = depth_ts[end, :]
        final_vol   = vol_ts[end, :]
        println("  Duration: $(round(times[end]/3600, digits=2)) h")
        println("  Final domain volume: $(round(sum(final_vol)/1e6, digits=3)) × 10⁶ m³")
        println("  Final max depth    : $(round(maximum(final_depth), digits=3)) m")
        println("  Final mean depth   : $(round(mean(final_depth[final_depth .> 0.001]), digits=3)) m (wet cells)")

        return (cell_ids=cell_ids, lons=lons, lats=lats, elevs=elevs,
                times=times, t_hr=t_hr,
                depth=depth_ts, volume=vol_ts, saturation=sat_ts, velocity=vel_ts,
                n_cells=n_cells, n_frames=n_frames)
    end
end

# ── Cell lookup ────────────────────────────────────────────────────────────────
function find_cell(cell_ids::Vector, target_id::String)
    # Normalise: zero-pad to 16 chars, lowercase
    norm(s) = lowercase(lpad(lstrip(string(s), '0'), 16, '0'))
    t = norm(target_id)
    for (i, cid) in enumerate(cell_ids)
        norm(cid) == t && return i
    end
    return nothing
end

# ── Variable selector ──────────────────────────────────────────────────────────
function get_var(data, varname::String)
    varname == "depth"      && return data.depth,      "Water depth (m)",         "turbo"
    varname == "volume"     && return data.volume,     "Volume (m³)",             "turbo"
    varname == "saturation" && return data.saturation, "Saturation (0–1)",        "blues"
    varname == "velocity"   && return data.velocity,   "Velocity (m/s)",          "plasma"
    error("Unknown variable '$varname'. Choose: depth, volume, saturation, velocity")
end

# ── Main ───────────────────────────────────────────────────────────────────────
h5file, cell_id, varname, time_req = parse_args(ARGS)
data = load_h5(h5file)
var_data, var_label, colormap = get_var(data, varname)

# Select map frame
map_fi, map_t_actual = select_frame(data.times, time_req)
map_t_s  = data.times[map_fi]
map_t_hr = map_t_s / 3600.0
if isnothing(time_req)
    map_title_suffix = "t = $(round(map_t_hr, digits=2)) h (final frame)"
else
    map_title_suffix = "t = $(round(map_t_hr, digits=2)) h ($(round(map_t_s, digits=0))s)"
end

# ── Cell time series ───────────────────────────────────────────────────────────
cell_idx = nothing
if !isnothing(cell_id)
    cell_idx = find_cell(data.cell_ids, cell_id)
    if isnothing(cell_idx)
        @warn "Cell ID '$cell_id' not found in mesh. Available IDs (first 5): $(data.cell_ids[1:min(5,end)])"
    else
        println("\nCell $cell_id (index $cell_idx):")
        println("  Lon/Lat  : $(round(data.lons[cell_idx],digits=5)), $(round(data.lats[cell_idx],digits=5))")
        println("  Elevation: $(round(data.elevs[cell_idx],digits=2)) m")
        println("  Final depth : $(round(data.depth[end, cell_idx],   digits=3)) m")
        println("  Final volume: $(round(data.volume[end, cell_idx],  digits=1)) m³")
        println("  Final sat   : $(round(data.saturation[end, cell_idx], digits=3))")
        println("  Peak depth  : $(round(maximum(data.depth[:, cell_idx]), digits=3)) m at t=$(round(data.t_hr[argmax(data.depth[:, cell_idx])], digits=2)) h")
        println("  Peak volume : $(round(maximum(data.volume[:, cell_idx]), digits=1)) m³ at t=$(round(data.t_hr[argmax(data.volume[:, cell_idx])], digits=2)) h")
    end
end

# ── Build figure ───────────────────────────────────────────────────────────────
fig = Figure(size = (1400, isnothing(cell_idx) ? 700 : 950))

# Top row: domain scatter map at selected frame
ax_map = Axis(fig[1, 1],
    title  = "$(var_label) — $map_title_suffix",
    xlabel = "Longitude", ylabel = "Latitude",
    aspect = DataAspect())

map_vals = var_data[map_fi, :]
vmin, vmax = minimum(map_vals), maximum(map_vals)
vmax <= vmin && (vmax = vmin + 1e-6)

sc = scatter!(ax_map, data.lons, data.lats;
    color      = map_vals,
    colormap   = colormap,
    colorrange = (vmin, vmax),
    markersize = 8)
Colorbar(fig[1, 2], sc, label = var_label)

# Highlight queried cell
if !isnothing(cell_idx)
    scatter!(ax_map, [data.lons[cell_idx]], [data.lats[cell_idx]];
        color = :red, marker = :star5, markersize = 18)
    text!(ax_map, data.lons[cell_idx], data.lats[cell_idx];
        text = "  $(cell_id[1:8])…", fontsize = 10, color = :red)
end

# Bottom row: time series for queried cell
if !isnothing(cell_idx)
    ax_ts = Axis(fig[2, 1:2],
        title  = "Cell $cell_id — time series",
        xlabel = "Time (h)", ylabel = var_label)

    cell_series = var_data[:, cell_idx]
    lines!(ax_ts, data.t_hr, cell_series; color = :steelblue, linewidth = 2)
    scatter!(ax_ts, data.t_hr, cell_series; color = :steelblue, markersize = 5)

    # Vertical line marking the map snapshot time
    vlines!(ax_ts, [map_t_hr]; color = :darkorange, linewidth = 1.5, linestyle = :dash)
    text!(ax_ts, map_t_hr, maximum(cell_series) * 0.05;
        text = " map\n $(round(map_t_s))s", fontsize = 9, color = :darkorange)

    # Annotate peak
    pk = argmax(cell_series)
    scatter!(ax_ts, [data.t_hr[pk]], [cell_series[pk]];
        color = :red, markersize = 10)
    text!(ax_ts, data.t_hr[pk], cell_series[pk];
        text = " peak\n $(round(cell_series[pk], digits=2))",
        fontsize = 10, color = :red)

    # Domain total on secondary y (volume or depth mean) for context
    ax_dom = Axis(fig[2, 1:2],
        ylabel    = "Domain total volume (×10⁶ m³)",
        yaxisposition = :right,
        ylabelcolor   = :gray50,
        yticklabelcolor = :gray50)
    hidespines!(ax_dom)
    hidexdecorations!(ax_dom)
    domain_vol = [sum(data.volume[fi, :]) for fi in 1:data.n_frames]
    lines!(ax_dom, data.t_hr, domain_vol ./ 1e6;
        color = :gray50, linewidth = 1.5, linestyle = :dash)

    println("\nTime series summary ($varname in cell $cell_id):")
    println("  (→ = map snapshot frame)")
    for (fi, t) in enumerate(data.t_hr)
        marker = fi == map_fi ? "→" : " "
        @printf("  %s t=%6.2f h  %-12s = %10.3f   domain_vol = %.3f ×10⁶ m³\n",
                marker, t, varname, cell_series[fi], domain_vol[fi]/1e6)
    end
end

display(fig)
println("\nClose the window to exit.")
wait(display(fig))   # keep window open until user closes
