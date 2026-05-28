# benchmark_sim.jl — FloodA5 simulation performance benchmark
#
# Usage (from project root):
#   julia --threads auto benchmark_sim.jl [options]
#   julia --threads 1    benchmark_sim.jl --out results_1t.csv
#
# Options:
#   --config FILE   JSON config  (default: benchmark_config.json)
#   --out    FILE   CSV output   (default: benchmark_results.csv)
#   --warmup N      JIT warmup steps before timing (default: 10)
#   --steps  N      Timed steps per run            (default: 200)
#   --dt     S      Fixed timestep in seconds      (default: 30.0)
#   --help          Print this message
#
# World-age design note:
#   All functions that reference FloodModel symbols (A5Grid, SGSFlow, etc.)
#   are defined INSIDE the top-level `if abspath(PROGRAM_FILE) == @__FILE__`
#   block, AFTER include("FloodModel.jl") has run.  This guarantees they are
#   compiled in the same world as the symbols they call — no invokelatest needed.
#   Functions above the include (CLI parsing, config loading, helpers) touch
#   only Julia stdlib types and are safe to define at parse time.

using Dates, Statistics, Printf

# ---------------------------------------------------------------------------
# CLI parsing  (no FloodModel symbols — safe at parse time)
# ---------------------------------------------------------------------------

function _parse_bench_args(args)
    cfg = Dict{String,Any}(
        "config"  => "benchmark_config.json",
        "out"     => "benchmark_results.csv",
        "warmup"  => 10,
        "steps"   => 200,
        "dt"      => 30.0,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--help", "-h")
            println("""
benchmark_sim.jl — FloodA5 simulation performance benchmark

Usage:
  julia --threads N benchmark_sim.jl [options]

Options:
  --config FILE   JSON config  (default: benchmark_config.json)
  --out    FILE   CSV output   (default: benchmark_results.csv)
  --warmup N      JIT warmup steps before timing (default: 10)
  --steps  N      Timed steps per run            (default: 200)
  --dt     S      Fixed timestep in seconds      (default: 30.0)
  --help          Print this message
""")
            exit(0)
        elseif a == "--config" && i < length(args)
            cfg["config"] = args[i+1];  i += 2
        elseif a == "--out" && i < length(args)
            cfg["out"] = args[i+1];     i += 2
        elseif a == "--warmup" && i < length(args)
            cfg["warmup"] = parse(Int, args[i+1]);    i += 2
        elseif a == "--steps" && i < length(args)
            cfg["steps"]  = parse(Int, args[i+1]);    i += 2
        elseif a == "--dt" && i < length(args)
            cfg["dt"] = parse(Float64, args[i+1]);    i += 2
        else
            @warn "Unknown argument: $a"
            i += 1
        end
    end
    return cfg
end

# ---------------------------------------------------------------------------
# Default run list  (no FloodModel symbols — safe at parse time)
# ---------------------------------------------------------------------------

function _default_runs()
    return [
        Dict{String,Any}(
            "name"       => "Flat res14 standard",
            "meshload"   => "test/flat_mesh_res14.parquet",
            "flow_model" => "standard",
            "rainfall"   => 50.0,
            "rainpoints" => String[]),
        Dict{String,Any}(
            "name"       => "Carlisle res14 standard",
            "meshload"   => "test/carlisle/carlisle_mesh14_standard.parquet",
            "flow_model" => "standard",
            "rainfall"   => 50.0,
            "rainpoints" => String[]),
        Dict{String,Any}(
            "name"       => "Carlisle res17 standard",
            "meshload"   => "test/carlisle/carlisle_mesh17_standard.parquet",
            "flow_model" => "standard",
            "rainfall"   => 0.0,
            "rainpoints" => ["54.908,-2.896,1000.0"]),
        Dict{String,Any}(
            "name"       => "Carlisle res14 SGS",
            "meshload"   => "test/carlisle/carlisle_mesh14_sgs.parquet",
            "flow_model" => "sgs",
            "rainfall"   => 50.0,
            "rainpoints" => String[]),
    ]
end

# ---------------------------------------------------------------------------
# Config loading  (no FloodModel symbols — safe at parse time)
# ---------------------------------------------------------------------------

function _load_runs(cfg)
    path = cfg["config"]
    if isfile(path)
        raw     = read(path, String)
        cleaned = replace(raw, r"//[^\n]*" => "")
        try
            data = Main.JSON3.read(cleaned)
            if haskey(data, :runs)
                return [Dict{String,Any}(string(k) => v for (k,v) in r)
                        for r in data[:runs]]
            end
        catch e
            @warn "Could not parse $path ($e) — using default run list"
        end
    else
        @info "No $(path) found — using built-in default run list"
    end
    return _default_runs()
end

# ---------------------------------------------------------------------------
# SI formatter  (no FloodModel symbols — safe at parse time)
# ---------------------------------------------------------------------------

function format_si(n::Int)::String
    n >= 1_000_000_000 && return @sprintf("%.2fG", n / 1e9)
    n >= 1_000_000     && return @sprintf("%.2fM", n / 1e6)
    n >= 1_000         && return @sprintf("%.1fk", n / 1e3)
    return string(n)
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
# Everything below runs only when this file is the main script.
# FloodModel.jl is included FIRST; all functions that call FloodModel symbols
# are defined AFTER the include, so they are compiled in the correct world.
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__

    # ── Step 1: parse CLI and load run config (no FloodModel needed yet) ──
    cfg  = _parse_bench_args(ARGS)
    runs = _load_runs(cfg)

    n_threads = Threads.nthreads()
    dt        = Float64(cfg["dt"])
    n_warmup  = Int(cfg["warmup"])
    n_steps   = Int(cfg["steps"])
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    hostname  = try chomp(read(`hostname`, String)) catch _ "unknown" end

    println("\n", "="^72)
    println("FloodA5 Simulation Benchmark  —  ", timestamp)
    @printf("  Julia threads : %d\n", n_threads)
    @printf("  Warmup steps  : %d\n", n_warmup)
    @printf("  Timed steps   : %d\n", n_steps)
    @printf("  Fixed dt      : %.0f s\n", dt)
    @printf("  Host          : %s\n", hostname)
    println("="^72)

    # ── Step 2: include FloodModel.jl at top level ────────────────────────
    # This defines A5Grid, SGSFlow, StandardFlow, initialise_flow_model,
    # step_standard!, step_sgs!, _find_nearest_cell, etc. in the current world.
    project_dir      = dirname(abspath(@__FILE__))
    flood_model_path = joinpath(project_dir, "FloodModel.jl")
    isfile(flood_model_path) ||
        error("FloodModel.jl not found at $flood_model_path")

    # Set env flag BEFORE include so FloodModel's entry-point guard sees it
    # and skips calling its own main(ARGS).  Cleared immediately after.
    ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
    @info "Loading FloodModel.jl..."
    include(flood_model_path)
    @info "FloodModel.jl loaded"
    delete!(ENV, "FLOODMODEL_INCLUDE_ONLY")

    # Suppress the step_standard! debug logging (first-2-calls counter) so it
    # doesn't fire during warmup and inflate timing or clutter output.
    _step_debug_count[] = typemax(Int)

    # ── Step 3: define benchmark functions IN THIS BLOCK, after include ───
    # Any function defined here is compiled AFTER the include, so it lives
    # in the same world as the FloodModel symbols it calls.  This is the key
    # fix: function definitions at the top of the file are compiled at parse
    # time (before include runs) regardless of where include appears.

    """
        _time_steps(state, flow_model, dt, rainfall_rate, n_warmup, n_timed)

    Run `n_warmup + n_timed` simulation steps at a fixed `dt`.

    Returns a named tuple with timing and flow-model diagnostic fields:

    **Timing**
    - `median_ms`, `min_ms`, `max_ms` — per-step wall-clock over timed steps (ms)
    - `first_step_ms` — wall time of step 1 including JIT compilation

    **Flow model metrics** (computed over the `n_timed` timed steps only)
    - `max_depth_m`       — peak water depth seen in any cell at any timed step (m)
    - `max_wet_area_m2`   — peak flood extent: Σ(saturation_i × cell_area_i) (m²),
                            SGS-correct (uses hypsometric wetted fraction per cell)
    - `final_wet_area_m2` — flood extent at the last timed step (m²)
    - `final_vol_m3`      — total stored volume at the last timed step (m³)
    - `mb_err_m3`         — mass balance error over timed steps: expected input
                            minus actual domain-volume change (m³); ~0 = perfect
    - `mb_err_pct`        — same as a percentage of total injected volume
    """
    function _time_steps(state, flow_model::Symbol,
                         dt::Float64, rainfall_rate::Float64,
                         n_warmup::Int, n_timed::Int)

        # Pre-sync water_depth before first step (mirrors run_simulation!)
        if flow_model == :standard
            for i in eachindex(state.cell_ids)
                state.cell_area[i] >= 1.0 &&
                    (state.water_depth[i] = state.volume[i] / state.cell_area[i])
            end
        end

        step_fn!      = flow_model == :standard ? step_standard! : step_sgs!
        times         = Vector{Float64}(undef, n_timed)
        first_step_ms = 0.0

        # Flow model tracking (timed steps only)
        max_depth_m      = 0.0
        max_wet_area_m2  = 0.0
        final_wet_area_m2 = 0.0
        vol_start        = sum(state.volume)   # volume at start of timed phase
        input_vol_timed  = 0.0                 # cumulative injection during timed steps

        # Helper: SGS-correct wetted area (fractional for SGS, binary for standard)
        function _wet_area(state)
            if isempty(state.sgs_tables)
                return sum(state.cell_area[i]
                           for i in eachindex(state.cell_ids)
                           if state.water_depth[i] > 1e-4; init=0.0)
            else
                return sum(begin
                    tbl    = state.sgs_tables[i]
                    w_area = wetted_area_from_wse(tbl,
                                 wse_from_volume(tbl, state.volume[i]))
                    clamp(w_area, 0.0, tbl.cell_area)
                end for i in eachindex(state.cell_ids); init=0.0)
            end
        end

        for pass in 1:(n_warmup + n_timed)
            if rainfall_rate > 0.0
                for i in eachindex(state.cell_ids)
                    state.volume[i] += rainfall_rate * dt * state.cell_area[i]
                end
            end

            t0         = time_ns()
            step_fn!(state, dt)
            elapsed_ms = (time_ns() - t0) / 1.0e6

            if pass == 1
                first_step_ms = elapsed_ms
            end

            if pass > n_warmup
                times[pass - n_warmup] = elapsed_ms

                # Accumulate input volume for this timed step
                input_vol_timed += rainfall_rate * dt *
                    sum(a for a in state.cell_area if a >= 1.0; init=0.0)

                # Track peak depth
                for d in state.water_depth
                    isfinite(d) && d > max_depth_m && (max_depth_m = d)
                end

                # Track peak and final wetted area
                wa = _wet_area(state)
                wa > max_wet_area_m2 && (max_wet_area_m2 = wa)
                if pass == n_warmup + n_timed
                    final_wet_area_m2 = wa
                end
            end
        end

        final_vol_m3 = sum(state.volume)
        vol_change   = final_vol_m3 - vol_start
        mb_err_m3    = input_vol_timed - vol_change
        mb_err_pct   = input_vol_timed > 0.0 ?
                       100.0 * mb_err_m3 / input_vol_timed : 0.0

        return (
            median_ms        = median(times),
            min_ms           = minimum(times),
            max_ms           = maximum(times),
            first_step_ms    = first_step_ms,
            max_depth_m      = max_depth_m,
            max_wet_area_m2  = max_wet_area_m2,
            final_wet_area_m2 = final_wet_area_m2,
            final_vol_m3     = final_vol_m3,
            mb_err_m3        = mb_err_m3,
            mb_err_pct       = mb_err_pct,
        )
    end

    # ── Step 4: open CSV output (append if file exists, create if not) ───────
    csv_path   = cfg["out"]
    out_dir    = dirname(abspath(csv_path))
    out_dir != "" && !isdir(out_dir) && mkpath(out_dir)
    file_exists = isfile(csv_path)
    csv_io      = open(csv_path, file_exists ? "a" : "w")
    # Write header only when creating a new file; appending skips it so that
    # multiple julia --threads N runs all land in the same CSV cleanly.
    if !file_exists
        println(csv_io,
            "timestamp,hostname,threads,run_name,flow_model,n_cells,n_edges," *
            "mesh_res,dt_s,warmup_steps,timed_steps," *
            "median_ms_per_step,min_ms_per_step,max_ms_per_step," *
            "first_step_ms,jit_overhead_x," *
            "cells_per_second,steps_per_sim_hour," *
            "mesh_load_s,model_init_s,startup_s,sim_only_s,total_wall_s," *
            "state_memory_mb,bytes_per_cell," *
            "mb_err_m3,mb_err_pct," *
            "max_depth_m,max_wet_area_m2,final_wet_area_m2,final_vol_m3")
    end

    # ── Step 5: run each benchmark ─────────────────────────────────────────
    for run in runs
        name       = get(run, "name",       "unnamed")
        meshfile   = get(run, "meshload",   "")   # key is "meshload" (matches launch.json convention)
        flow_model = Symbol(get(run, "flow_model", "standard"))
        rainfall   = Float64(get(run, "rainfall",  0.0))
        rainpoints = get(run, "rainpoints", String[])

        if !isfile(meshfile)
            @warn "[$name] Mesh not found: $meshfile — skipping"
            continue
        end

        println("\n── $name ──")
        @printf("  Loading mesh from %s\n", meshfile)
        t_load  = @elapsed mesh = A5Grid.load_mesh_geoparquet(meshfile)
        n_cells = length(mesh.cells)
        res     = mesh.resolution
        @printf("  %d cells (res %d) loaded in %.2f s\n", n_cells, res, t_load)

        @printf("  Initialising %s flow model...\n", flow_model)
        GC.gc()   # clean baseline before measuring state footprint
        mem_before = Base.gc_live_bytes()
        t_init = @elapsed begin
            method = flow_model == :sgs ? SGSFlow() : StandardFlow()
            state  = initialise_flow_model(mesh, method)
        end
        GC.gc()   # collect any temporaries from init
        mem_after      = Base.gc_live_bytes()
        state_mem_mb   = (mem_after - mem_before) / 1024^2
        bytes_per_cell = state_mem_mb * 1024^2 / n_cells
        n_edges        = state.edges.n_edges
        @printf("  %d edges ready in %.2f s  (state footprint: %.1f MB, %.0f B/cell)\n",
                n_edges, t_init, state_mem_mb, bytes_per_cell)

        # Seed rainpoints: inject warmup-equivalent volume so cells are wet
        for rp in rainpoints
            parts = split(string(rp), ",")
            if length(parts) == 3
                lat       = parse(Float64, strip(parts[1]))
                lon       = parse(Float64, strip(parts[2]))
                rate_mmhr = parse(Float64, strip(parts[3]))
                # _find_nearest_cell signature: (mesh, lon, lat) → (idx, cell_id, dist_m)
                idx, _, dist_m = _find_nearest_cell(mesh, lon, lat)
                if dist_m < 2000.0
                    seed_vol = (rate_mmhr / 3_600_000.0) *
                               state.cell_area[idx] * n_warmup * dt
                    state.volume[idx] += seed_vol
                    @printf("  Seeded cell %d with %.1f m³\n", idx, seed_vol)
                else
                    @warn "Rainpoint (lat=$lat, lon=$lon) is $(round(dist_m/1000,digits=1)) km from mesh — not seeding"
                end
            end
        end

        rainfall_rate = rainfall / 3_600_000.0

        @printf("  Timing (%d warmup + %d steps, dt=%.0f s)...\n",
                n_warmup, n_steps, dt)
        flush(stdout)

        r = _time_steps(state, flow_model, dt, rainfall_rate, n_warmup, n_steps)

        # Pure computation time: median step × n timed steps (no warmup/JIT/IO)
        sim_only_s      = n_steps * r.median_ms / 1000.0

        cells_per_sec   = round(Int, n_cells / (r.median_ms / 1000.0))
        steps_per_simhr = 3600.0 / dt / (r.median_ms / 1000.0)
        jit_overhead_x  = r.median_ms > 0.0 ? r.first_step_ms / r.median_ms : 0.0
        startup_s       = t_load + t_init
        total_wall_s    = startup_s + sim_only_s

        @printf("  step timing:  median %8.3f ms   min %8.3f   max %8.3f\n",
                r.median_ms, r.min_ms, r.max_ms)
        @printf("  first step:   %8.3f ms  (JIT overhead: %.1fx median)\n",
                r.first_step_ms, jit_overhead_x)
        @printf("  throughput:   %-10s cells/s   %.1f steps per sim-hour\n",
                format_si(cells_per_sec), steps_per_simhr)
        @printf("  sim only:     %.3f s for %d steps  (%.1f%% of total wall)\n",
                sim_only_s, n_steps, 100.0 * sim_only_s / max(total_wall_s, 1e-9))
        @printf("  startup:      %.2f s  (load %.2f s + init %.2f s)\n",
                startup_s, t_load, t_init)
        @printf("  total wall:   %.2f s\n", total_wall_s)
        @printf("  memory:       %.1f MB state  (%.0f B/cell)\n",
                state_mem_mb, bytes_per_cell)
        @printf("  mass balance: err = %+.4f m³  (%+.6f%%)\n",
                r.mb_err_m3, r.mb_err_pct)
        @printf("  flood extent: max = %.0f m²  final = %.0f m²  (%.2f km²)\n",
                r.max_wet_area_m2, r.final_wet_area_m2, r.max_wet_area_m2 / 1e6)
        @printf("  max depth:    %.3f m   final domain vol: %.1f m³\n",
                r.max_depth_m, r.final_vol_m3)

        println(csv_io, join([
            timestamp, hostname, n_threads,
            "\"$name\"", flow_model,
            n_cells, n_edges, res, dt,
            n_warmup, n_steps,
            round(r.median_ms,       digits=4),
            round(r.min_ms,          digits=4),
            round(r.max_ms,          digits=4),
            round(r.first_step_ms,   digits=4),
            round(jit_overhead_x,    digits=2),
            cells_per_sec,
            round(steps_per_simhr,   digits=2),
            round(t_load,            digits=3),
            round(t_init,            digits=3),
            round(startup_s,         digits=3),
            round(sim_only_s,        digits=3),
            round(total_wall_s,      digits=3),
            round(state_mem_mb,      digits=2),
            round(bytes_per_cell,    digits=0),
            round(r.mb_err_m3,       digits=4),
            round(r.mb_err_pct,      digits=6),
            round(r.max_depth_m,     digits=4),
            round(r.max_wet_area_m2, digits=0),
            round(r.final_wet_area_m2, digits=0),
            round(r.final_vol_m3,    digits=2),
        ], ","))
        flush(csv_io)
    end

    close(csv_io)
    println("\n", "="^72)
    println("Results written to: $csv_path")
    println("="^72, "\n")

end  # if abspath(PROGRAM_FILE) == @__FILE__
