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

    Returns `(median_ms, min_ms, max_ms, first_step_ms)` where:
    - `median/min/max_ms` are per-step wall-clock times over the `n_timed` timed steps.
    - `first_step_ms` is the wall-clock time of the very first step (pass 1), which
      includes JIT compilation cost.  Comparing this to `median_ms` gives the JIT
      overhead multiplier.
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
                first_step_ms = elapsed_ms   # JIT-inclusive first call
            end
            if pass > n_warmup
                times[pass - n_warmup] = elapsed_ms
            end
        end

        return median(times), minimum(times), maximum(times), first_step_ms
    end

    # ── Step 4: open CSV output ────────────────────────────────────────────
    csv_path = cfg["out"]
    out_dir  = dirname(abspath(csv_path))
    out_dir != "" && !isdir(out_dir) && mkpath(out_dir)
    csv_io   = open(csv_path, "w")
    println(csv_io,
        "timestamp,hostname,threads,run_name,flow_model,n_cells,n_edges," *
        "mesh_res,dt_s,warmup_steps,timed_steps," *
        "median_ms_per_step,min_ms_per_step,max_ms_per_step," *
        "first_step_ms,jit_overhead_x," *
        "cells_per_second,steps_per_sim_hour," *
        "mesh_load_s,model_init_s,startup_s,sim_only_s,total_wall_s," *
        "state_memory_mb,bytes_per_cell")

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

        med_ms, min_ms, max_ms, first_step_ms =
            _time_steps(state, flow_model, dt, rainfall_rate, n_warmup, n_steps)

        # Pure computation time: median step time × number of timed steps.
        # This excludes warmup (JIT), I/O, and mesh initialisation — it is the
        # best estimate of how long the solver itself would take in production
        # for a run of the same length.
        sim_only_s      = n_steps * med_ms / 1000.0

        cells_per_sec   = round(Int, n_cells / (med_ms / 1000.0))
        steps_per_simhr = 3600.0 / dt / (med_ms / 1000.0)
        jit_overhead_x  = med_ms > 0.0 ? first_step_ms / med_ms : 0.0
        startup_s       = t_load + t_init
        total_wall_s    = startup_s + sim_only_s

        @printf("  step timing:  median %8.3f ms   min %8.3f   max %8.3f\n",
                med_ms, min_ms, max_ms)
        @printf("  first step:   %8.3f ms  (JIT overhead: %.1fx median)\n",
                first_step_ms, jit_overhead_x)
        @printf("  throughput:   %-10s cells/s   %.1f steps per sim-hour\n",
                format_si(cells_per_sec), steps_per_simhr)
        @printf("  sim only:     %.3f s for %d steps  (%.1f%% of total wall)\n",
                sim_only_s, n_steps, 100.0 * sim_only_s / max(total_wall_s, 1e-9))
        @printf("  startup:      %.2f s  (load %.2f s + init %.2f s)\n",
                startup_s, t_load, t_init)
        @printf("  total wall:   %.2f s\n", total_wall_s)
        @printf("  memory:       %.1f MB state  (%.0f B/cell)\n",
                state_mem_mb, bytes_per_cell)

        println(csv_io, join([
            timestamp, hostname, n_threads,
            "\"$name\"", flow_model,
            n_cells, n_edges, res, dt,
            n_warmup, n_steps,
            round(med_ms,          digits=4),
            round(min_ms,          digits=4),
            round(max_ms,          digits=4),
            round(first_step_ms,   digits=4),
            round(jit_overhead_x,  digits=2),
            cells_per_sec,
            round(steps_per_simhr, digits=2),
            round(t_load,          digits=3),
            round(t_init,          digits=3),
            round(startup_s,       digits=3),
            round(sim_only_s,      digits=3),
            round(total_wall_s,    digits=3),
            round(state_mem_mb,    digits=2),
            round(bytes_per_cell,  digits=0),
        ], ","))
        flush(csv_io)
    end

    close(csv_io)
    println("\n", "="^72)
    println("Results written to: $csv_path")
    println("="^72, "\n")

end  # if abspath(PROGRAM_FILE) == @__FILE__
