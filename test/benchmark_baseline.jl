# test/benchmark_baseline.jl
# --------------------------
# Performance baseline for the step functions (WP0, 2026-10). Freeze these
# numbers BEFORE the WP4 performance refactor so every later change can be
# measured against them (FloodA5_Review_and_ForwardPlan_2026-10.md §7.4).
#
# Uses an in-process Cartesian mesh (no data files, no Python mesh bridge) and
# times step_standard! directly at a fixed dt after a warm-up, so the numbers
# reflect solver cost only (no CFL, sources, logging or I/O). Reports
# ns per cell·step, million cell·steps per second, and bytes allocated per step.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t 1 test/benchmark_baseline.jl
#   FLOODA5_NO_VIS=1 julia --project=. -t 8 test/benchmark_baseline.jl
#   FLOODA5_NO_VIS=1 julia --project=. -t 8 test/benchmark_baseline.jl --n 400 --steps 200
#
# For end-to-end numbers (including SGS, CFL, sources and HDF5 output) use the
# CLI timing summary on the Carlisle res-18 meshes — see the plan §7.4.

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Logging, Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : parse(Int, ARGS[i+1]))
const NSIDE  = _arg("--n", 200)       # grid is NSIDE × NSIDE cells
const NSTEPS = _arg("--steps", 100)

function _bench_state(; wet::Symbol, kwargs...)
    dx = 10.0
    dd = rad2deg(dx / A5Grid._EARTH_R)
    mesh = generate_cartesian_mesh_bounds(0.0, 0.0, (NSIDE - 0.5) * dd, (NSIDE - 0.5) * dd,
                                          dx; lat_ref = 0.0)
    n = length(mesh.cells)
    mesh.static_vars["elevation"] = [0.0005 * dx * (mesh.cells[i].center_lon / dd) +
                                     0.01 * sin(0.37 * i) for i in 1:n]
    st = initialise_flow_model(mesh, StandardFlow(); kwargs...)
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    if wet === :all
        for i in 1:n                                   # 0.3 m everywhere ("rain run")
            st.volume[i] = 0.3 * st.cell_area[i]; st.water_depth[i] = 0.3
        end
    else                                               # ~1% wet patch ("point source")
        r = max(1, round(Int, 0.05 * NSIDE))
        for i in 1:n
            row, col = divrem(i - 1, NSIDE)
            if abs(row - div(NSIDE, 2)) <= r && abs(col - div(NSIDE, 2)) <= r
                st.volume[i] = 0.5 * st.cell_area[i]; st.water_depth[i] = 0.5
            end
        end
    end
    return st
end

function _time_steps!(st, nsteps; dt = 1.0)
    for _ in 1:5; step_standard!(st, dt); end          # warm-up / compile
    GC.gc()
    bytes = @allocated step_standard!(st, dt)
    t0 = time_ns()
    for _ in 1:nsteps; step_standard!(st, dt); end
    el = (time_ns() - t0) / 1e9
    return el, bytes
end

configs = [
    ("std legacy (default)", (;)),
    ("std WLSQ corrected",   (; gradient_correction = true)),
    ("std diamond",          (; gradient_correction = true, face_flux_method = :diamond)),
    ("std diamond + cell",   (; gradient_correction = true, face_flux_method = :diamond,
                                momentum_model = :cell)),
]

println("FloodA5 step benchmark — $(NSIDE)×$(NSIDE) Cartesian cells, $NSTEPS steps, " *
        "threads = $(Threads.nthreads()), Julia $(VERSION)")
@printf("%-24s %-6s %12s %12s %14s\n", "configuration", "wet", "ns/cell·step", "M cell·st/s", "alloc/step")
with_logger(ConsoleLogger(stderr, Logging.Warn)) do
    for (label, kw) in configs, wet in (:all, :patch)
        st = _bench_state(; wet = wet, kw...)
        n  = length(st.cell_ids)
        el, bytes = _time_steps!(st, NSTEPS)
        ns = el * 1e9 / (n * NSTEPS)
        @printf("%-24s %-6s %12.1f %12.2f %11.1f KB\n", label, string(wet), ns,
                n * NSTEPS / el / 1e6, bytes / 1024)
    end
end
