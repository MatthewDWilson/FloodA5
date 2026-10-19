# test/check_determinism.jl
# -------------------------
# Thread-count determinism check (WP0, 2026-10).
#
# Runs a short fixed point-source simulation on an in-process Cartesian mesh
# for several solver configurations and prints one fingerprint line per
# configuration (a hash of the final volume and momentum arrays). Run it with
# different thread counts and compare the output: the lines must be
# IDENTICAL. Before the WP0 race fix, the edge-momentum (Q-centred)
# configurations were not deterministic across thread counts.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t 1 test/check_determinism.jl > det_t1.txt
#   FLOODA5_NO_VIS=1 julia --project=. -t 4 test/check_determinism.jl > det_t4.txt
#   diff det_t1.txt det_t4.txt        # (PowerShell: Compare-Object (gc det_t1.txt) (gc det_t4.txt))
#
# Only lines starting with "FINGERPRINT" are meant for comparison (log lines
# also go to stderr).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Logging

function _det_case(; nsteps::Int = 300, kwargs...)
    R    = A5Grid._EARTH_R
    dx   = 20.0
    dd   = rad2deg(dx / R)
    nrc  = 41
    mesh = generate_cartesian_mesh_bounds(0.0, 0.0, (nrc - 0.5) * dd, (nrc - 0.5) * dd,
                                          dx; lat_ref = 0.0)
    n = length(mesh.cells)
    # Gentle tilted bed with a small deterministic roughness pattern, so the
    # flow is genuinely 2-D and touches every code path.
    mesh.static_vars["elevation"] =
        [0.001 * (mesh.cells[i].center_lon / dd) * dx + 0.002 * sin(1.3 * i) for i in 1:n]
    state = initialise_flow_model(mesh, StandardFlow(); kwargs...)
    for k in eachindex(state.ghost_cell_bc)
        state.ghost_cell_bc[k] = Closed
    end
    src = div(n, 2) + 1                       # centre cell
    for _ in 1:nsteps
        dt = 0.5
        state.volume[src] += 2.0 * dt         # 2 m³/s
        for i in 1:n
            state.water_depth[i] = state.volume[i] / state.cell_area[i]
        end
        step_standard!(state, dt)
    end
    return hash((state.volume, state.edges.flux, state.qvec_u, state.qvec_v))
end

configs = [
    ("edge_theta0.9",        (;)),
    ("edge_theta1.0",        (; q_centre_theta = 1.0)),
    ("wlsq_edge",            (; gradient_correction = true)),
    ("diamond_cell",         (; gradient_correction = true, face_flux_method = :diamond,
                                momentum_model = :cell)),
]

with_logger(ConsoleLogger(stderr, Logging.Warn)) do
    for (label, kw) in configs
        h = _det_case(; kw...)
        println("FINGERPRINT $label $(string(h, base = 16))")
    end
end
println(stderr, "threads = $(Threads.nthreads())")
