# test/benchmark_backend.jl
# -------------------------
# WP5: steady-state step throughput of the SWE solver on each backend, after
# warm-up (compilation excluded), on square Cartesian grids of increasing
# size: a dam break over a sloping bed, all cells wet, order 2. Reports
# cell·sub-steps per second and the GPU/CPU speed-up, and checks that the
# device result matches the CPU one.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t auto test/benchmark_backend.jl [--sizes 128,256,512,1024] [--steps 20] [--gpu]

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
using Logging, Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const SIZES = parse.(Int, split(_arg("--sizes", "128,256,512"), ","))
const NSTEP = parse(Int, _arg("--steps", "20"))
const DX = 10.0
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))
const _R = A5Grid._EARTH_R

function state_for(n)
    d = rad2deg(DX / _R)
    mesh = _quiet(() -> generate_cartesian_mesh_bounds(0.0, 0.0, (n - 0.5) * d, (n - 0.5) * d, DX; lat_ref = 0.0))
    x = [c.center_lon for c in mesh.cells] ./ ((n - 0.5) * d)
    mesh.static_vars["elevation"] = 0.5 .* (1 .- x)
    s = _quiet(() -> initialise_flow_model(mesh, SWEFlow(); manning_n = 0.03))
    fill!(s.ghost_cell_bc, ZeroGradient)
    s.volume .= ([xi < 0.3 ? 2.0 : 1.0 for xi in x] .+ 0.5 .* (1 .- x)) .* s.cell_area .- s.elevation .* s.cell_area
    s.volume .= max.(s.volume, 0.0); s.water_depth .= s.volume ./ s.cell_area
    return s
end

# a fixed number of sub-steps at a fixed dt (the first CFL step)
function run_steps!(be, a, dt, nstep)
    removed = Ref(0.0)
    for _ in 1:nstep
        removed[] += _swe_substep!(be, a, dt)
    end
    synchronize!(be)
    return removed[]
end

gpu = "--gpu" in ARGS ? _cuda_backend() : nothing
@printf("%8s %10s %12s %12s %9s %12s\n", "grid", "cells", "CPU Mc·s/s", "GPU Mc·s/s", "speed-up", "max |ΔV|/V")
for n in SIZES
    s = state_for(n)
    a = _swe_args(s); be = CPUThreads()
    dt = _swe_cfl(be, a, 0.5)
    sg = deepcopy(s)
    run_steps!(be, a, dt, 2)                                       # warm-up
    tc = @elapsed run_steps!(be, a, dt, NSTEP)
    cps_cpu = length(s.volume) * NSTEP / tc / 1e6
    cps_gpu = NaN; err = NaN
    if gpu !== nothing
        ag = to_device(gpu, _swe_args(sg))
        Base.invokelatest(run_steps!, gpu, ag, dt, 2)
        tg = @elapsed Base.invokelatest(run_steps!, gpu, ag, dt, NSTEP)
        cps_gpu = length(s.volume) * NSTEP / tg / 1e6
        Vg = Array(ag.V)
        err = maximum(abs.(Vg .- s.volume) ./ max.(s.volume, 1.0))
    end
    @printf("%8s %10d %12.1f %12.1f %9.2f %12.2e\n", "$(n)×$(n)", length(s.volume), cps_cpu, cps_gpu,
            cps_gpu / cps_cpu, err)
end
