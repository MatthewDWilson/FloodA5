# test/benchmark_backend.jl
# -------------------------
# WP5: steady-state step throughput on each backend, after warm-up
# (compilation excluded), on square Cartesian grids of increasing size, all
# cells wet (a dam break over a sloping bed, open boundaries). Reports cell
# steps per second, and with both backends the speed-up and the largest
# relative volume difference after the same steps.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t auto test/benchmark_backend.jl
#       [--solver swe|standard] [--sizes 256,512,1024] [--steps 20] [--cpu] [--gpu]
#
# --cpu / --gpu choose the backends (default: both); the CPU path is the
# normal step (step_standard!, or step_swe!'s sub-stepping on CPUThreads).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
using Logging, Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const SOLVER = Symbol(_arg("--solver", "swe"))
const SIZES  = parse.(Int, split(_arg("--sizes", "256,512"), ","))
const NSTEP  = parse(Int, _arg("--steps", "20"))
const DO_CPU = "--cpu" in ARGS || !("--gpu" in ARGS)
const DO_GPU = "--gpu" in ARGS || !("--cpu" in ARGS)
const DX = 10.0
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))

function state_for(n)
    d = rad2deg(DX / A5Grid._EARTH_R)
    mesh = _quiet(() -> generate_cartesian_mesh_bounds(0.0, 0.0, (n - 0.5) * d, (n - 0.5) * d, DX; lat_ref = 0.0))
    x = [c.center_lon for c in mesh.cells] ./ ((n - 0.5) * d)
    mesh.static_vars["elevation"] = 0.5 .* (1 .- x)
    method = SOLVER === :swe ? SWEFlow() : StandardFlow()
    kw = SOLVER === :swe ? (;) : scheme_options(:recommended, :standard)
    s = _quiet(() -> initialise_flow_model(mesh, method; manning_n = 0.03, kw...))
    fill!(s.ghost_cell_bc, ZeroGradient)
    s.volume .= max.(([xi < 0.3 ? 2.0 : 1.0 for xi in x] .+ 0.5 .* (1 .- x) .- s.elevation), 0.0) .* s.cell_area
    s.water_depth .= s.volume ./ s.cell_area
    return s, method
end

# the CPU step and the device step, at a fixed dt (SWE: both sub-step to dt
# at Courant 0.5, with the CFL reduction and depth sync, as step_swe!)
cpu_step!(s, ::SWEFlow, dt) = _swe_advance!(CPUThreads(), _swe_args(s), dt, 0.5, Ref(0.0))
cpu_step!(s, ::StandardFlow, dt) = step_standard!(s, dt; velocity = false)
dev_step!(dev::DeviceRun, dt) = (_device_step!(dev, dt, 0.5); nothing)

function bench_cpu(s, method, dt)
    for _ in 1:2; cpu_step!(s, method, dt); end
    t = @elapsed (for _ in 1:NSTEP; cpu_step!(s, method, dt); end)
    return length(s.volume) * NSTEP / t / 1e6
end
function bench_gpu(be, s, method, dt)
    dev = DeviceRun(s, method, be, AbstractSource[])
    for _ in 1:2; dev_step!(dev, dt); end; synchronize!(be)
    t = @elapsed (for _ in 1:NSTEP; dev_step!(dev, dt); end; synchronize!(be))
    _device_sync!(dev, s)
    return length(s.volume) * NSTEP / t / 1e6
end

gpu = DO_GPU ? _cuda_backend() : nothing
@printf("solver %s, %d steps, %d CPU threads\n", SOLVER, NSTEP, Threads.nthreads())
@printf("%10s %10s %12s %12s %9s %12s\n", "grid", "cells", "CPU Mc·s/s", "GPU Mc·s/s", "speed-up", "max |ΔV|/V")
for n in SIZES
    s, method = state_for(n)
    dt = 0.5 * (SOLVER === :swe ? _swe_cfl(CPUThreads(), _swe_args(s), 1.0) : _cfl_dt(s, method; courant = 1.0))
    sg = deepcopy(s)
    cps_cpu = DO_CPU ? bench_cpu(s, method, dt) : NaN
    cps_gpu = gpu === nothing ? NaN : Base.invokelatest(bench_gpu, gpu, sg, method, dt)
    err = DO_CPU && gpu !== nothing ? maximum(abs.(sg.volume .- s.volume) ./ max.(s.volume, 1.0)) : NaN
    @printf("%10s %10d %12.1f %12.1f %9.2f %12.2e\n", "$(n)×$(n)", length(s.volume), cps_cpu, cps_gpu,
            cps_gpu / cps_cpu, err)
end
