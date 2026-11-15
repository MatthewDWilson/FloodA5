# test/benchmark_backend.jl
# -------------------------
# WP5: steady-state step throughput on each backend and precision, after
# warm-up (compilation excluded), on square Cartesian grids of increasing
# size, all cells wet (a dam break over a sloping bed, open boundaries).
# Reports million cell-steps per second and, against the CPU-thread run of
# the same steps, the speed-up and the largest depth difference.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t auto test/benchmark_backend.jl
#       [--solver swe|standard] [--sizes 256,512,1024] [--steps 20]
#       [--min-cellsteps N] [--configs cpu,gpu64,gpu32] [--json FILE]
#
# --min-cellsteps N (default 0) raises the step count of a grid until at least
#   N cell-steps are timed: 20 steps of a small grid take a GPU a tenth of a
#   second, too short to time (512², standard: 25–38 Mc·s/s over 20 steps, 49
#   over 200). The benchmark suite uses 1e8.
# --configs (default cpu,gpu64,gpu32), in the order given:
#     cpu      CPU threads — the normal step (step_standard!, or step_swe!'s
#              sub-stepping on CPUThreads), Float64: the reference
#     gpu64    --backend gpu
#     gpu32    --backend gpu --precision float32 (0.6.0)
#     kacpu64, kacpu32   the device code path on the CPU (no GPU needed)
#   (--cpu / --gpu, from before 0.6.0: cpu only / gpu64,gpu32 only)
# --json FILE appends one JSON object per (size, config) to FILE
#   (benchmarks/BenchSuite.jl collects these; notebooks/05_performance.ipynb).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
using Logging, Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const SOLVER = Symbol(_arg("--solver", "swe"))
const SIZES  = parse.(Int, split(_arg("--sizes", "256,512"), ","))
const NSTEP  = parse(Int, _arg("--steps", "20"))
const MIN_CELLSTEPS = parse(Float64, _arg("--min-cellsteps", "0"))
nsteps(cells) = max(NSTEP, ceil(Int, MIN_CELLSTEPS / cells))
const CONFIGS = let c = _arg("--configs", "")
    c != "" ? String.(split(c, ",")) :
    "--cpu" in ARGS && !("--gpu" in ARGS) ? ["cpu"] :
    "--gpu" in ARGS && !("--cpu" in ARGS) ? ["gpu64", "gpu32"] : ["cpu", "gpu64", "gpu32"]
end
const JSON_OUT = _arg("--json", "")
const DX = 10.0
SOLVER in (:swe, :standard) || error("--solver must be swe or standard")
all(in(("cpu", "gpu64", "gpu32", "kacpu64", "kacpu32")), CONFIGS) ||
    error("--configs: cpu, gpu64, gpu32, kacpu64, kacpu32 (got $(join(CONFIGS, ",")))")
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
cpu_step!(s, ::SWEFlow, dt) = (r = Ref(s.vol_removed); _swe_advance!(CPUThreads(), _swe_args(s), dt, 0.5, r);
                               s.vol_removed = r[])
cpu_step!(s, ::StandardFlow, dt) = step_standard!(s, dt; velocity = false)
dev_step!(dev::DeviceRun, dt) = (_device_step!(dev, dt, 0.5); nothing)

# → (million cell-steps per second, seconds for set-up and the two warm-up steps)
function bench_cpu(s, method, dt, nstep)
    tw = @elapsed (for _ in 1:2; cpu_step!(s, method, dt); end)
    t = @elapsed (for _ in 1:nstep; cpu_step!(s, method, dt); end)
    return length(s.volume) * nstep / t / 1e6, tw
end
function bench_dev(be, s, method, dt, nstep)
    local dev
    tw = @elapsed begin
        dev = DeviceRun(s, method, be, AbstractSource[])
        for _ in 1:2; dev_step!(dev, dt); end; synchronize!(be)
    end
    t = @elapsed (for _ in 1:nstep; dev_step!(dev, dt); end; synchronize!(be))
    _device_sync!(dev, s)
    return length(s.volume) * nstep / t / 1e6, tw
end

function backend_for(c)
    c == "gpu64"   && return _cuda_backend()
    c == "gpu32"   && return _cuda_backend(; precision = Float32)
    c == "kacpu64" && return KADevice(KernelAbstractions.CPU())
    c == "kacpu32" && return KADevice(KernelAbstractions.CPU(); precision = Float32)
    return CPUThreads()
end
const BACKENDS = Dict(c => backend_for(c) for c in CONFIGS)
_json(d) = "{" * join(("\"$k\":" * (v isa AbstractString ? "\"$v\"" : v isa Real && !isfinite(v) ? "null" : string(v))
                       for (k, v) in d), ",") * "}"

# the standard solver's open boundaries as the command line sets them (--scheme recommended)
BOUNDARY_POSITIVITY[] = SOLVER === :standard

@printf("solver %s, %s steps after 2 warm-up steps, %d CPU threads\n", SOLVER,
        join(nsteps.(SIZES .^ 2), "/"), Threads.nthreads())
@printf("%10s %10s %8s %12s %9s %10s %11s %11s\n", "grid", "cells", "config", "Mc·s/s", "speed-up",
        "set-up s", "max |Δh| m", "mass err")
for n in SIZES
    s0, method = state_for(n)
    dt = 0.5 * (SOLVER === :swe ? _swe_cfl(CPUThreads(), _swe_args(s0), 1.0) : _cfl_dt(s0, method; courant = 1.0))
    v0 = sum(s0.volume); nstep = nsteps(length(s0.volume))
    ref = nothing; mcs_ref = NaN
    for c in CONFIGS
        s = deepcopy(s0)
        mcs, tw = c == "cpu" ? bench_cpu(s, method, dt, nstep) :
                  Base.invokelatest(bench_dev, BACKENDS[c], s, method, dt, nstep)
        mass = abs(sum(s.volume) + s.vol_removed - v0) / v0          # no sources: what is stored or left
        dh = ref === nothing ? NaN : maximum(abs.(s.water_depth .- ref.water_depth))
        rms = ref === nothing ? NaN : sqrt(sum(abs2, s.water_depth .- ref.water_depth) / length(s.volume))
        @printf("%10s %10d %8s %12.1f %9.2f %10.2f %11.2e %11.2e\n", "$(n)×$(n)", length(s.volume), c, mcs,
                mcs / mcs_ref, tw, dh, mass)
        if JSON_OUT != ""
            open(JSON_OUT, "a") do io
                println(io, _json(["solver" => string(SOLVER), "n" => n, "cells" => length(s.volume), "config" => c,
                                   "threads" => Threads.nthreads(), "steps" => nstep, "mcs" => mcs,
                                   "setup_s" => tw, "max_dh" => dh, "rms_dh" => rms, "mass_rel" => mass]))
            end
        end
        c == "cpu" && (ref = s; mcs_ref = mcs)
    end
end
