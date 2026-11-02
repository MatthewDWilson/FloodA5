# test/profile_backend.jl
# -----------------------
# WP5: time per phase of the standard step (0.3.0 default options) on a
# backend, synchronising after each phase, on an n×n Cartesian grid with all
# cells wet. For finding where a device step spends its time; the totals are
# slower than an unsynchronised step.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t auto test/profile_backend.jl [--n 360] [--steps 50] [--gpu]

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
using Logging, Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const N = parse(Int, _arg("--n", "360")); const NSTEP = parse(Int, _arg("--steps", "50"))
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))
const DX = 5.5
d = rad2deg(DX / A5Grid._EARTH_R)
mesh = _quiet(() -> generate_cartesian_mesh_bounds(0.0, 0.0, (N - 0.5) * d, (N - 0.5) * d, DX; lat_ref = 0.0))
x = [c.center_lon for c in mesh.cells] ./ ((N - 0.5) * d)
mesh.static_vars["elevation"] = 0.5 .* (1 .- x)
s = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); scheme_options(:recommended, :standard)...))
fill!(s.ghost_cell_bc, ZeroGradient)
s.volume .= (1.0 .- s.elevation .+ 0.2) .* s.cell_area; s.water_depth .= s.volume ./ s.cell_area
BOUNDARY_POSITIVITY[] = true

be = "--gpu" in ARGS ? _cuda_backend() : KADevice(KernelAbstractions.CPU())
function profile(be, s)
    dev = DeviceRun(s, StandardFlow(), be, AbstractSource[])
    a = dev.a; dt = 0.5
    phases = [
        ("cfl (incl. host sync)", () -> _std_cfl(be, a, 0.7)),
        ("rain",            () -> _device_sources!(dev, 1e-6, AbstractSource[], 0.0, dt)),
        ("depth sync",      () -> _device_depth!(dev)),
        ("wse + gradients", () -> (foreach_index!(be, a.n, _wse_k!, a); foreach_index!(be, a.n, _grad_wse_k!, a))),
        ("diamond vertices", () -> a.use_diamond && (foreach_index!(be, a.nv, _veta_k!, a); a.wetdry && foreach_index!(be, a.nv, _vwet_k!, a))),
        ("flux_prev copy",  () -> copyto!(a.flux_prev, a.flux)),
        ("Phase A faces",   () -> foreach_index!(be, a.ne, _face_flux_std!, a.opts, a, dt)),
        ("positivity",      () -> (foreach_index!(be, a.n, _pos_out_k!, a); foreach_index!(be, a.ne, _pos_scale_k!, a))),
        ("gather",          () -> foreach_index!(be, a.n, _gather_std_k!, a)),
        ("Phase D ghosts ($(length(dev.glptr) - 1) levels)", () -> _ghost_pass!(be, _ghost_std_k!, a, dev.glptr, dt)),
        ("η filter",        () -> begin
            foreach_index!(be, a.n, _ef_eta_k!, a); foreach_index!(be, a.n, _ef_r_k!, a)
            foreach_index!(be, a.ne, _ef_flux_k!, a, a.eta_cs, dt); foreach_index!(be, a.n, _ef_out_k!, a)
            foreach_index!(be, a.ne, _ef_scale_k!, a); foreach_index!(be, a.n, _ef_apply_k!, a) end),
        ("Phase F qvec",    () -> foreach_index!(be, a.n, _qvec_k!, a, Val(false))),
    ]
    for (_, f) in phases; f(); end; _std_step!(be, a, dt, dev.glptr); synchronize!(be)   # compile
    tot = 0.0
    @printf("%-28s %10s\n", "phase", "ms/step")
    for (name, f) in phases
        synchronize!(be); t = time_ns()
        for _ in 1:NSTEP; f(); end
        synchronize!(be); ms = (time_ns() - t) / 1e6 / NSTEP; tot += ms
        @printf("%-28s %10.3f\n", name, ms)
    end
    @printf("%-28s %10.3f   (%d cells, %d faces, %d ghost faces)\n", "total", tot, a.n, a.ne, a.ng)
    synchronize!(be); t = time_ns()
    for _ in 1:NSTEP; _std_step!(be, a, dt, dev.glptr); end
    synchronize!(be)
    @printf("%-28s %10.3f\n", "_std_step! unsynchronised", (time_ns() - t) / 1e6 / NSTEP)
end
Base.invokelatest(profile, be, s)
