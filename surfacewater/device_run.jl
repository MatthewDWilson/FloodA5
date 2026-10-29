# surfacewater/device_run.jl
# --------------------------
# WP5 (2026-10): the time loop's per-step work on a compute device
# (`run_simulation!(…; backend = KADevice(…))`, `--backend gpu`).
#
# The solver state lives on the device for the whole run: Δt, rainfall,
# point sources, the depth sync and the step itself are device kernels. The
# host keeps the FlowState for I/O and copies the volume, depth and
# discharge vectors back only when it needs them (output frames, the
# visualiser, the log lines; `_device_sync!`). Point-source volumes are
# evaluated on the host each step (hydrograph interpolation) and added on the
# device in source order by a one-thread kernel, so cells fed by several
# sources get the same sum as on the CPU.
#
# Stage 1: `--flow-model swe`.

struct DeviceRun{B<:ComputeBackend, AT, VI, VF}
    be         :: B
    a          :: AT                      # the solver's arrays on the device
    src_cell   :: VI                      # point sources: cell index (0 = nothing this step)
    src_inc    :: VF                      # point sources: volume added this step (m³)
    src_cell_h :: Vector{Int}             # host staging for the two above
    src_inc_h  :: Vector{Float64}
    removed    :: Base.RefValue{Float64}  # boundary outflow so far (state.vol_removed)
end

"""
    DeviceRun(state, method, backend, sources) → DeviceRun

Copy the solver's arrays to the device. Supports `SWEFlow` (WP5 stage 1).
"""
function DeviceRun(state::FlowState, method::FlowMethod, be::ComputeBackend, sources)
    method isa SWEFlow ||
        throw(ArgumentError("the GPU backend supports --flow-model swe so far (got $(typeof(method)))"))
    for s in sources
        s isa Union{InflowPoint, InjectionPoint, RainPoint} ||
            throw(ArgumentError("source type $(typeof(s)) is not supported on a device backend"))
    end
    sw = state.swe
    @inbounds for k in eachindex(sw.gbc); sw.gbc[k] = _swe_bc_code(state.ghost_cell_bc[k]); end
    a = to_device(be, _swe_args(state))
    ns = length(sources)
    return DeviceRun(be, a, to_device(be, zeros(Int, ns)), to_device(be, zeros(ns)),
                     zeros(Int, ns), zeros(ns), Ref(state.vol_removed))
end

@inline function _rain_k!(i, V, A, r)
    @inbounds V[i] += r * A[i]
    return nothing
end

@inline function _point_src_k!(k, V, cell, inc)
    @inbounds begin
        c = cell[k]
        c > 0 && (V[c] += inc[k])
    end
    return nothing
end

# the volume a point source adds over [t, t + dt], or nothing (as apply_source!)
_src_volume(s::InflowPoint, t, dt) = (Q = _interp_hydrograph(s.t_s, s.Q_m3s, t); Q > 0.0 ? Q * dt : nothing)
_src_volume(s::Union{InjectionPoint, RainPoint}, t, dt) = s.rate_m3s * dt

"Uniform rainfall and point sources over [t, t + dt], as in run_simulation!."
function _device_sources!(dev::DeviceRun, rainfall_rate::Float64, sources, t::Float64, dt::Float64)
    a = dev.a
    rainfall_rate > 0.0 && foreach_index!(dev.be, a.n, _rain_k!, a.V, a.A, rainfall_rate * dt)
    isempty(sources) && return nothing
    for (k, s) in enumerate(sources)
        v = _src_volume(s, t, dt)
        dev.src_cell_h[k] = v === nothing ? 0 : s.cell_index
        dev.src_inc_h[k]  = v === nothing ? 0.0 : v
    end
    copyto!(dev.src_cell, dev.src_cell_h); copyto!(dev.src_inc, dev.src_inc_h)
    on_host_order!(dev.be, length(sources), _point_src_k!, a.V, dev.src_cell, dev.src_inc)
    return nothing
end

_device_cfl(dev::DeviceRun, courant::Float64) = _swe_cfl(dev.be, dev.a, courant)
_device_depth!(dev::DeviceRun) = foreach_index!(dev.be, dev.a.n, _swe_depth_k!, dev.a)
_device_step!(dev::DeviceRun, dt::Float64, courant::Float64) =
    _swe_advance!(dev.be, dev.a, dt, courant, dev.removed)

"Copy the state the host reads (volume, depth, discharge vectors, outflow) back."
function _device_sync!(dev::DeviceRun, state::FlowState)
    synchronize!(dev.be)
    a = dev.a
    copyto!(state.volume, a.V); copyto!(state.water_depth, a.d)
    copyto!(state.qvec_u, a.qx); copyto!(state.qvec_v, a.qy)
    copyto!(state.swe.clamp_vol, a.clamp)
    state.vol_removed = dev.removed[]
    return nothing
end
