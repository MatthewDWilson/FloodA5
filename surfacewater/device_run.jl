# surfacewater/device_run.jl
# --------------------------
# WP5 (2026-10): the time loop's per-step work on a compute device
# (`run_simulation!(…; backend = KADevice(…))`, `--backend gpu`).
#
# The solver state lives on the device for the whole run: Δt, rainfall,
# point sources, the depth sync and the step itself are device kernels. The
# host keeps the FlowState for I/O and copies the state back only when it
# needs it (output frames, the visualiser, the log lines; `_device_sync!`).
# Point-source volumes are evaluated on the host each step (hydrograph
# interpolation) and added on the device in source order by a one-thread
# kernel, so cells fed by several sources get the same sum as on the CPU.
#
# Solvers: `--flow-model swe` (stage 1), `standard` (stage 2; the options
# listed in surfacewater/std_device.jl) and `sgs` (stage 3, local CFL).
# Not on a device: --limiter-stats, legacy CFL, cell/blend momentum.

struct DeviceRun{K, B<:ComputeBackend, AT, VI, VF}
    be         :: B
    a          :: AT                      # the solver's arrays on the device
    src_cell   :: VI                      # point sources: cell index (0 = nothing this step)
    src_inc    :: VF                      # point sources: volume added this step (m³)
    src_cell_h :: Vector{Int}             # host staging for the two above
    src_inc_h  :: Vector{Float64}
    removed    :: Base.RefValue{Float64}  # SWE: boundary outflow so far (state.vol_removed)
end

"""
    DeviceRun(state, method, backend, sources) → DeviceRun

Copy the solver's arrays to the device. `SWEFlow`, `StandardFlow` (edge
momentum, local CFL) and `SGSFlow` (local CFL) are supported.
"""
function DeviceRun(state::FlowState, method::FlowMethod, be::ComputeBackend, sources)
    for s in sources
        s isa Union{InflowPoint, InjectionPoint, RainPoint} ||
            throw(ArgumentError("source type $(typeof(s)) is not supported on a device backend"))
    end
    ns = length(sources)
    stage(K, a) = DeviceRun{K, typeof(be), typeof(a), typeof(to_device(be, zeros(Int, 0))),
                            typeof(to_device(be, zeros(0)))}(
        be, a, to_device(be, zeros(Int, ns)), to_device(be, zeros(ns)),
        zeros(Int, ns), zeros(ns), Ref(state.vol_removed))
    if method isa SWEFlow
        sw = state.swe
        @inbounds for k in eachindex(sw.gbc); sw.gbc[k] = _swe_bc_code(state.ghost_cell_bc[k]); end
        return stage(:swe, to_device(be, _swe_args(state)))
    elseif method isa StandardFlow
        state.momentum_model === :edge ||
            throw(ArgumentError("a device backend supports --momentum-model edge only"))
        state.cfl_mode === :local ||
            throw(ArgumentError("a device backend supports --cfl-mode local only"))
        return stage(:std, to_device(be, _std_device_args(state)))
    elseif method isa SGSFlow
        state.cfl_mode === :local ||
            throw(ArgumentError("a device backend supports --cfl-mode local only"))
        return stage(:sgs, to_device(be, _sgs_device_args(state)))
    end
    throw(ArgumentError("flow method $(typeof(method)) is not supported on a device backend"))
end

_dev_V(dev::DeviceRun{:swe}) = dev.a.V
_dev_V(dev::DeviceRun) = dev.a.volume
_dev_A(dev::DeviceRun{:swe}) = dev.a.A
_dev_A(dev::DeviceRun) = dev.a.cell_area

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
    V = _dev_V(dev)
    rainfall_rate > 0.0 && foreach_index!(dev.be, length(V), _rain_k!, V, _dev_A(dev), rainfall_rate * dt)
    isempty(sources) && return nothing
    for (k, s) in enumerate(sources)
        v = _src_volume(s, t, dt)
        dev.src_cell_h[k] = v === nothing ? 0 : s.cell_index
        dev.src_inc_h[k]  = v === nothing ? 0.0 : v
    end
    copyto!(dev.src_cell, dev.src_cell_h); copyto!(dev.src_inc, dev.src_inc_h)
    on_host_order!(dev.be, length(sources), _point_src_k!, V, dev.src_cell, dev.src_inc)
    return nothing
end

_device_cfl(dev::DeviceRun{:swe}, courant::Float64) = _swe_cfl(dev.be, dev.a, courant)
_device_cfl(dev::DeviceRun{:std}, courant::Float64) = _std_cfl(dev.be, dev.a, courant)
_device_cfl(dev::DeviceRun{:sgs}, courant::Float64) = _sgs_cfl(dev.be, dev.a, courant)
_device_depth!(dev::DeviceRun{:swe}) = foreach_index!(dev.be, dev.a.n, _swe_depth_k!, dev.a)
_device_depth!(dev::DeviceRun{:std}) = foreach_index!(dev.be, dev.a.n, _depth_std_k!, dev.a)
_device_depth!(dev::DeviceRun{:sgs}) = nothing        # step_sgs! keeps the depth in sync itself
_device_step!(dev::DeviceRun{:swe}, dt::Float64, courant::Float64) =
    _swe_advance!(dev.be, dev.a, dt, courant, dev.removed)
_device_step!(dev::DeviceRun{:std}, dt::Float64, courant::Float64) = _std_step!(dev.be, dev.a, dt)
_device_step!(dev::DeviceRun{:sgs}, dt::Float64, courant::Float64) = _sgs_step!(dev.be, dev.a, dt)

"Copy the state the host reads back (volume, depth, discharges, outflow)."
function _device_sync!(dev::DeviceRun{:swe}, state::FlowState)
    synchronize!(dev.be)
    a = dev.a
    copyto!(state.volume, a.V); copyto!(state.water_depth, a.d)
    copyto!(state.qvec_u, a.qx); copyto!(state.qvec_v, a.qy)
    copyto!(state.swe.clamp_vol, a.clamp)
    state.vol_removed = dev.removed[]
    return nothing
end
function _device_sync!(dev::DeviceRun{K}, state::FlowState) where {K}   # :std, :sgs
    synchronize!(dev.be)
    a = dev.a
    copyto!(state.volume, a.volume); copyto!(state.water_depth, a.water_depth)
    K === :sgs ? copyto!(state.edges.flux_Q, a.flux_Q) : copyto!(state.edges.flux, a.flux)
    copyto!(state.qvec_u, a.qvec_u); copyto!(state.qvec_v, a.qvec_v)
    state.vol_removed = Array(a.removed)[1]
    gf = Array(a.g_fprev)                     # ghost-edge momentum, for a later CPU step
    for (k, g) in enumerate(state.ghost_edges)
        g.flux_prev = gf[k]
        K === :std && (g.is_Q_flux = false)
    end
    return nothing
end
