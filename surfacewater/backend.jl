# surfacewater/backend.jl
# -----------------------
# WP5 (2026-10): compute backends — one code base for CPU threads and GPUs.
#
# Every per-step pass is a kernel body `body(i, args...)` for one cell, face
# or boundary face, writing only its own slots (no atomics, no shared writes;
# surfacewater/kernels.jl, swe_reference.jl). `foreach_index!` launches a
# body over 1:n on a backend:
#
#   CPUThreads()     Polyester `@batch per=core` — the default, and the
#                    reference path: identical arithmetic to the loops it
#                    replaced, so results are bit-identical to earlier versions
#   KADevice(ka)     a KernelAbstractions backend: `KADevice(CUDABackend())`
#                    runs on an NVIDIA GPU, `KADevice(KernelAbstractions.CPU())`
#                    runs the device code path on the CPU (tests without a GPU)
#
# The arrays a pass needs are gathered into a NamedTuple; `to_device` adapts
# it to the backend (CuArrays for CUDA; the same Arrays for the CPU backends).
# Reductions go through `_minimum` / `_dot_sum`, which keep the serial
# summation order on host arrays (so KADevice(CPU()) matches CPUThreads
# exactly) and use the device's own reductions on GPU arrays (round-off-level
# differences in sums; minima are exact either way).
#
# CUDA.jl is loaded only when a GPU run is requested (`_cuda_backend()`), so
# CPU runs and CI need neither CUDA nor a GPU.

using KernelAbstractions
import Adapt

abstract type ComputeBackend end

"CPU threads via Polyester `@batch per=core`: the default and reference path."
struct CPUThreads <: ComputeBackend end

"""
    KADevice(ka)

A KernelAbstractions backend: `KADevice(CUDABackend())` for an NVIDIA GPU;
`KADevice(KernelAbstractions.CPU())` runs the device code path on the CPU.
"""
struct KADevice{B} <: ComputeBackend
    ka :: B
end

@kernel function _foreach_kernel!(body, args)
    i = @index(Global, Linear)
    body(i, args...)
end

"""
    foreach_index!(backend, n, body, args...)

Run `body(i, args...)` for i = 1:n. The body must write only slots owned by
`i`, and must not read what another `i` of the same launch writes.
"""
@inline function foreach_index!(::CPUThreads, n::Int, body::F, args...) where {F}
    @batch per=core for i in 1:n
        body(i, args...)
    end
    return nothing
end
function foreach_index!(be::KADevice, n::Int, body::F, args...) where {F}
    n > 0 && _foreach_kernel!(be.ka, 256)(body, args; ndrange = n)
    return nothing
end

"""
    on_host_order!(backend, n, body, args...)

Run `body(i, args...)` for i = 1:n in order, on one thread: a plain loop on
the CPU, a single-work-item kernel on a device (for short loops whose order
matters, e.g. point sources added to cells in a fixed order).
"""
on_host_order!(::CPUThreads, n::Int, body::F, args...) where {F} =
    (for i in 1:n; body(i, args...); end; nothing)
function on_host_order!(be::KADevice, n::Int, body::F, args...) where {F}
    n > 0 && _serial_kernel!(be.ka, 1)(body, n, args; ndrange = 1)
    return nothing
end
@kernel function _serial_kernel!(body, n, args)
    for i in 1:n
        body(i, args...)
    end
end

synchronize!(::CPUThreads) = nothing
synchronize!(be::KADevice) = KernelAbstractions.synchronize(be.ka)

"Adapt a NamedTuple (or array) of host arrays to the backend's memory."
to_device(::CPUThreads, x) = x
to_device(be::KADevice, x) = Adapt.adapt(be.ka, x)

# Reductions: serial on host arrays (bit-identical to the loops they replace),
# the device's reductions otherwise.
function _minimum(x::Array{Float64})
    m = Inf
    @inbounds for v in x
        m = min(m, v)
    end
    return m
end
_minimum(x::AbstractArray) = minimum(x; init = Inf)
function _maximum(x::Array{Float64})
    m = 0.0
    @inbounds for v in x
        m = max(m, v)
    end
    return m
end
_maximum(x::AbstractArray) = maximum(x; init = 0.0)     # for non-negative x

"Σ_k dt·x[k]·y[k], in index order on host arrays."
function _dot_sum(dt::Float64, x::Array{Float64}, y::Array{Float64})
    s = 0.0
    @inbounds for k in eachindex(x)
        s += dt * x[k] * y[k]
    end
    return s
end
_dot_sum(dt::Float64, x::AbstractArray, y::AbstractArray) =
    dt * mapreduce(*, +, x, y; init = 0.0)

"""
    _cuda_backend() → KADevice

Load CUDA.jl (on first use) and return `KADevice(CUDABackend())`; errors if
no functional GPU is found. Call the GPU run through `Base.invokelatest`, so
that it sees the methods CUDA.jl has just defined.
"""
function _cuda_backend()
    if !isdefined(Main, :CUDA)
        try
            @eval Main using CUDA
        catch err
            error("--backend gpu needs CUDA.jl (in the full project environment): $err")
        end
    end
    ok = Base.invokelatest(() -> Main.CUDA.functional())
    ok || error("--backend gpu: CUDA.jl found no functional GPU")
    return KADevice(Base.invokelatest(() -> Main.CUDA.CUDABackend()))
end
