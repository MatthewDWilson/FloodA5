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
#
# 0.6.0: a KADevice also carries the working precision (`--precision`):
# Float64 (default), or Float32 for the fluxes and momentum with the volume
# kept in Float64 — surfacewater/precision.jl.

using KernelAbstractions
import Adapt

abstract type ComputeBackend end

"CPU threads via Polyester `@batch per=core`: the default and reference path."
struct CPUThreads <: ComputeBackend end

"""
    KADevice(ka; precision = Float64, levels = <by precision>)

A KernelAbstractions backend: `KADevice(CUDABackend())` for an NVIDIA GPU;
`KADevice(KernelAbstractions.CPU())` runs the device code path on the CPU.
`precision` is the working float type of the solver arrays (Float64, or
Float32 with the volume kept in Float64); `levels` the water-level
formulation, `AbsLevels()` for Float64 and `RelLevels()` for Float32
(`RelLevels()` with Float64 is for tests) — surfacewater/precision.jl.
"""
struct KADevice{B, T<:AbstractFloat, LV} <: ComputeBackend
    ka  :: B
    lev :: LV
end
function KADevice(ka; precision::Type = Float64,
                  levels = precision === Float64 ? AbsLevels() : RelLevels())
    precision in (Float64, Float32) ||
        throw(ArgumentError("precision must be Float64 or Float32 (got $precision)"))
    levels isa Union{AbsLevels, RelLevels} ||
        throw(ArgumentError("levels must be AbsLevels() or RelLevels()"))
    precision === Float32 && levels isa AbsLevels &&
        throw(ArgumentError("Float32 needs RelLevels(): absolute water levels are not resolved in Float32"))
    return KADevice{typeof(ka), precision, typeof(levels)}(ka, levels)
end

"Working float type of a backend's solver arrays."
_precision(::ComputeBackend) = Float64
_precision(::KADevice{B, T}) where {B, T} = T
"Water-level formulation of a backend (surfacewater/precision.jl)."
_levels(::ComputeBackend) = AbsLevels()
_levels(be::KADevice) = be.lev

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
_minimum(x::AbstractArray) = minimum(x; init = typemax(eltype(x)))
function _maximum(x::Array{Float64})
    m = 0.0
    @inbounds for v in x
        m = max(m, v)
    end
    return m
end
_maximum(x::AbstractArray) = maximum(x; init = zero(eltype(x)))     # for non-negative x

"Σ_k dt·x[k]·y[k], in index order on host arrays."
function _dot_sum(dt::Float64, x::Array{Float64}, y::Array{Float64})
    s = 0.0
    @inbounds for k in eachindex(x)
        s += dt * x[k] * y[k]
    end
    return s
end
# (each product in the arrays' type, as the cell update forms it; the sum in Float64)
_dot_sum(dt::Float64, x::AbstractArray, y::AbstractArray) =
    dt * mapreduce((a, b) -> Float64(a * b), +, x, y; init = 0.0)

"""
    _cuda_backend(; precision = Float64) → KADevice

Load CUDA.jl (on first use) and return `KADevice(CUDABackend(); precision)`;
errors if no functional GPU is found. Call the GPU run through
`Base.invokelatest`, so that it sees the methods CUDA.jl has just defined.
"""
function _cuda_backend(; precision::Type = Float64)
    if !isdefined(Main, :CUDA)
        try
            @eval Main using CUDA
        catch err
            error("--backend gpu needs CUDA.jl (in the full project environment): $err")
        end
    end
    ok = Base.invokelatest(() -> Main.CUDA.functional())
    ok || error("--backend gpu: CUDA.jl found no functional GPU")
    return KADevice(Base.invokelatest(() -> Main.CUDA.CUDABackend()); precision = precision)
end
