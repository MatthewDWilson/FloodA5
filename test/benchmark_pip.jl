"""
benchmark_pip.jl
----------------
Benchmarks the point-in-polygon (PIP) sampling step of the A5 flood model
across three backends:

  1. GPU    — CUDA kernel  (if CUDA.jl is available and GPU is functional)
  2. CPU-∥  — Threads.@threads  (Julia multi-threaded, --threads auto)
  3. CPU-1  — single-threaded loop (baseline)

Also benchmarks the full mesh_for_aoi pipeline (PIP + Python bridge) so you
can see how much of the wall time is PIP vs Python overhead.

Usage
-----
    julia --threads auto benchmark_pip.jl

Output
------
Prints a Markdown table of median times ± IQR for each (backend, point-count)
combination, then a summary table for the full pipeline at each resolution.

Requires
--------
    - CUDA.jl  (optional; GPU benchmarks skipped if not available)
    - BenchmarkTools.jl  (Pkg.add("BenchmarkTools"))
    - A5Grid.jl and a5_bridge.py live in mesh/ relative to the project root.

Notes
-----
The GPU benchmark measures TWO scenarios:
  (a) transfer-included: full round-trip Host → Device → Kernel → Host
  (b) kernel-only:       kernel launch + sync, data already on device
This lets you see whether the bottleneck is PCIe bandwidth or compute.
"""

using BenchmarkTools
using Printf
using Dates
using Statistics

# CUDA is loaded at top level so CuArray etc. are available in benchmark closures.
# The actual GPU benchmarks are skipped at runtime if no GPU is present.
const _HAVE_CUDA = try
    using CUDA
    true
catch
    false
end

# ── Load A5Grid (mesh/ relative to project root) ─────────────────────────
const SCRIPT_DIR  = @__DIR__
const PROJECT_DIR = dirname(SCRIPT_DIR)
push!(LOAD_PATH, joinpath(PROJECT_DIR, "mesh"))
include(joinpath(PROJECT_DIR, "mesh", "A5Grid.jl"))
using .A5Grid

# ── Geometry helpers (duplicated here so benchmark is self-contained) ──────

@inline function _pip_cpu_single(lon::Float64, lat::Float64,
                                  ring_xs::Vector{Float64},
                                  ring_ys::Vector{Float64})::Bool
    inside = false
    n = length(ring_xs)
    j = n
    @inbounds for i in 1:n
        xi, yi = ring_xs[i], ring_ys[i]
        xj, yj = ring_xs[j], ring_ys[j]
        if ((yi > lat) != (yj > lat)) &&
           (lon < (xj - xi) * (lat - yi) / (yj - yi + 1e-15) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function pip_serial(lons::Vector{Float64}, lats::Vector{Float64},
                    ring_xs::Vector{Float64}, ring_ys::Vector{Float64})::Vector{Bool}
    mask = Vector{Bool}(undef, length(lons))
    @inbounds for i in eachindex(lons)
        mask[i] = _pip_cpu_single(lons[i], lats[i], ring_xs, ring_ys)
    end
    return mask
end

function pip_parallel(lons::Vector{Float64}, lats::Vector{Float64},
                      ring_xs::Vector{Float64}, ring_ys::Vector{Float64})::Vector{Bool}
    mask = Vector{Bool}(undef, length(lons))
    Threads.@threads for i in eachindex(lons)
        mask[i] = _pip_cpu_single(lons[i], lats[i], ring_xs, ring_ys)
    end
    return mask
end

# ── GPU kernel — defined at TOP LEVEL so @benchmark can see it ──────────────
# BenchmarkTools generates a separate sample function; closures defined inside
# @benchmark blocks are not in scope there.  Top-level functions are.
function gpu_pip_kernel(lons, lats, rxs, rys, mask, nv)
    idx = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    idx > length(lons) && return nothing
    inside = false; j = nv
    lon = lons[idx]; lat = lats[idx]
    @inbounds for i in 1:nv
        xi, yi = rxs[i], rys[i]; xj, yj = rxs[j], rys[j]
        if ((yi > lat) != (yj > lat)) &&
           (lon < (xj - xi) * (lat - yi) / (yj - yi + 1e-15) + xi)
            inside = !inside
        end
        j = i
    end
    mask[idx] = inside
    return nothing
end

# ── Sample data: circular AOI centred on Christchurch ─────────────────────

function make_circle_ring(cx::Float64, cy::Float64, r::Float64, n::Int=64)
    xs = [cx + r * cos(2π*i/n) for i in 0:n-1]
    ys = [cy + r * sin(2π*i/n) for i in 0:n-1]
    push!(xs, xs[1]); push!(ys, ys[1])   # close ring
    return xs, ys
end

function make_sample_grid(n_side::Int)
    # Dense grid over Christchurch bounding box
    lons = Float64[]; lats = Float64[]
    for i in 1:n_side, j in 1:n_side
        push!(lons, 171.5 + (i-1) * (2.0 / n_side))
        push!(lats, -44.0 + (j-1) * (2.0 / n_side))
    end
    return lons, lats
end

# ── Benchmark runner ───────────────────────────────────────────────────────

struct BenchResult
    label::String
    n_points::Int
    median_ms::Float64
    iqr_ms::Float64
    throughput_Mpts_s::Float64
end

function run_pip_benchmarks()
    println("\n", "="^72)
    println("  PIP Benchmark — $(Dates.now())")
    println("  Julia threads: $(Threads.nthreads())")
    println("  CUDA available: $(A5Grid._cuda_available[])")
    println("="^72, "\n")

    # Ring: circle ~0.5° radius centred on Christchurch
    ring_xs, ring_ys = make_circle_ring(172.6, -43.5, 0.5, 128)

    # Point counts to test
    n_sides   = [64, 128, 256, 384, 512]   # n_side² points
    n_points_list = n_sides .^ 2            # 4K → 262K

    results = BenchResult[]

    for n_side in n_sides
        lons, lats = make_sample_grid(n_side)
        N = length(lons)
        @printf "── %6d points ─────────────────────────────────────────\n" N

        # ── CPU serial ──
        b = @benchmark pip_serial($lons, $lats, $ring_xs, $ring_ys) samples=20 evals=3
        med = median(b).time / 1e6      # ns → ms
        iqr_ = (quantile(b.times, 0.75) - quantile(b.times, 0.25)) / 1e6
        tp  = N / (median(b).time / 1e9) / 1e6
        push!(results, BenchResult("CPU-serial", N, med, iqr_, tp))
        @printf "  CPU-serial  : %7.2f ms  ±%-6.2f ms  (%6.1f M pts/s)\n" med iqr_ tp

        # ── CPU parallel ──
        b = @benchmark pip_parallel($lons, $lats, $ring_xs, $ring_ys) samples=20 evals=3
        med = median(b).time / 1e6
        iqr_ = (quantile(b.times, 0.75) - quantile(b.times, 0.25)) / 1e6
        tp  = N / (median(b).time / 1e9) / 1e6
        push!(results, BenchResult("CPU-parallel($(Threads.nthreads())t)", N, med, iqr_, tp))
        @printf "  CPU-parallel: %7.2f ms  ±%-6.2f ms  (%6.1f M pts/s)\n" med iqr_ tp

        # ── GPU (transfer-included) ──
        if A5Grid._cuda_available[]
            # Kernels referenced by name — defined at module scope above.
            b = @benchmark begin
                _dl = CUDA.CuArray($lons)
                _da = CUDA.CuArray($lats)
                _dx = CUDA.CuArray($ring_xs)
                _dy = CUDA.CuArray($ring_ys)
                _dm = CUDA.zeros(Bool, length($lons))
                _nv = length($ring_xs)
                CUDA.@cuda threads=256 blocks=cld(length($lons),256) gpu_pip_kernel(
                    _dl,_da,_dx,_dy,_dm,_nv)
                CUDA.synchronize()
                Array(_dm)
            end samples=20 evals=3
            med = median(b).time / 1e6
            iqr_ = (quantile(b.times, 0.75) - quantile(b.times, 0.25)) / 1e6
            tp  = N / (median(b).time / 1e9) / 1e6
            push!(results, BenchResult("GPU (incl. transfer)", N, med, iqr_, tp))
            @printf "  GPU (full)  : %7.2f ms  ±%-6.2f ms  (%6.1f M pts/s)\n" med iqr_ tp

            # ── GPU kernel-only (data already on device) ──
            d_lons = CUDA.CuArray(lons)
            d_lats = CUDA.CuArray(lats)
            d_rxs  = CUDA.CuArray(ring_xs)
            d_rys  = CUDA.CuArray(ring_ys)
            d_mask = CUDA.zeros(Bool, N)
            nv     = length(ring_xs)
            b = @benchmark begin
                CUDA.@cuda threads=256 blocks=cld($N,256) gpu_pip_kernel(
                    $d_lons,$d_lats,$d_rxs,$d_rys,$d_mask,$nv)
                CUDA.synchronize()
            end samples=50 evals=5
            med = median(b).time / 1e6
            iqr_ = (quantile(b.times, 0.75) - quantile(b.times, 0.25)) / 1e6
            tp  = N / (median(b).time / 1e9) / 1e6
            push!(results, BenchResult("GPU (kernel-only)", N, med, iqr_, tp))
            @printf "  GPU (kernel): %7.2f ms  ±%-6.2f ms  (%6.1f M pts/s)\n" med iqr_ tp
        end
        println()
    end

    return results
end

function run_pipeline_benchmarks()
    println("\n", "="^72)
    println("  Full Pipeline Benchmark (PIP + Python bridge)")
    println("  Resolutions 10, 12, 14, 16  ×  christchurch_aoi.geojson")
    println("="^72, "\n")

    aoi_path = joinpath(SCRIPT_DIR, "christchurch_aoi.geojson")
    isfile(aoi_path) || (println("  ⚠  christchurch_aoi.geojson not found — skipping"); return)

    results = []
    for res in [10, 12, 14, 16]
        @printf "── Resolution %d ───────────────────────────────────────────\n" res
        t0 = time_ns()
        mesh = A5Grid.mesh_for_aoi(aoi_path, res)
        elapsed_ms = (time_ns() - t0) / 1e6
        backend = A5Grid._cuda_available[] ? "GPU" : "CPU-∥"
        @printf "  Resolution %2d: %5d cells  %7.0f ms  [%s]\n" res length(mesh) elapsed_ms backend
        push!(results, (resolution=res, cells=length(mesh), ms=elapsed_ms, backend=backend))
    end
    return results
end

function print_summary(pip_results)
    println("\n", "="^72)
    println("  SUMMARY TABLE")
    println("="^72)
    @printf "\n  %-28s %8s %10s %10s %10s\n" "Backend" "N points" "Median ms" "± IQR ms" "M pts/s"
    println("  ", "-"^68)
    prev_n = 0
    for r in pip_results
        r.n_points != prev_n && (println(); prev_n = r.n_points)
        @printf "  %-28s %8d %10.2f %10.2f %10.1f\n" r.label r.n_points r.median_ms r.iqr_ms r.throughput_Mpts_s
    end
    println()

    # Cross-backend speedup at largest N
    backends = ["CPU-serial", "CPU-parallel", "GPU (incl. transfer)", "GPU (kernel-only)"]
    max_n = maximum(r.n_points for r in pip_results)
    at_max = filter(r -> r.n_points == max_n, pip_results)
    if length(at_max) >= 2
        serial_ms = first(r.median_ms for r in at_max if startswith(r.label, "CPU-serial"))
        println("  Speedup at N=$max_n (vs CPU-serial):")
        for r in at_max
            r.label == "CPU-serial" && continue
            @printf "    %-30s  %.1f×\n" r.label (serial_ms / r.median_ms)
        end
    end
    println()
end

# ── Main ──────────────────────────────────────────────────────────────────

pip_results      = run_pip_benchmarks()
pipeline_results = run_pipeline_benchmarks()
print_summary(pip_results)

println("""
Interpretation guide
────────────────────
• GPU (incl. transfer) ≈ GPU (kernel-only):
    Compute dominates — GPU is genuinely faster. Consider pre-uploading
    the sample grid once and reusing across resolutions.

• GPU (incl. transfer) >> GPU (kernel-only):
    PCIe bandwidth is the bottleneck (common for small N).
    Crossover point (where GPU beats CPU-∥) is typically 500K–2M points
    for a PIP kernel. Below that, multi-threaded CPU wins on total latency.

• CPU-∥ speedup ≈ thread count:
    Good scaling — no false-sharing or contention issues.

• CPU-∥ speedup << thread count:
    Memory-bandwidth limited. Check NUMA topology or reduce working set.

Recommendation for resolution 14 (~147K points)
────────────────────────────────────────────────
If GPU (incl. transfer) > CPU-∥ at 147K points, switch the default
backend to CPU for resolutions ≤ 15, GPU for ≥ 16. Set in A5Grid.jl:

    _cuda_threshold_pts = 500_000   # tune from benchmark results
    use_gpu = _cuda_available[] && length(lons) > _cuda_threshold_pts
""")
