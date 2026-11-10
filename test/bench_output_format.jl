# test/bench_output_format.jl — output format 1 (≤ 0.3.x, one group per frame)
# against format 2 (0.4.0, one dataset per variable): write time, file size and
# read times, using the frames of an existing output file.
#
#   julia --project=. test/bench_output_format.jl RUN.h5 [n_series_cells]
#
# The frames of RUN.h5 (either format) are written again in both layouts
# (both with deflate 4; format 1 with the 0.3.x writer's settings: per-frame groups, 4096-cell
# chunks, deflate 4), each file reopened per frame as the model does. Reads:
# one whole frame (the last), the time series of n cells (default 20, spread
# over the domain), and one variable for all frames.

isdefined(Main, :FloodA5Output) || include(joinpath(@__DIR__, "..", "io", "FloodA5Output.jl"))
const FO = FloodA5Output
using HDF5, Printf

src = ARGS[1]
nser = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20
const VARS = ("water_depth", "volume", "saturation", "velocity", "vel_u", "vel_v")

ct = FO.cell_table(src)
n = length(ct.ids)
ts = FO.frame_times(src)
nf = length(ts)
frames = [Dict(v => FO.read_frame(src, v, k) for v in VARS if FO.has_var(src, v)) for k in 1:nf]
vars = [v for v in VARS if haskey(frames[1], v)]
@printf("%s: %d cells, %d frames, %d variables\n", basename(src), n, nf, length(vars))

dir = mktempdir()
p1 = joinpath(dir, "v1.h5"); p2 = joinpath(dir, "v2.h5")

# ── format 1, as written by FloodA5 0.3.x ────────────────────────────────────
h5open(p1, "w") do f
    g = create_group(f, "mesh"); g["cell_ids"] = ct.ids
    create_group(f, "frames")
end
t1 = @elapsed for k in 1:nf
    h5open(p1, "r+") do f
        g = create_group(f["frames"], lpad(string(k), 6, '0'))
        g["t"] = ts[k]
        for v in vars
            ds = create_dataset(g, v, Float64, (n,); chunk = (min(n, 4096),), deflate = 4)
            write(ds, frames[k][v])
        end
    end
end

# ── format 2 ────────────────────────────────────────────────────────────────
cells = (ids = ct.ids, center_lons = ct.lon, center_lats = ct.lat, elevations = ct.elevation,
         static_vars = Dict{String,Vector{Float64}}(), id_scheme = "bench")
FO.create_output(p2; cells = cells)
t2 = @elapsed for k in 1:nf
    FO.append_frame!(p2, ts[k], [v => frames[k][v] for v in vars])
end

cellsel = round.(Int, range(1, n; length = nser))
function reads(p)
    h5open(p, "r") do f
        FO.read_frame(f, "water_depth", 1); FO.read_series(f, "water_depth", [1]); # warm-up
        a = @elapsed FO.read_frame(f, "water_depth", nf)
        b = @elapsed FO.read_series(f, "water_depth", cellsel)
        c = @elapsed FO.read_var(f, "water_depth")
        (a, b, c)
    end
end
r1 = reads(p1); r2 = reads(p2)
@assert FO.read_var(p1, "water_depth") == FO.read_var(p2, "water_depth")

println()
@printf("%-34s %12s %12s\n", "", "format 1", "format 2")
@printf("%-34s %12.2f %12.2f\n", "write, all frames (s)", t1, t2)
@printf("%-34s %12.1f %12.1f\n", "file size (MB)", filesize(p1) / 1e6, filesize(p2) / 1e6)
@printf("%-34s %12.4f %12.4f\n", "read one frame (s)", r1[1], r2[1])
@printf("%-34s %12.4f %12.4f\n", "read time series, cells (s)", r1[2], r2[2])
@printf("%-34s %12.4f %12.4f\n", "read one variable, all frames (s)", r1[3], r2[3])
println("(time series: $nser cells; identical data in both files)")
