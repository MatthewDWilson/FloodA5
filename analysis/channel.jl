# analysis/channel.jl — straight-channel representation series (HYDRAULICS.md §15.12).
#
# The series (test/channel_bench/make_channel_meshes.jl): a straight
# rectangular channel (width W, bank depth D, bed slope S, floodplain cross
# slope SX, bed Z0 at the inflow) at angle α to east and offset `off` across
# the cells, inflow Q, Manning n. Metrics from a run's last output frame over
# the measurement reach s ∈ [−0.3 L, 0.05 L]:
#   deta     mean water-surface error (m) in the cells that overlap the
#            channel, against the 2-D normal-flow surface Z0 − S·s + h_n.
#            Lost conveyance shows as a positive error.
#   spill    wet (> 5 cm) floodplain area per metre of reach (m) in cells
#            wholly outside the channel; 0 for in-bank flow
#   unsteady max |Δη| over the channel cells between the last two frames (m)
#   outflow  discharge leaving the domain over the last frame interval,
#            as a fraction of the inflow (1 at steady state), from the run log

export CHANNEL, channel_normal_depth, channel_coords, channel_masks, channel_metrics,
       parse_log_volumes, log_volume_at, log_outflow, read_channel_manifest,
       channel_run_metrics, channel_series_metrics

"""
    CHANNEL

The channel series' constants (make_channel_meshes.jl): length `L`, width
`W`, bank depth `D` (m), bed slope `S`, floodplain cross slope `SX`, bed `Z0`
(m) at the inflow, inflow `Q` (m³/s), Manning `NMAN`, origin (`LON0`, `LAT0`).
Every channel function takes these as keyword overrides.
"""
const CHANNEL = (L = 2000.0, W = 20.0, D = 2.0, S = 0.001, SX = 0.005, Z0 = 10.0,
                 Q = 20.0, NMAN = 0.03, LON0 = 0.0, LAT0 = 51.0)

"Normal depth (m) of a wide rectangular channel: (Q n / (W √S))^(3/5)."
channel_normal_depth(; Q = CHANNEL.Q, NMAN = CHANNEL.NMAN, W = CHANNEL.W, S = CHANNEL.S) =
    (Q * NMAN / (W * sqrt(S)))^(3 / 5)

"""
    channel_coords(lon, lat; angle, offset, LON0, LAT0) -> (s, n)
    channel_coords(mesh; angle, offset, ...)

Along-channel `s` and across-channel `n` coordinates (m) of points (cell
centres): local x east, y north about (LON0, LAT0) (`local_xy`), rotated by
`angle` (degrees from east) and shifted across by `offset` (m).
"""
function channel_coords(lon::AbstractVector, lat::AbstractVector; angle::Real, offset::Real,
                        LON0 = CHANNEL.LON0, LAT0 = CHANNEL.LAT0, kw...)
    x, y = local_xy(lon, lat, LON0, LAT0)
    α = angle
    s = [x[i] * cosd(α) + y[i] * sind(α) for i in eachindex(x)]
    n = [-x[i] * sind(α) + y[i] * cosd(α) - offset for i in eachindex(x)]
    return s, n
end
channel_coords(m::MeshGeom; kw...) = channel_coords(m.lon, m.lat; kw...)

"""
    channel_masks(s, n; cell, L, W) -> (zone, chan, flood)

Cell masks: `zone` the measurement reach −0.3 L ≤ s ≤ 0.05 L; `chan` cells
there that overlap the channel (|n| ≤ W/2 + cell/2); `flood` cells there
wholly outside it (|n| > W/2 + cell). `cell`: the nominal cell size (m).
"""
function channel_masks(s::AbstractVector, n::AbstractVector; cell::Real,
                       L = CHANNEL.L, W = CHANNEL.W, kw...)
    zone = [-0.3L <= si <= 0.05L for si in s]
    chan = [zone[i] && abs(n[i]) <= W / 2 + 0.5cell for i in eachindex(s)]
    flood = [zone[i] && abs(n[i]) > W / 2 + cell for i in eachindex(s)]
    return zone, chan, flood
end

"""
    channel_metrics(mesh, h0, h1; angle, offset, cell, z = mesh.elevation,
                    area = mesh.area, L, W, S, Z0, Q, NMAN, ...) -> (deta, spill, unsteady)

Channel metrics from the last two frames' depths `h0`, `h1` (mesh cell
order): `deta` the mean of h + z − (Z0 − S s + h_n) over the channel cells
with h1 > 1 cm (m), `spill` the area of floodplain cells with h1 > 5 cm per
metre of the 0.35 L reach (m; needs areas), `unsteady` the largest |h1 − h0|
over those channel cells (m). NaN for `deta`/`unsteady` if no channel cell is
wet.
"""
function channel_metrics(mesh::MeshGeom, h0::AbstractVector, h1::AbstractVector;
                         angle::Real, offset::Real, cell::Real,
                         z::AbstractVector = mesh.elevation, area::AbstractVector = mesh.area,
                         L = CHANNEL.L, W = CHANNEL.W, S = CHANNEL.S, Z0 = CHANNEL.Z0,
                         Q = CHANNEL.Q, NMAN = CHANNEL.NMAN,
                         LON0 = CHANNEL.LON0, LAT0 = CHANNEL.LAT0, kw...)
    nc = length(mesh)
    length(h0) == length(h1) == length(z) == nc ||
        error("channel_metrics: depths and elevations must cover the mesh's $nc cells")
    length(area) == nc || error("channel_metrics: need cell areas (format 1: pass the mesh's parquet)")
    HN = channel_normal_depth(; Q, NMAN, W, S)
    s, n = channel_coords(mesh; angle, offset, LON0, LAT0)
    _, chan, flood = channel_masks(s, n; cell, L, W)
    wetc = [i for i in findall(chan) if h1[i] > 0.01]
    deta = isempty(wetc) ? NaN : mean(h1[i] + z[i] - (Z0 - S * s[i] + HN) for i in wetc)
    spill = sum(area[i] for i in findall(flood) if h1[i] > 0.05; init = 0.0) / (0.35L)
    unst = isempty(wetc) ? NaN : maximum(abs(h1[i] - h0[i]) for i in wetc)
    return (deta = deta, spill = spill, unsteady = unst)
end

"""
    parse_log_volumes(logpath) -> Vector{(t, vol_removed)}

The (time s, cumulative volume removed at the boundaries m³) of every step
line of a FloodModel run log (`t=…s … vol_removed=…m³`).
"""
function parse_log_volumes(logpath::AbstractString)
    out = NamedTuple{(:t, :vol_removed),Tuple{Float64,Float64}}[]
    for l in eachline(logpath)
        m = match(r"t=([0-9.]+)s.*vol_removed=([-0-9.e]+)m³", l)
        m === nothing && continue
        push!(out, (t = parse(Float64, m[1]), vol_removed = parse(Float64, m[2])))
    end
    return out
end

"""
    log_volume_at(entries, t) -> (t_line, vol_removed)

The `parse_log_volumes` entry nearest time `t` (the first if tied); (NaN, NaN)
if there are none.
"""
function log_volume_at(entries::AbstractVector, t::Real)
    best = (Inf, NaN, NaN)
    for e in entries
        abs(e.t - t) < best[1] && (best = (abs(e.t - t), e.t, e.vol_removed))
    end
    return best[2], best[3]
end
log_volume_at(logpath::AbstractString, t::Real) = log_volume_at(parse_log_volumes(logpath), t)

"""
    log_outflow(logpath, t0, t1; Q = CHANNEL.Q) -> fraction

Outflow between the log's step lines nearest `t0` and `t1`, as a fraction of
the inflow `Q`: Δvol_removed / (Q Δt). NaN if the log has no such interval.
"""
function log_outflow(logpath::AbstractString, t0::Real, t1::Real; Q = CHANNEL.Q, kw...)
    e = parse_log_volumes(logpath)
    (ta, va), (tb, vb) = log_volume_at(e, t0), log_volume_at(e, t1)
    return tb > ta ? (vb - va) / (Q * (tb - ta)) : NaN
end

"""
    read_channel_manifest(path) -> Vector{NamedTuple}

The rows of a make_channel_meshes.jl `manifest.csv`: `case`, `mesh`, `res`,
`cell` (m), `angle` (°), `offset` (m), `mesh_file`.
"""
function read_channel_manifest(path::AbstractString)
    map(readlines(path)[2:end]) do l
        r = split(l, ",")
        (case = String(r[1]), mesh = String(r[2]), res = parse(Int, r[3]),
         cell = parse(Float64, r[4]), angle = parse(Float64, r[5]),
         offset = parse(Float64, r[6]), mesh_file = String(r[10]))
    end
end

"""
    channel_run_metrics(r::Run; angle, offset, cell, log = <r.path with .log>, kw...)
        -> (deta, spill, unsteady, outflow, t0, t1)

`channel_metrics` from the run's last two frames (times `t0`, `t1`), and the
outflow fraction over [t0, t1] from its log (NaN if there is no log). `kw`:
the `CHANNEL` constants.
"""
function channel_run_metrics(r::Run; angle::Real, offset::Real, cell::Real,
                             log::AbstractString = replace(r.path, ".h5" => ".log"), kw...)
    nf = nframes(r)
    nf >= 2 || error("$(r.label): channel metrics need ≥ 2 output frames")
    t0, t1 = r.t[nf-1], r.t[nf]
    m = channel_metrics(r.mesh, depth(r, nf - 1), depth(r, nf); angle, offset, cell, kw...)
    of = isfile(log) ? log_outflow(log, t0, t1; kw...) : NaN
    return (deta = m.deta, spill = m.spill, unsteady = m.unsteady, outflow = of, t0 = t0, t1 = t1)
end

"""
    channel_series_metrics(meshdir, rundir; tag = "chan", solvers = ("swe", "li", "limd"),
                           io = stdout, csv = joinpath(rundir, "metrics.csv"), match = :id, kw...)
        -> Vector{NamedTuple}

Metrics of every run `rundir/<tag>_<case>_<solver>.h5` of the cases in
`meshdir/manifest.csv` (a missing `limd` run is skipped silently, others with
a warning; runs with < 2 frames with a warning). Writes the rows to `csv`
(`nothing`: don't) as test/channel_bench/analyse_channel.jl always has, and
prints one line per run to `io` (`nothing`: don't). Format-1 runs take their
geometry from the manifest's `mesh_file`, cells matched by `match` (see
`mesh_geometry`). Returns the rows.
"""
function channel_series_metrics(meshdir::AbstractString, rundir::AbstractString;
                                tag::AbstractString = "chan",
                                solvers = ("swe", "li", "limd"),
                                io = stdout, csv = joinpath(rundir, "metrics.csv"),
                                match::Symbol = :id, kw...)
    rows = NamedTuple[]
    out = csv === nothing ? nothing : open(csv, "w")
    try
        out === nothing ||
            println(out, "case,mesh,res,cell_m,angle,offset_m,solver,deta,spill,unsteady,outflow")
        for c in read_channel_manifest(joinpath(meshdir, "manifest.csv"))
            for solver in solvers     # limd: li + --momentum-damping 0.05, if run
                p = joinpath(rundir, "$(tag)_$(c.case)_$(solver).h5")
                isfile(p) || (solver == "limd" || @warn "missing $p"; continue)
                nfr = try FO.n_frames(p) catch; 0 end
                nfr >= 2 || (@warn "no output frames in $p (failed run?)"; continue)
                run = load_run(p; parquet = FO.format_version(p) >= 2 ? nothing : c.mesh_file, match)
                m = channel_run_metrics(run; angle = c.angle, offset = c.offset, cell = c.cell, kw...)
                row = (case = c.case, mesh = c.mesh, res = c.res, cell_m = c.cell, angle = c.angle,
                       offset_m = c.offset, solver = String(solver), deta = m.deta,
                       spill = m.spill, unsteady = m.unsteady, outflow = m.outflow)
                push!(rows, row)
                out === nothing ||
                    @printf(out, "%s,%s,%d,%.3f,%.1f,%.3f,%s,%.5f,%.3f,%.5f,%.4f\n", row.case,
                            row.mesh, row.res, row.cell_m, row.angle, row.offset_m, row.solver,
                            row.deta, row.spill, row.unsteady, row.outflow)
                io === nothing ||
                    @printf(io, "%-18s %-3s Δη %+7.3f m  spill %7.1f m  unsteady %.4f  outflow %.3f\n",
                            row.case, row.solver, row.deta, row.spill, row.unsteady, row.outflow)
            end
        end
    finally
        out === nothing || close(out)
    end
    return rows
end
