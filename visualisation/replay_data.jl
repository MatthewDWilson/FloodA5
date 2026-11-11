# visualisation/replay_data.jl — FloodA5Replay: data for replaying an output file
#
# The Makie-free half of the viewer's replay mode (`MakieVisualiser.replay`,
# `visualisation/replay.jl`): everything that turns an output file into what
# the live window is fed, so CI can test it without GLMakie
# (test/test_replay.jl).
#
#   include("visualisation/replay_data.jl"); const RP = FloodA5Replay
#   rd = RP.open_replay("run.h5")                    # format 2
#   rd = RP.open_replay("old.h5"; parquet = "m.parquet")   # format 1
#   fr = rd.frame(k)        # depth, saturation, volume, velocity, vel_u, vel_v
#   rd.series               # per-frame volume budget, wet cells (all frames)
#   close(rd)
#
# What is read and what is computed:
#   - map fields: read from the file, frame by frame (one frame cached);
#     a variable the file lacks is shown as zeros and listed in `rd.notes`.
#   - domain volume Σ `volume` and wet cells (depth > 1e-4 m, as the live
#     viewer): computed from the file for every frame.
#   - input volume: reconstructed from the run record (format 2): uniform
#     rainfall × t × Σ area (cells ≥ 1 m²), point sources, and the embedded
#     hydrograph / BCI / BDY files, integrated as the model's mass-balance
#     logger (`total_cumulative_input`). It is checked against the recorded
#     `/summary` total at the last frame and dropped (NaN, with a note) if it
#     does not match or cannot be built — never guessed.
#   - removed volume: the file stores it only at the end of the run. Runs
#     with `--closed-boundaries` remove nothing (0 at every frame); for open
#     boundaries the viewer shows added − domain (= removed + mass-balance
#     error) and the recorded end-of-run values as markers.

module FloodA5Replay

using HDF5, Printf

if isdefined(Main, :FloodA5Analysis)
    const FA = Main.FloodA5Analysis
else
    include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
    const FA = FloodA5Analysis
end
const FO = FA.FO

# The model's own time-series parsers (CSV, LISFLOOD-FP .bdy/.bci). They need
# `BCIEntry` from boundaryinputs/boundary_conditions.jl, which pulls in the
# mesh code; this is a field-for-field copy (test/test_replay.jl checks it).
struct BCIEntry
    boundary_type :: Char
    x1            :: Float64
    y1            :: Float64
    x2            :: Float64
    y2            :: Float64
    bc_code       :: String
    bc_value      :: String
end
include(joinpath(@__DIR__, "..", "boundaryinputs", "timeseries_io.jl"))

export ReplayData, open_replay, mesh_adaptor, read_frame, input_model, input_volume

# ─────────────────────────────── mesh adaptor ──────────────────────────────

"One cell as `MakieVisualiser.start` reads it (`id`, `center_lon/lat`, `boundary`)."
struct ReplayCell
    id         :: String
    center_lon :: Float64
    center_lat :: Float64
    boundary   :: Vector{Vector{Float64}}     # closed ring of [lon, lat]
end

"The mesh object `MakieVisualiser.start` takes: `cells` and `resolution`."
struct ReplayMesh
    cells      :: Vector{ReplayCell}
    resolution :: Int
end
Base.length(m::ReplayMesh) = length(m.cells)

"Mean Earth surface area (m²) on the sphere of `FA.R_EARTH`."
const _EARTH_AREA = 4π * FA.R_EARTH^2

"""
    a5_resolution(g::FA.MeshGeom) -> Int

A5 resolution level of the mesh from its median cell area (an A5 cell at
level r ≥ 1 is 1/(60·4^(r-1)) of the sphere), or 0 for Cartesian meshes and
when the areas match no level within 25% — as `mesh.resolution` of a live
run (0 for Cartesian). Only shown in the diagnostics panel.
"""
function a5_resolution(g::FA.MeshGeom)
    g.id_scheme == "cartesian" && return 0
    a = filter(x -> isfinite(x) && x > 0, g.area)
    isempty(a) && return 0
    med = sort(a)[cld(length(a), 2)]
    r = 1 + log(4, _EARTH_AREA / (60 * med))
    rr = round(Int, r)
    return (rr >= 1 && abs(4.0^(r - rr) - 1) < 0.25) ? rr : 0
end

"""
    mesh_adaptor(g::FA.MeshGeom; resolution = a5_resolution(g)) -> ReplayMesh

Cells in the output file's order with their polygons (closed rings). Needs
the polygons: format-2 files carry them; format-1 files need `parquet`.
"""
function mesh_adaptor(g::FA.MeshGeom; resolution::Integer = a5_resolution(g))
    isempty(g.polygons) &&
        error("no cell polygons: a format-1 output file needs its mesh (parquet = …, --mesh)")
    cells = Vector{ReplayCell}(undef, length(g))
    for i in 1:length(g)
        P = g.polygons[i]
        ring = [[P[1, j], P[2, j]] for j in axes(P, 2)]
        (isempty(ring) || ring[1] == ring[end]) || push!(ring, ring[1])
        cells[i] = ReplayCell(g.ids[i], g.lon[i], g.lat[i], ring)
    end
    return ReplayMesh(cells, Int(resolution))
end

"""
    adjacency(g::FA.MeshGeom) -> Dict{String,Vector{String}}

Cell id → neighbour ids (the live viewer's `flow_state.adjacency`, used for
its ring-index mode); empty if the file has no neighbour table.
"""
function adjacency(g::FA.MeshGeom)
    adj = Dict{String,Vector{String}}()
    size(g.neighbours, 1) == 0 && return adj
    for i in 1:length(g)
        adj[g.ids[i]] = [g.ids[j] for j in g.neighbours[:, i] if j > 0]
    end
    return adj
end

"The model's `_find_nearest_cell`: nearest centre in plain degrees."
nearest_cell(g::FA.MeshGeom, lon::Real, lat::Real) =
    argmin((g.lon .- lon) .^ 2 .+ (g.lat .- lat) .^ 2)

_f(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)

# ────────────────────────────── recorded inputs ────────────────────────────

"""
    InputModel

The run's water inputs as the model's mass-balance logger integrates them
(`total_cumulative_input`, boundaryinputs/sources.jl): a constant rate
(injection points, rain points, BCI QFIX), hydrographs (inflow points, BCI
QVAR; trapezoidal, flat beyond the knots) and uniform rainfall over the
cells of area ≥ 1 m².
"""
struct InputModel
    const_rate  :: Float64
    hydrographs :: Vector{Tuple{Vector{Float64},Vector{Float64}}}
    rain_rate   :: Float64            # m/s
    rain_area   :: Float64            # m²
end

"`cumulative_volume(::InflowPoint, t)` (boundaryinputs/sources.jl)."
function _hydrograph_volume(t_s::Vector{Float64}, Q::Vector{Float64}, t::Float64)
    isempty(t_s) && return 0.0
    t <= 0.0 && return 0.0
    t <= t_s[1] && return max(0.0, Q[1]) * t
    vol = 0.0
    t_s[1] > 0.0 && (vol += max(0.0, Q[1]) * t_s[1])
    for k in 1:length(t_s)-1
        t1, t2 = t_s[k], t_s[k+1]
        Q1, Q2 = max(0.0, Q[k]), max(0.0, Q[k+1])
        if t2 <= t
            vol += 0.5 * (Q1 + Q2) * (t2 - t1)
        else
            Qt = Q1 + (t - t1) / (t2 - t1) * (Q2 - Q1)
            return vol + 0.5 * (Q1 + Qt) * (t - t1)
        end
    end
    return vol + max(0.0, Q[end]) * (t - t_s[end])
end

"Cumulative input volume (m³) from t = 0 to `t`."
function input_volume(m::InputModel, t::Real)
    t = Float64(t)
    vol = m.const_rate * t
    for (ts, q) in m.hydrographs
        vol += _hydrograph_volume(ts, q, t)
    end
    m.rain_rate > 0.0 && (vol += m.rain_rate * t * m.rain_area)
    return vol
end

# Path of the run's input file recorded under `flag` (and `argv_value`, if
# given): the embedded copy written into `dir`, or the original file if it
# was not embedded and is unchanged. `nothing` if there is none.
function _input_file(f::HDF5.File, info, flag, dir; argv_value = nothing)
    gi = f["run/inputs"]
    for k in sort(collect(keys(gi)))
        a = FO.attrs_dict(gi[k])
        a["flag"] == flag || continue
        argv_value === nothing || strip(a["argv_value"]) == strip(argv_value) || continue
        if a["embedded"] && haskey(gi[k], "bytes")
            d = joinpath(dir, k); mkpath(d)
            p = joinpath(d, basename(replace(a["original_path"], '\\' => '/')))
            write(p, read(gi[k]["bytes"]))
            return p
        end
        p = a["original_path"]
        (isfile(p) && FO.sha256_file(p) == a["sha256"]) && return p
        return nothing
    end
    return nothing
end

"""
    input_model(r::FA.Run) -> (InputModel or nothing, note::String)

Rebuild the run's inputs from its run record (format 2). `nothing` with the
reason if that is not possible (format 1, an input file neither embedded
nor still on disk unchanged, an unreadable series).
"""
function input_model(r::FA.Run)
    r.format >= 2 || return nothing, "format-1 file: no run record, inputs unknown"
    o = get(r.info, "options", Dict{String,Any}())
    isempty(o) && return nothing, "no options in the run record"
    g = r.mesh
    rate = 0.0
    hyd = Tuple{Vector{Float64},Vector{Float64}}[]
    for s in something(get(o, "injection_points", nothing), [])   # (lon, lat, m³/s)
        rate += _f(s[3])
    end
    for s in something(get(o, "rainpoints", nothing), [])         # (lon, lat, mm/h)
        isempty(g.area) && return nothing, "rain points need cell areas"
        rate += _f(s[3]) / 3_600_000.0 * g.area[nearest_cell(g, _f(s[1]), _f(s[2]))]
    end
    rain = _f(something(get(o, "rainfall_rate_m_per_s", 0.0), 0.0))
    area = 0.0
    if rain > 0.0
        isempty(g.area) && return nothing, "rainfall needs cell areas"
        area = sum(a for a in g.area if a >= 1.0; init = 0.0)
    end
    dir = mktempdir()
    try
        h5open(r.path, "r") do f
            for s in something(get(o, "inflow_points", nothing), [])   # (lat, lon, file[, label])
                path = String(s[3])
                p = _input_file(f, r.info, "--inflow-point", dir; argv_value = path)
                p === nothing && error("inflow hydrograph $path is neither embedded nor on disk unchanged")
                label = length(s) >= 4 ? String(s[4]) : splitext(basename(path))[1]
                series = read_timeseries(auto_timeseries_reader(p))
                key = haskey(series, label) ? label : first(keys(series))
                push!(hyd, series[key])
            end
            bci = get(o, "inflow_bci", nothing)
            if bci !== nothing
                pb = _input_file(f, r.info, "--inflow-bci", dir)
                pb === nothing && error("BCI file $bci is neither embedded nor on disk unchanged")
                # the .bdy (given or found next to the .bci) goes next to the .bci copy
                pd = _input_file(f, r.info, "--inflow-bdy", dirname(pb))
                entries, series = parse_bci_file(pb; bdy_path = pd)
                for e in entries
                    e.boundary_type == 'P' || continue
                    if e.bc_code == "QVAR" && haskey(series, e.bc_value)
                        push!(hyd, series[e.bc_value])
                    elseif e.bc_code == "QFIX"
                        q = tryparse(Float64, e.bc_value)
                        q === nothing || (rate += q)
                    end
                end
            end
        end
    catch err
        return nothing, "inputs not rebuilt: " * sprint(showerror, err)
    finally
        rm(dir; recursive = true, force = true)
    end
    return InputModel(rate, hyd, rain, area), "inputs rebuilt from the run record"
end

"""
    source_cells(r::FA.Run) -> (indices, note)

Cells of the run's point sources (injection points, rain points, inflow
points, BCI point entries with QVAR/QFIX), found as the model finds them.
Points given in projected coordinates (`--bc-epsg`) are skipped (converting
them needs PROJ); the note says so.
"""
function source_cells(r::FA.Run)
    r.format >= 2 || return Int[], ""
    o = get(r.info, "options", Dict{String,Any}())
    g = r.mesh
    idx = Int[]
    for s in something(get(o, "injection_points", nothing), [])
        push!(idx, nearest_cell(g, _f(s[1]), _f(s[2])))
    end
    for s in something(get(o, "rainpoints", nothing), [])
        push!(idx, nearest_cell(g, _f(s[1]), _f(s[2])))
    end
    projected = get(o, "bc_epsg", nothing) !== nothing
    note = ""
    inflow = something(get(o, "inflow_points", nothing), [])
    if projected && (!isempty(inflow) || get(o, "inflow_bci", nothing) !== nothing)
        note = "inflow points in projected coordinates (--bc-epsg): not highlighted"
    elseif !projected
        for s in inflow
            push!(idx, nearest_cell(g, _f(s[2]), _f(s[1])))
        end
        if get(o, "inflow_bci", nothing) !== nothing
            dir = mktempdir()
            try
                h5open(r.path, "r") do f
                    pb = _input_file(f, r.info, "--inflow-bci", dir)
                    pb === nothing && return
                    entries, _ = parse_bci_file(pb; bdy_path = _input_file(f, r.info, "--inflow-bdy", dirname(pb)))
                    for e in entries
                        e.boundary_type == 'P' && e.bc_code in ("QVAR", "QFIX") &&
                            push!(idx, nearest_cell(g, e.x1, e.y1))
                    end
                end
            catch
            finally
                rm(dir; recursive = true, force = true)
            end
        end
    end
    return idx, note
end

# ───────────────────────────────── series ──────────────────────────────────

const WET_DEPTH = 1e-4      # m, the live viewer's wet-cell threshold

"""
    replay_series(r::FA.Run, f::HDF5.File) -> (series::NamedTuple, notes)

The time-series strip for every frame: `t`, `vol_added`, `vol_domain`,
`vol_removed`, `mb_err`, `n_wet` (NaN where not known), `removed_label`,
`mb_title`, `end_mark` (`nothing`, or the recorded end-of-run
`(t, vol_removed, mb_err)` when the per-frame values are not known),
`inputs_known`, `removed` (`:closed`, `:end_only` or `:unknown`).
"""
function replay_series(r::FA.Run, f::HDF5.File)
    nf = length(r.t)
    notes = String[]
    vol_domain = zeros(nf); n_wet = zeros(nf)
    hasdepth = FO.has_var(f, "water_depth")
    hasvol = FO.has_var(f, "volume")
    for k in 1:nf
        if hasvol
            v = FO.read_frame(f, "volume", k)
            vol_domain[k] = sum(x for x in v if isfinite(x); init = 0.0)
        else
            vol_domain[k] = NaN
        end
        hasdepth && (n_wet[k] = count(>(WET_DEPTH), FO.read_frame(f, "water_depth", k)))
    end
    hasvol || push!(notes, "no `volume` in the file: domain volume unknown")

    summ = get(r.info, "summary", Dict{String,Any}())
    model, note = input_model(r)
    vol_added = fill(NaN, nf)
    if model !== nothing
        vol_added .= input_volume.(Ref(model), r.t)
        if haskey(summ, "vol_added_m3") && get(summ, "completed", false) && nf > 0
            rec = Float64(summ["vol_added_m3"])
            if !isapprox(vol_added[end], rec; rtol = 1e-6, atol = 1e-6)
                vol_added .= NaN
                note = @sprintf("rebuilt inputs (%.6g m³ at the last frame) do not match the recorded %.6g m³: not shown",
                                input_volume(model, r.t[end]), rec)
            else
                note *= " (match the recorded total)"
            end
        end
    end
    push!(notes, note)

    opts = get(r.info, "options", Dict{String,Any}())
    closed = get(opts, "closed_boundaries", false) === true
    removed = closed ? :closed : (haskey(summ, "vol_removed_m3") ? :end_only : :unknown)
    vol_removed = closed ? zeros(nf) : fill(NaN, nf)
    end_mark = nothing
    removed_label = "Removed"
    mb_title = "Mass balance error"
    if closed
        mb_err = vol_added .- vol_domain
    else
        # removed volume is not stored per frame: show added − domain
        # (= removed + mass-balance error) in its place
        vol_removed = vol_added .- vol_domain
        removed_label = "Removed + MB err"
        mb_err = fill(NaN, nf)
        mb_title = "Mass balance error (end of run)"
        if removed === :end_only && nf > 0
            end_mark = (Float64(r.t[end]), Float64(summ["vol_removed_m3"]),
                        Float64(get(summ, "mass_balance_error_m3", NaN)))
        end
        push!(notes, "removed volume is stored only at the end of the run" *
                     (removed === :end_only ? " (marker)" : " (not recorded)") *
                     "; the strip shows added − domain")
    end
    series = (t = copy(r.t), vol_added = vol_added, vol_domain = vol_domain,
              vol_removed = vol_removed, mb_err = mb_err, n_wet = n_wet,
              removed_label = removed_label, mb_title = mb_title, end_mark = end_mark,
              inputs_known = model !== nothing && all(isfinite, vol_added),
              removed = removed)
    return series, notes
end

# ───────────────────────────────── frames ──────────────────────────────────

const FRAME_VARS = ("water_depth", "saturation", "volume", "velocity", "vel_u", "vel_v")

"""
    ReplayData

An output file opened for replay: `run` (`FA.Run`), `mesh` (`ReplayMesh`),
`t`, `series` (see `replay_series`), `notes`, `source_indices`, `adjacency`,
`label`, `missing_vars`, and `frame(k)` → the fields of frame k. Holds the
HDF5 file open: `close(rd)` when done.
"""
mutable struct ReplayData
    run            :: FA.Run
    mesh           :: ReplayMesh
    t              :: Vector{Float64}
    series         :: NamedTuple
    notes          :: Vector{String}
    source_indices :: Vector{Int}
    adjacency      :: Dict{String,Vector{String}}
    label          :: String
    missing_vars   :: Vector{String}
    h5             :: HDF5.File
    cache_k        :: Int
    cache          :: Any
    frame          :: Function
end

Base.close(rd::ReplayData) = (isopen(rd.h5) && close(rd.h5); nothing)
Base.show(io::IO, rd::ReplayData) = print(io, "ReplayData(\"", rd.label, "\", ",
    length(rd.t), " frames, ", length(rd.mesh), " cells)")

"""
    read_frame(rd, k) -> NamedTuple

Frame `k`'s `depth`, `saturation`, `volume`, `velocity`, `vel_u`, `vel_v`
(Float64, file cell order). Variables the file lacks are zeros (`velocity`
falls back to |(vel_u, vel_v)|). The last frame read is cached.
"""
function read_frame(rd::ReplayData, k::Integer)
    1 <= k <= length(rd.t) || error("frame $k out of range 1:$(length(rd.t))")
    rd.cache_k == k && return rd.cache
    n = length(rd.mesh)
    get1(v) = v in rd.missing_vars ? zeros(n) : Float64.(FO.read_frame(rd.h5, v, k))
    u = get1("vel_u"); v = get1("vel_v")
    vel = "velocity" in rd.missing_vars ? hypot.(u, v) : get1("velocity")
    fr = (depth = get1("water_depth"), saturation = get1("saturation"),
          volume = get1("volume"), velocity = vel, vel_u = u, vel_v = v)
    rd.cache_k = k; rd.cache = fr
    return fr
end

"""
    open_replay(path; parquet = nothing, match = :id, label = basename(path)) -> ReplayData

Open an output file for `MakieVisualiser.replay`. Format-2 files carry the
mesh; a format-1 file needs `parquet` (its GeoParquet mesh; `match` as in
`FA.mesh_geometry`). Reads the per-frame series of the whole run (one pass
over `volume` and `water_depth`).
"""
function open_replay(path::AbstractString; parquet = nothing, match::Symbol = :id,
                     label::AbstractString = basename(path))
    r = FA.load_run(path; parquet, match, label)
    length(r.t) > 0 || error("$path has no frames")
    mesh = mesh_adaptor(r.mesh)
    h5 = h5open(path, "r")
    try
        series, notes = replay_series(r, h5)
        missing_vars = [v for v in FRAME_VARS if !FO.has_var(h5, v)]
        isempty(missing_vars) ||
            push!(notes, "not in the file (shown as 0): " * join(missing_vars, ", "))
        src, snote = source_cells(r)
        isempty(snote) || push!(notes, snote)
        rd = ReplayData(r, mesh, copy(r.t), series, filter(!isempty, notes), src,
                        adjacency(r.mesh), String(label), missing_vars, h5, 0, nothing,
                        identity)
        rd.frame = k -> read_frame(rd, k)
        return rd
    catch
        close(h5)
        rethrow()
    end
end

end # module
