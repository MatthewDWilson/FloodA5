# analysis/runs.jl — Run and MeshGeom: loading, frames, fields, coordinates.

export MeshGeom, Run, load_run, mesh_geometry, nframes, times, frame_at,
       matched_frames, field, has_field, series, depth, volume, speed, wse,
       storage, stored_volume, local_xy, source_cell, check_frame, print_frames,
       R_EARTH

"Mean Earth radius (m), as `A5Grid._EARTH_R`."
const R_EARTH = 6_371_000.0

"""
    MeshGeom

Mesh geometry in the output file's cell order.

- `ids`, `lon`, `lat`, `elevation`: always present.
- `area` (m²): empty if unknown (format-1 file without a GeoParquet).
- `polygons`: per cell a 2 × nv matrix of (lon, lat) vertices; empty if unknown.
- `neighbours`: max_nb × n, 1-based, 0 = none; 0 × n if unknown.
- `static`: other per-cell scalar variables (e.g. `sgs_z_min`, `friction`).
- `id_scheme`: "a5-v0", "a5-v1", "cartesian", … or "".
"""
struct MeshGeom
    ids::Vector{String}            # hex cell ids, as in the output file
    lon::Vector{Float64}
    lat::Vector{Float64}
    elevation::Vector{Float64}
    area::Vector{Float64}
    polygons::Vector{Matrix{Float64}}
    neighbours::Matrix{Int32}
    static::Dict{String,Vector{Float64}}
    id_scheme::String
end

Base.length(m::MeshGeom) = length(m.ids)
Base.show(io::IO, m::MeshGeom) = print(io, "MeshGeom(", length(m), " cells",
    isempty(m.area) ? "" : ", areas", isempty(m.polygons) ? "" : ", polygons",
    size(m.neighbours, 1) == 0 ? "" : ", neighbours", ")")

"""
    Run

A loaded output file: frame times, mesh geometry and the run record. Frame
data are read on demand (`field`, `depth`, …).

- `format`: output format (1 = FloodA5 ≤ 0.3.x, 2 = 0.4.0+).
- `info`: `FO.run_info` for format 2 (argv, options, summary, …); for format 1
  only what the file has.
- `flow_model`: "sgs", "standard", "swe" or "" (unknown).
"""
struct Run
    path::String
    label::String
    format::Int
    t::Vector{Float64}
    mesh::MeshGeom
    info::Dict{String,Any}
    flow_model::String
end

Base.show(io::IO, r::Run) = print(io, "Run(\"", r.label, "\", format ", r.format, ", ",
    length(r.t), " frames, t = ", isempty(r.t) ? "–" : "$(r.t[1])…$(r.t[end]) s", ", ",
    length(r.mesh), " cells", isempty(r.flow_model) ? "" : ", " * r.flow_model, ")")

"""
    mesh_geometry(h5path; parquet = nothing, match = :id) -> MeshGeom

Format 2: everything from the file's `/mesh`. Format 1: the cell table from
the file; area, polygons, neighbours and static variables from `parquet`
(a GeoParquet path, or a `MeshGeom` already loaded with `load_mesh_parquet`
to reuse it across runs) if given. Cells are matched to the file by id
(`match = :id`, an error if the sets differ) or by row (`match = :row`,
for runs on a re-labelled copy of the mesh, e.g. pya5 0.x vs 1.x ids).
"""
function mesh_geometry(h5path::AbstractString; parquet = nothing, match::Symbol = :id)
    h5open(h5path, "r") do f
        ct = FO.cell_table(f)
        ids = String.(ct.ids)
        scheme = haskey(attributes(f["mesh"]), "id_scheme") ?
                 String(read(attributes(f["mesh"])["id_scheme"])) : ""
        g = f["mesh"]
        # per-cell scalar static variables (both formats store them in /mesh)
        reserved = ("cell_ids", "center_lons", "center_lats", "elevations", "cell_area",
                    "polygon_offsets", "polygon_lon", "polygon_lat", "neighbours",
                    "edge_cells", "source")
        static = Dict{String,Vector{Float64}}()
        for k in keys(g)
            k in reserved && continue
            g[k] isa HDF5.Dataset || continue
            v = read(g[k])
            v isa AbstractVector{<:Real} && length(v) == length(ids) &&
                (static[k] = Float64.(v))
        end
        if FO.format_version(f) >= 2
            area = haskey(g, "cell_area") ? Float64.(FO.cell_area(f)) : Float64[]
            polys = haskey(g, "polygon_offsets") ? FO.cell_polygons(f) : Matrix{Float64}[]
            nb = haskey(g, "neighbours") ? Int32.(FO.cell_neighbours(f)) :
                 zeros(Int32, 0, length(ids))
            return MeshGeom(ids, Float64.(ct.lon), Float64.(ct.lat), Float64.(ct.elevation),
                            area, polys, nb, static, scheme)
        end
        # format 1: the model's cell area is sgs_cell_area when the mesh has it
        # (as load_mesh_parquet); polygons and neighbours need the GeoParquet
        area = haskey(static, "sgs_cell_area") ? copy(static["sgs_cell_area"]) : Float64[]
        base = MeshGeom(ids, Float64.(ct.lon), Float64.(ct.lat), Float64.(ct.elevation),
                        area, Matrix{Float64}[], zeros(Int32, 0, length(ids)), static, scheme)
        parquet === nothing && return base
        pm = parquet isa MeshGeom ? parquet : load_mesh_parquet(String(parquet))
        return _merge_parquet(base, pm, match)
    end
end

"Attach a GeoParquet mesh's geometry to a format-1 cell table (see `mesh_geometry`)."
function _merge_parquet(base::MeshGeom, pm::MeshGeom, match::Symbol)
    n = length(base)
    length(pm) == n || error("mesh has $(length(pm)) cells, the output file $n")
    if match === :id
        pos = Dict(id => i for (i, id) in enumerate(pm.ids))
        perm = [get(pos, id, 0) for id in base.ids]       # file cell i = parquet row perm[i]
        any(==(0), perm) && error("the output's cell ids are not those of the mesh " *
                                  "(a re-labelled copy? use match = :row)")
    elseif match === :row
        perm = collect(1:n)
    else
        error("match must be :id or :row")
    end
    inv = invperm(perm)
    nb = similar(pm.neighbours)
    for i in 1:n, s in axes(nb, 1)
        j = pm.neighbours[s, perm[i]]
        nb[s, i] = j == 0 ? Int32(0) : Int32(inv[j])
    end
    return MeshGeom(base.ids, base.lon, base.lat, base.elevation,
                    isempty(pm.area) ? Float64[] : pm.area[perm],
                    isempty(pm.polygons) ? Matrix{Float64}[] : pm.polygons[perm], nb,
                    Dict(k => v[perm] for (k, v) in pm.static), base.id_scheme)
end

"""
    load_run(path; label = basename(path), parquet = nothing, match = :id,
             flow_model = nothing) -> Run

Load an output file's frame times, mesh geometry and run record. `parquet`
and `match`: see `mesh_geometry` (format 1 only). `flow_model` overrides the
recorded one (format-1 files do not record it).
"""
function load_run(path::AbstractString; label::AbstractString = basename(path),
                  parquet = nothing, match::Symbol = :id, flow_model = nothing)
    isfile(path) || error("no such file: $path")
    fmt = FO.format_version(path)
    t = Float64.(FO.frame_times(path))
    mesh = mesh_geometry(path; parquet, match)
    info = fmt >= 2 ? Dict{String,Any}(FO.run_info(path)) : Dict{String,Any}()
    fm = flow_model !== nothing ? String(flow_model) :
         String(get(info, "flow_model", get(get(info, "options", Dict()), "flow_model", "")))
    return Run(String(path), String(label), fmt, t, mesh, info, fm)
end

nframes(r::Run) = length(r.t)
times(r::Run) = r.t

"""
    frame_at(r, t; tol = Inf) -> k

Index of the frame nearest time `t`; an error if it is more than `tol` s away.
"""
function frame_at(r::Run, t::Real; tol::Real = Inf)
    nframes(r) > 0 || error("$(r.label): no frames")
    k = argmin(abs.(r.t .- t))
    abs(r.t[k] - t) <= tol ||
        error("$(r.label): nearest frame to t = $t s is at $(r.t[k]) s (tol $tol s)")
    return k
end

"""
    matched_frames(a, b; tol = 60.0) -> Vector{Tuple{Int,Int}}

Frame pairs `(ka, kb)` at matched simulation times: for each frame of `a`,
the nearest frame of `b` if it is within `tol` s.
"""
function matched_frames(a::Run, b::Run; tol::Real = 60.0)
    out = Tuple{Int,Int}[]
    nframes(b) == 0 && return out
    for (ka, t) in enumerate(a.t)
        kb = argmin(abs.(b.t .- t))
        abs(b.t[kb] - t) <= tol && push!(out, (ka, kb))
    end
    return out
end

"Per-cell variable `var` at frame `k` (default: the last), as Float64, file cell order."
field(r::Run, var::AbstractString, k::Integer = nframes(r)) =
    Float64.(FO.read_frame(r.path, var, k))

has_field(r::Run, var::AbstractString) = FO.has_var(r.path, var)

"Time series of `var` at `cells` (n_frames × length(cells))."
series(r::Run, var::AbstractString, cells) = FO.read_series(r.path, var, cells)

depth(r::Run, k::Integer = nframes(r))  = field(r, "water_depth", k)
volume(r::Run, k::Integer = nframes(r)) = field(r, "volume", k)

"Cell speed (m/s): |(vel_u, vel_v)| if written, else the `velocity` variable."
function speed(r::Run, k::Integer = nframes(r))
    has_field(r, "vel_u") || return field(r, "velocity", k)
    return hypot.(field(r, "vel_u", k), field(r, "vel_v", k))
end

"""
    wse(r, k; sgs = (r.flow_model == "sgs")) -> Vector

Water-surface elevation. SGS output measures depth from the cell's lowest DEM
point, so WSE = `sgs_z_min` + depth there; otherwise mean elevation + depth.
"""
function wse(r::Run, k::Integer = nframes(r); sgs::Bool = r.flow_model == "sgs")
    h = depth(r, k)
    if sgs
        haskey(r.mesh.static, "sgs_z_min") ||
            error("$(r.label): SGS WSE needs sgs_z_min (format 1: pass the mesh's parquet)")
        return r.mesh.static["sgs_z_min"] .+ h
    end
    return r.mesh.elevation .+ h
end

"""
    storage(r, k) -> m³

Water stored in the domain. Standard and SWE runs: Σ depth·area over cells
with a finite depth (needs cell areas). SGS runs: Σ `volume` — SGS depth is
measured from the cell's lowest DEM point, so depth·area overstates the
store (planar embankment: 263 800 vs 212 990 m³).
"""
function storage(r::Run, k::Integer = nframes(r))
    r.flow_model == "sgs" && return stored_volume(r, k)
    isempty(r.mesh.area) && error("$(r.label): no cell areas (format 1: pass the mesh's parquet)")
    h = depth(r, k)
    return sum(h[i] * r.mesh.area[i] for i in eachindex(h) if isfinite(h[i]); init = 0.0)
end

"Σ of the `volume` variable over cells where it is finite (m³)."
stored_volume(r::Run, k::Integer = nframes(r)) = sum(filter(isfinite, volume(r, k)); init = 0.0)

"""
    local_xy(lon, lat, lon0, lat0; R = R_EARTH) -> (x, y)
    local_xy(mesh, lon0, lat0; R = R_EARTH)

Equirectangular coordinates (m) about (lon0, lat0): x east, y north.
"""
function local_xy(lon::AbstractVector, lat::AbstractVector, lon0::Real, lat0::Real; R::Real = R_EARTH)
    x = deg2rad.(lon .- lon0) .* R .* cosd(lat0)
    y = deg2rad.(lat .- lat0) .* R
    return x, y
end
local_xy(m::MeshGeom, lon0::Real, lat0::Real; R::Real = R_EARTH) = local_xy(m.lon, m.lat, lon0, lat0; R)

"""
    source_cell(mesh, lon, lat) -> i

The cell whose centre is nearest (lon, lat) (equirectangular distance) — the
cell where FloodA5 injects a point source at that location.
"""
source_cell(m::MeshGeom, lon::Real, lat::Real) =
    argmin((m.lon .- lon) .^ 2 .* cosd(lat)^2 .+ (m.lat .- lat) .^ 2)

# ─────────────────────────── frame helpers ─────────────────────────────────

"""
    check_frame(r, k; hint = "")

Error unless 1 ≤ k ≤ nframes(r); the message lists the available frames by
their old zero-padded /frames/NNNNNN labels, followed by `hint`.
"""
function check_frame(r::Run, k::Integer; hint::AbstractString = "")
    names = [lpad(string(j), 6, '0') for j in 1:nframes(r)]
    1 <= k <= length(names) ||
        error("Frame $(lpad(string(k), 6, '0')) not found in $(r.path). Available: " *
              "$(first(names, 3))…$(last(names, 3)) ($(length(names)) total)." * hint)
    return nothing
end

"""
    print_frames([io,] r)

Table of the frames: index (zero-padded), time, wet cells (saturation ≥ 1)
and maximum depth (NaN without `water_depth`) — for choosing a frame.
"""
function print_frames(io::IO, r::Run)
    println(io, "Frames in $(r.path):")
    @printf(io, "%8s  %10s  %8s  %10s\n", "frame#", "t (s)", "n_wet", "max_depth")
    hasdepth = has_field(r, "water_depth")
    for k in 1:nframes(r)
        n_wet = count(>=(1.0), field(r, "saturation", k))
        d = hasdepth ? field(r, "water_depth", k) : Float64[]
        md = isempty(d) ? NaN : maximum(d)
        @printf(io, "%8s  %10.1f  %8d  %10.4f\n", lpad(string(k), 6, '0'), r.t[k], n_wet, md)
    end
end
print_frames(r::Run) = print_frames(stdout, r)
