# analysis/parquet_mesh.jl — format-1 fallback: mesh geometry from GeoParquet.
#
# Pure Julia, no Python: the file is read with GDAL's (Arrow) Parquet driver
# through ArchGDAL, which FloodA5 already depends on (root and ci projects).
# Parquet2.jl cannot be used: it does not support nested columns (the
# `neighbours` list<string> column, the SGS list<double> curves) and its
# schema reader loses the columns after the first nested one (`geometry`).
#
# What it reproduces of the model's loader (A5Grid._load_from_geoparquet,
# FloodModel.initialise_flow_model) — and so of a format-2 file's /mesh:
#   ids        the `cell_id` column as stored (the output's `cell_ids`)
#   polygons   the WKB exterior ring as stored, i.e. shapely's `exterior`
#              coordinates: closed, first vertex repeated (= `cell.boundary`
#              of a loaded mesh, = format-2 `FO.cell_polygons`)
#   neighbours the `neighbours` column (ids, normalised to 16-digit hex as
#              `_norm_cell_id`), slot order as stored; ids not in the mesh
#              leave their slot 0. Without the column: shared-vertex
#              adjacency (two shared vertices at 1e-7°), sorted by row.
#   area       the model's `state.cell_area` (format 2 `cell_area`):
#              `sgs_cell_area` if the mesh has it (written by SGS builds and
#              by every --meshout since 2026-10), else `polygon_area_m2` of
#              the polygon. The older bridge's `cell_area` column is ignored,
#              as by the model (it may not be the per-cell area).
#   static     scalar numeric columns other than cell_id, resolution,
#              center_lon/lat, elevation, cell_area (null → NaN); list
#              columns (SGS curves, edge sills) are skipped.

using ArchGDAL
using JSON3

export load_mesh_parquet, polygon_area_m2

"""
    polygon_area_m2(poly::AbstractMatrix) -> Float64

Plan area (m²) of a 2 × nv (lon, lat) polygon: shoelace on an
equirectangular projection about the mean vertex latitude. Bit-identical to
`A5Grid._polygon_area_m2(boundary)` for `boundary = [poly[:, k] for k in 1:nv]`
(the closing vertex, if present, counts in the mean as it does there).
"""
function polygon_area_m2(poly::AbstractMatrix)
    n = size(poly, 2)
    n < 3 && return 0.0
    lat0 = mean(poly[2, k] for k in 1:n)
    cos_lat = cosd(lat0)
    area = 0.0
    j = n
    for i in 1:n
        xi = deg2rad(poly[1, i]) * R_EARTH * cos_lat
        yi = deg2rad(poly[2, i]) * R_EARTH
        xj = deg2rad(poly[1, j]) * R_EARTH * cos_lat
        yj = deg2rad(poly[2, j]) * R_EARTH
        area += (xj + xi) * (yj - yi)
        j = i
    end
    return abs(area) / 2.0
end

"`_norm_cell_id` (FloodModel.jl): 16-digit zero-padded hex if the id parses as hex."
function _norm_id(id::AbstractString)
    u = tryparse(UInt64, id; base = 16)
    return u === nothing ? String(id) : string(u; base = 16, pad = 16)
end

const _PQ_BASE = ("cell_id", "resolution", "center_lon", "center_lat")
const _PQ_NOT_STATIC = ("elevation", "cell_area", "cell_ids", "center_lons", "center_lats",
                        "elevations", "polygon_offsets", "polygon_lon", "polygon_lat",
                        "neighbours", "edge_cells", "source")

"""
    load_mesh_parquet(path) -> MeshGeom

Mesh geometry from a FloodA5 GeoParquet mesh, in the file's row order: ids,
centres, elevation (0 if the mesh has none, as the model runs it), cell areas, polygons (closed
rings, as format-2 output), neighbours (max_nb × n, 1-based, 0 = none),
scalar static columns, and the id scheme. Areas are what the model uses:
`sgs_cell_area` if present, else `polygon_area_m2`. Reads with GDAL
(ArchGDAL), no Python; list columns (SGS curves) are not read.
"""
function load_mesh_parquet(path::AbstractString)
    isfile(path) || error("no such file: $path")
    # GDAL's Arrow reader deadlocks inside a Julia process (Windows, GDAL 3.12)
    # unless it decodes on the calling thread: GDAL_NUM_THREADS = 1 while reading.
    old = ArchGDAL.getconfigoption("GDAL_NUM_THREADS", "")
    ArchGDAL.setconfigoption("GDAL_NUM_THREADS", "1")
    try
        return _load_mesh_parquet(path)
    finally
        isempty(old) ? ArchGDAL.clearconfigoption("GDAL_NUM_THREADS") :
                       ArchGDAL.setconfigoption("GDAL_NUM_THREADS", old)
    end
end

function _load_mesh_parquet(path::AbstractString)
    ArchGDAL.read(path) do ds
        lyr = ArchGDAL.getlayer(ds, 0)
        fd = ArchGDAL.layerdefn(lyr)
        nf = ArchGDAL.nfield(fd)
        names = [ArchGDAL.getname(ArchGDAL.getfielddefn(fd, i)) for i in 0:nf-1]
        types = [ArchGDAL.gettype(ArchGDAL.getfielddefn(fd, i)) for i in 0:nf-1]
        col(name) = (i = findfirst(==(name), names); i === nothing ? -1 : i - 1)
        ic, ilon, ilat, ires = col("cell_id"), col("center_lon"), col("center_lat"), col("resolution")
        (ic < 0 || ilon < 0 || ilat < 0) &&
            error("$path: not a FloodA5 mesh (needs cell_id, center_lon, center_lat)")
        inb = col("neighbours")
        inb >= 0 && types[inb + 1] != ArchGDAL.OFTStringList && (inb = -1)
        numeric = (ArchGDAL.OFTReal, ArchGDAL.OFTInteger, ArchGDAL.OFTInteger64)
        scal = [i - 1 for i in 1:nf if types[i] in numeric && !(names[i] in _PQ_BASE)]
        # Do not decode the (large) list columns we skip.
        skip = [names[i] for i in 1:nf if types[i] in (ArchGDAL.OFTRealList,
                ArchGDAL.OFTIntegerList, ArchGDAL.OFTInteger64List) ||
                (types[i] == ArchGDAL.OFTStringList && i - 1 != inb)]
        _ignore_fields!(lyr, skip)

        n = ArchGDAL.nfeature(lyr)
        ids = Vector{String}(undef, n); lon = Vector{Float64}(undef, n); lat = similar(lon)
        res = zeros(Int, n)
        vals = Dict(i => fill(NaN, n) for i in scal)
        polys = Vector{Matrix{Float64}}(undef, n)
        nbids = Vector{Vector{String}}(undef, n)
        getd(f, i) = ArchGDAL.isfieldsetandnotnull(f, i) ? ArchGDAL.asdouble(f, i) : NaN
        r = 0
        ArchGDAL.resetreading!(lyr)
        for f in lyr
            r += 1
            ids[r] = ArchGDAL.asstring(f, ic)
            lon[r] = getd(f, ilon); lat[r] = getd(f, ilat)
            ires >= 0 && ArchGDAL.isfieldsetandnotnull(f, ires) && (res[r] = ArchGDAL.asint64(f, ires))
            for i in scal
                vals[i][r] = getd(f, i)
            end
            nbids[r] = inb >= 0 && ArchGDAL.isfieldsetandnotnull(f, inb) ?
                       copy(ArchGDAL.asstringlist(f, inb)) : String[]
            g = ArchGDAL.getgeom(f)
            ring = ArchGDAL.getgeom(g, 0)
            nv = ArchGDAL.ngeom(ring)
            P = Matrix{Float64}(undef, 2, nv)
            for k in 1:nv
                P[1, k] = ArchGDAL.getx(ring, k - 1); P[2, k] = ArchGDAL.gety(ring, k - 1)
            end
            polys[r] = P
        end
        r == n || error("$path: read $r of $n rows")

        static = Dict{String,Vector{Float64}}(names[i+1] => vals[i] for i in scal)
        elevation = haskey(static, "elevation") ? static["elevation"] : zeros(n)   # no DEM: the model runs at z = 0
        area = haskey(static, "sgs_cell_area") ? copy(static["sgs_cell_area"]) :
               [polygon_area_m2(P) for P in polys]
        for k in _PQ_NOT_STATIC
            delete!(static, k)
        end
        nb = inb >= 0 ? _neighbours_from_ids(ids, nbids) : _neighbours_shared_vertices(polys)
        scheme = _parquet_id_scheme(lyr, res)
        return MeshGeom(ids, lon, lat, elevation, area, polys, nb, static, scheme)
    end
end

"Ask GDAL not to decode the fields `names` (OGR_L_SetIgnoredFields; a speed-up only)."
function _ignore_fields!(lyr, names::Vector{String})
    isempty(names) && return
    GC.@preserve names begin
        ptrs = [Base.unsafe_convert(Cstring, Base.cconvert(Cstring, s)) for s in names]
        push!(ptrs, Cstring(C_NULL))
        ArchGDAL.GDAL.ogr_l_setignoredfields(lyr, ptrs)
    end
    return
end

"Neighbour matrix from per-row neighbour id lists (as `_build_adjacency_matrix!`)."
function _neighbours_from_ids(ids::Vector{String}, nbids::Vector{Vector{String}})
    pos = Dict(_norm_id(id) => i for (i, id) in enumerate(ids))
    max_nb = maximum(length, nbids; init = 0)
    nb = zeros(Int32, max_nb, length(ids))
    for i in eachindex(nbids), (s, id) in enumerate(nbids[i])
        nb[s, i] = get(pos, _norm_id(id), 0)
    end
    return nb
end

"Shared-vertex adjacency (as `_build_adjacency_shared_vertices`), neighbours sorted by row."
function _neighbours_shared_vertices(polys::Vector{Matrix{Float64}})
    key(P, k) = (round(Int64, P[1, k] * 1e7), round(Int64, P[2, k] * 1e7))
    verts(P) = (nv = size(P, 2); nv > 1 && P[:, 1] ≈ P[:, nv] ? (1:nv-1) : (1:nv))
    v2c = Dict{Tuple{Int64,Int64},Vector{Int}}()
    for (i, P) in enumerate(polys), k in verts(P)
        push!(get!(v2c, key(P, k), Int[]), i)
    end
    lists = Vector{Vector{Int}}(undef, length(polys))
    for (i, P) in enumerate(polys)
        shared = Dict{Int,Int}()
        for k in verts(P), j in get(v2c, key(P, k), Int[])
            j == i || (shared[j] = get(shared, j, 0) + 1)
        end
        lists[i] = sort([j for (j, c) in shared if c == 2])
    end
    nb = zeros(Int32, maximum(length, lists; init = 0), length(polys))
    for i in eachindex(lists), (s, j) in enumerate(lists[i])
        nb[s, i] = j
    end
    return nb
end

"The mesh's id scheme: the `flooda5` file metadata, else as the model infers it."
function _parquet_id_scheme(lyr, res::Vector{Int})
    for m in ArchGDAL.metadata(lyr; domain = "_PARQUET_METADATA_")
        startswith(m, "flooda5=") || continue
        s = get(JSON3.read(m[9:end]), :a5_id_scheme, nothing)
        s === nothing || return String(s)
    end
    return !isempty(res) && all(==(0), res) ? "cartesian" : "a5-v0"
end
