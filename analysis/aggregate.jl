# analysis/aggregate.jl — a fine run on a coarse mesh's cells: like-for-like
# comparison across resolutions, and of SGS with cell-mean models.
#
# Why. A coarse model and a fine reference do not hold their water on the same
# cells, and an SGS model does not mean the same by "depth" as a cell-mean
# model: SGS depth is measured from the cell's lowest ground, and an SGS cell
# with a little water in a sub-grid channel has a depth, where a cell-mean
# model holding the same volume spreads it as a film over the whole cell.
# So the comparison is made per coarse cell, through quantities that mean the
# same in every model:
#   volume          stored water (m³)
#   wetted fraction the share of the cell's plan area where the local water
#                   depth exceeds a threshold — fine grid: the area of its wet
#                   cells; cell-mean model: 0 or 1; SGS: from its storage table
#   water level     fine grid: the area-weighted level of its wet cells; any
#                   model: its own level (cell-mean: bed + depth, SGS: η);
#                   compared only where both are wet
# Fine cells are assigned to the coarse polygon containing their centre.
# (Plan of 2026-11 §9.3a; first done in a scratch script, fine_compare.jl.)
#
# SGS storage tables are not in the output file (its /mesh carries the scalar
# columns only); `sgs_tables` reads them from the run's mesh GeoParquet.

export assign_to_cells, aggregate_fine, aggregate_onto, SGSAreaTable, sgs_wet_area,
       load_sgs_tables, load_mesh_columns, sgs_tables, wet_fraction_cells, wet_fraction_sgs, wetted_fraction,
       compare_aggregated, compare_to_fine

"""
    assign_to_cells(lon, lat, polygons; bucket = 50.0, lon0 = nothing, lat0 = nothing)
        -> Vector{Int}
    assign_to_cells(fine::MeshGeom, coarse::MeshGeom; bucket = 50.0)

For each point (lon, lat) (degrees), the index of the polygon that contains
it, 0 if none (even–odd rule; a point on a shared side goes to the first
polygon found). `polygons`: per cell a 2 × nv (lon, lat) matrix, closed or
not. Polygons and points are projected equirectangularly about (lon0, lat0)
(default: the mean polygon vertex) and the polygons sorted into square
buckets of `bucket` m, so each point is tested against a few polygons only.
"""
function assign_to_cells(lon::AbstractVector, lat::AbstractVector, polygons::AbstractVector;
                         bucket::Real = 50.0, lon0 = nothing, lat0 = nothing)
    length(lon) == length(lat) || error("assign_to_cells: lon and lat lengths differ")
    if lon0 === nothing || lat0 === nothing
        nv = sum(size(P, 2) for P in polygons; init = 0)
        nv > 0 || return zeros(Int, length(lon))
        lon0 = something(lon0, sum(sum(view(P, 1, :)) for P in polygons) / nv)
        lat0 = something(lat0, sum(sum(view(P, 2, :)) for P in polygons) / nv)
    end
    mx = R_EARTH * π / 180 * cosd(lat0); my = R_EARTH * π / 180
    B = Float64(bucket)
    xs = [(P[1, :] .- lon0) .* mx for P in polygons]
    ys = [(P[2, :] .- lat0) .* my for P in polygons]
    buckets = Dict{Tuple{Int,Int},Vector{Int}}()
    for k in eachindex(polygons)
        isempty(xs[k]) && continue
        for bx in floor(Int, minimum(xs[k]) / B):floor(Int, maximum(xs[k]) / B),
            by in floor(Int, minimum(ys[k]) / B):floor(Int, maximum(ys[k]) / B)
            push!(get!(buckets, (bx, by), Int[]), k)
        end
    end
    out = zeros(Int, length(lon))
    none = Int[]
    Threads.@threads for i in eachindex(lon)
        px = (lon[i] - lon0) * mx; py = (lat[i] - lat0) * my
        for k in get(buckets, (floor(Int, px / B), floor(Int, py / B)), none)
            if _in_polygon(px, py, xs[k], ys[k])
                out[i] = k
                break
            end
        end
    end
    return out
end
assign_to_cells(fine::MeshGeom, coarse::MeshGeom; bucket::Real = 50.0) =
    assign_to_cells(fine.lon, fine.lat, coarse.polygons; bucket)

"Even–odd point-in-polygon test (vertices x, y; closed or not)."
function _in_polygon(px::Real, py::Real, x::AbstractVector, y::AbstractVector)
    c = false; n = length(x)
    for i in 1:n
        j = i == 1 ? n : i - 1
        if ((y[i] > py) != (y[j] > py)) && (px < (x[j] - x[i]) * (py - y[i]) / (y[j] - y[i]) + x[i])
            c = !c
        end
    end
    return c
end

"""
    aggregate_fine(owner, coarse_area, V, h, η, A; wet = 0.1, cover = 0.9) -> NamedTuple

A fine field on coarse cells. `owner[i]`: the coarse cell of fine cell i (0:
none; `assign_to_cells`); `coarse_area`: coarse cell areas (m²); `V`, `h`,
`η`, `A`: the fine cells' volume (m³), depth and water level (m) and area
(m²). Fine cells with a non-finite volume take no part. Per coarse cell:

- `volume`: Σ V of its fine cells (m³);
- `fine_area`: Σ A of its fine cells (m²);
- `wet_area`: Σ A of its fine cells deeper than `wet` (m²);
- `wet_frac`: wet_area / fine_area (NaN where fine_area = 0);
- `level`: Σ A·η / Σ A over its wet fine cells (NaN if none);
- `covered`: the fine cells cover at least `cover` of the coarse cell's area
  (cells the fine grid does not cover are left out of comparisons).
"""
function aggregate_fine(owner::AbstractVector{<:Integer}, coarse_area::AbstractVector,
                        V::AbstractVector, h::AbstractVector, η::AbstractVector, A::AbstractVector;
                        wet::Real = 0.1, cover::Real = 0.9)
    length(owner) == length(V) == length(h) == length(η) == length(A) ||
        error("aggregate_fine: the fine vectors differ in length")
    n = length(coarse_area)
    volume = zeros(n); fine_area = zeros(n); wet_area = zeros(n); wl = zeros(n)
    for i in eachindex(owner)
        k = owner[i]
        (k == 0 || !isfinite(V[i])) && continue
        fine_area[k] += A[i]
        volume[k] += V[i]
        if isfinite(h[i]) && h[i] > wet && isfinite(η[i])
            wet_area[k] += A[i]; wl[k] += A[i] * η[i]
        end
    end
    wet_frac = [fine_area[k] > 0 ? wet_area[k] / fine_area[k] : NaN for k in 1:n]
    level = [wet_area[k] > 0 ? wl[k] / wet_area[k] : NaN for k in 1:n]
    covered = [coarse_area[k] > 0 && fine_area[k] >= cover * coarse_area[k] for k in 1:n]
    return (; volume, fine_area, wet_area, wet_frac, level, covered)
end

"""
    aggregate_onto(fine::Run, coarse::MeshGeom, k; wet = 0.1, cover = 0.9, owner = nothing)
        -> NamedTuple

`aggregate_fine` for frame `k` of `fine` on the cells of `coarse`, with the
fine run's `volume`, depth and `wse`. `owner` (from `assign_to_cells(fine.mesh,
coarse)`) can be passed to reuse the assignment across frames; the result
carries it (`owner`) and the fine frame's time (`t`).
"""
function aggregate_onto(fine::Run, coarse::MeshGeom, k::Integer; wet::Real = 0.1,
                        cover::Real = 0.9, owner = nothing)
    isempty(coarse.polygons) && error("aggregate_onto: the coarse mesh has no polygons")
    isempty(fine.mesh.area) && error("aggregate_onto: $(fine.label) has no cell areas")
    own = owner === nothing ? assign_to_cells(fine.mesh, coarse) : owner
    a = aggregate_fine(own, coarse.area, volume(fine, k), depth(fine, k), wse(fine, k),
                       fine.mesh.area; wet, cover)
    return merge(a, (owner = own, t = fine.t[k]))
end

# ─────────────────────────── SGS wetted area ───────────────────────────────

"""
    SGSAreaTable

The part of a cell's SGS storage table that gives its wetted area: the level
knots `elev_bins` and the cumulative wetted area `area_curve` at them (m²),
the cell's area and its lowest and highest ground (`A5Grid.SGSTable`).
"""
struct SGSAreaTable
    elev_bins::Vector{Float64}
    area_curve::Vector{Float64}
    cell_area::Float64
    z_min::Float64
    z_max::Float64
end

"""
    sgs_wet_area(t::SGSAreaTable, η) -> m²

Plan area of the cell with ground below level η, as the model computes it
(`A5Grid.wetted_area_from_wse`, the same arithmetic): 0 at or below the
lowest ground, the whole cell at or above the highest, linear between knots.
The area where the local depth exceeds d is `sgs_wet_area(t, η - d)`.
"""
function sgs_wet_area(t::SGSAreaTable, η::Real)
    η = Float64(η)
    η <= t.z_min && return 0.0
    η >= t.z_max && return t.cell_area
    eb = t.elev_bins; lo, hi = 1, length(eb)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        eb[mid] <= η ? lo = mid : hi = mid
    end
    frac = (η - eb[lo]) / (eb[hi] - eb[lo] + 1e-15)
    return t.area_curve[lo] + frac * (t.area_curve[hi] - t.area_curve[lo])
end

"""
    load_sgs_tables(parquet; ids = nothing) -> Vector{SGSAreaTable}

The SGS wetted-area tables of a mesh GeoParquet (columns `sgs_elev_bins`,
`sgs_area_curve`, `sgs_cell_area`, `sgs_z_min`, `sgs_z_max`), read with GDAL
(no Python), in the file's row order, or in the order of `ids` (hex cell
ids, e.g. a run's `r.mesh.ids`; an error if one is missing). Cells without
tables (no DEM) get NaN extremes.
"""
function load_sgs_tables(path::AbstractString; ids = nothing)
    cols = ("sgs_elev_bins", "sgs_area_curve", "sgs_cell_area", "sgs_z_min", "sgs_z_max")
    c = load_mesh_columns(path, cols; ids)
    missing_cols = [n for n in cols if !haskey(c, n)]
    isempty(missing_cols) ||
        error("$path: no SGS tables (missing $(join(missing_cols, ", "))); build the mesh with --flow-model sgs")
    return [SGSAreaTable(c["sgs_elev_bins"][i], c["sgs_area_curve"][i], c["sgs_cell_area"][i],
                         c["sgs_z_min"][i], c["sgs_z_max"][i]) for i in eachindex(c["sgs_z_min"])]
end

"""
    load_mesh_columns(parquet, names; ids = nothing) -> Dict(name => Vector)

Columns `names` of a mesh GeoParquet, read with GDAL (no Python): numeric
columns as `Vector{Float64}` (null → NaN), list columns (SGS curves and
per-face tables, e.g. `sgs_face_raise`, one entry per neighbour slot) as
`Vector{Vector{Float64}}` (null → empty). Rows in file order, or in the order
of `ids` (hex cell ids; an error if one is missing). Columns the file does
not have are left out of the result.
"""
function load_mesh_columns(path::AbstractString, names; ids = nothing)
    isfile(path) || error("no such file: $path")
    old = ArchGDAL.getconfigoption("GDAL_NUM_THREADS", "")
    ArchGDAL.setconfigoption("GDAL_NUM_THREADS", "1")      # see load_mesh_parquet
    fid, cols = try
        _read_mesh_columns(path, String.(collect(names)))
    finally
        isempty(old) ? ArchGDAL.clearconfigoption("GDAL_NUM_THREADS") :
                       ArchGDAL.setconfigoption("GDAL_NUM_THREADS", old)
    end
    ids === nothing && return cols
    pos = Dict(_norm_id(id) => i for (i, id) in enumerate(fid))
    perm = [get(pos, _norm_id(id), 0) for id in ids]
    j = findfirst(==(0), perm)
    j === nothing || error("$path: no cell $(ids[j])")
    return Dict{String,Any}(k => v[perm] for (k, v) in cols)
end

function _read_mesh_columns(path::AbstractString, want::Vector{String})
    ArchGDAL.read(path) do ds
        lyr = ArchGDAL.getlayer(ds, 0)
        fd = ArchGDAL.layerdefn(lyr)
        nf = ArchGDAL.nfield(fd)
        names = [ArchGDAL.getname(ArchGDAL.getfielddefn(fd, i)) for i in 0:nf-1]
        types = [ArchGDAL.gettype(ArchGDAL.getfielddefn(fd, i)) for i in 0:nf-1]
        ic = findfirst(==("cell_id"), names)
        ic === nothing && error("$path: not a FloodA5 mesh (no cell_id)")
        sel = [i for i in eachindex(names) if names[i] in want]
        _ignore_fields!(lyr, [names[i] for i in eachindex(names) if i != ic && !(i in sel)])
        islist(i) = types[i] in (ArchGDAL.OFTRealList, ArchGDAL.OFTIntegerList, ArchGDAL.OFTInteger64List)
        n = ArchGDAL.nfeature(lyr)
        fid = Vector{String}(undef, n)
        cols = Dict{String,Any}(names[i] => (islist(i) ? Vector{Vector{Float64}}(undef, n) : fill(NaN, n))
                                for i in sel)
        r = 0
        ArchGDAL.resetreading!(lyr)
        for f in lyr
            r += 1
            fid[r] = ArchGDAL.asstring(f, ic - 1)
            for i in sel
                set = ArchGDAL.isfieldsetandnotnull(f, i - 1)
                if islist(i)
                    cols[names[i]][r] = !set ? Float64[] :
                        types[i] == ArchGDAL.OFTRealList ? Float64.(ArchGDAL.asdoublelist(f, i - 1)) :
                        types[i] == ArchGDAL.OFTIntegerList ? Float64.(ArchGDAL.asintlist(f, i - 1)) :
                        Float64.(ArchGDAL.asint64list(f, i - 1))
                elseif set
                    cols[names[i]][r] = ArchGDAL.asdouble(f, i - 1)
                end
            end
        end
        r == n || error("$path: read $r of $n rows")
        return fid, cols
    end
end

"""
    sgs_tables(r::Run; parquet = nothing) -> Vector{SGSAreaTable}

The SGS wetted-area tables of an SGS run's cells, in the output's cell order,
from `parquet` or else the mesh the run loaded (`--meshload` in its recorded
command). An error if the file's lowest ground (`sgs_z_min`) differs from the
one the output recorded: the mesh has changed since the run.
"""
function sgs_tables(r::Run; parquet = nothing)
    path = parquet
    if path === nothing
        argv = get(r.info, "argv", String[])
        i = findfirst(==("--meshload"), argv)
        (i === nothing || i == length(argv)) &&
            error("$(r.label): no --meshload in the run's command; pass parquet = MESH")
        path = String(argv[i+1])
        isfile(path) || error("$(r.label): its mesh $path is not here; pass parquet = MESH")
    end
    t = load_sgs_tables(String(path); ids = r.mesh.ids)
    zr = get(r.mesh.static, "sgs_z_min", nothing)
    if zr !== nothing
        for i in eachindex(t)
            a, b = t[i].z_min, zr[i]
            (isnan(a) && isnan(b)) || abs(a - b) <= 1e-9 * max(1.0, abs(b)) ||
                error("$(r.label): $path is not the mesh of this run (cell $i: z_min $a, run $b)")
        end
    end
    return t
end

"""
    wet_fraction_cells(h; wet = 0.1) -> Vector

Wetted fraction of cell-mean cells: 1 where the depth exceeds `wet`, else 0.
"""
wet_fraction_cells(h::AbstractVector; wet::Real = 0.1) =
    [isfinite(x) && x > wet ? 1.0 : 0.0 for x in h]

"""
    wet_fraction_sgs(tables, η; wet = 0.1) -> Vector

Wetted fraction of SGS cells at levels η: the share of each cell's area where
η − ground > `wet`, from its table (0 for cells without a table or level).
"""
wet_fraction_sgs(tables::AbstractVector{SGSAreaTable}, η::AbstractVector; wet::Real = 0.1) =
    [(isfinite(η[i]) && isfinite(tables[i].z_min) && tables[i].cell_area > 0) ?
     sgs_wet_area(tables[i], η[i] - wet) / tables[i].cell_area : 0.0 for i in eachindex(η)]

"""
    wetted_fraction(r::Run, k; wet = 0.1, tables = nothing) -> Vector

The share of each cell's area wetter than `wet` (local depth) at frame `k`:
SGS runs from their tables (`tables`, default `sgs_tables(r)`), other runs
0 or 1 by the cell's depth.
"""
function wetted_fraction(r::Run, k::Integer; wet::Real = 0.1, tables = nothing)
    r.flow_model == "sgs" || return wet_fraction_cells(depth(r, k); wet)
    return wet_fraction_sgs(tables === nothing ? sgs_tables(r) : tables, wse(r, k); wet)
end

"""
    compare_aggregated(agg, V, η, fm, A; both_wet = 0.5, mean_depth = 0.1) -> NamedTuple

A coarse model against a fine run aggregated onto its cells (`agg`, from
`aggregate_fine`). `V`, `η`, `fm`, `A`: the model's volume, water level,
wetted fraction (`wetted_fraction`) and cell area. Over the cells the fine
grid covers where the model's volume is finite:

- `storage`, `storage_fine` (m³) and `storage_ratio` = storage / storage_fine;
- `wet_area`, `wet_area_fine`: Σ A·fraction (m²), and the fuzzy CSI
  Σ A·min(f, f_fine) / Σ A·max(f, f_fine);
- `level_median`, `level_rms`, `level_mean` (m): model level − fine level over
  the `n_level` cells where both wetted fractions are ≥ `both_wet`;
- `depth_rms` (m): rms of the mean-depth (V/A) difference over the
  `n_depth` cells where either holds ≥ `mean_depth` of mean depth.
"""
function compare_aggregated(agg, V::AbstractVector, η::AbstractVector, fm::AbstractVector,
                            A::AbstractVector; both_wet::Real = 0.5, mean_depth::Real = 0.1)
    n = length(A)
    length(V) == length(η) == length(fm) == length(agg.volume) == n ||
        error("compare_aggregated: lengths differ")
    ok = [agg.covered[i] && isfinite(V[i]) for i in 1:n]
    st = sum(V[ok]; init = 0.0); stf = sum(agg.volume[ok]; init = 0.0)
    f5 = [isfinite(x) ? x : 0.0 for x in agg.wet_frac]
    wa = sum(A[i] * fm[i] for i in 1:n if ok[i]; init = 0.0)
    waf = sum(A[i] * f5[i] for i in 1:n if ok[i]; init = 0.0)
    mx = sum(A[i] * max(fm[i], f5[i]) for i in 1:n if ok[i]; init = 0.0)
    mn = sum(A[i] * min(fm[i], f5[i]) for i in 1:n if ok[i]; init = 0.0)
    d = [η[i] - agg.level[i] for i in 1:n
         if ok[i] && f5[i] >= both_wet && fm[i] >= both_wet && isfinite(η[i]) && isfinite(agg.level[i])]
    dv = [(V[i] - agg.volume[i]) / A[i] for i in 1:n
          if ok[i] && A[i] > 0 && (V[i] / A[i] >= mean_depth || agg.volume[i] / A[i] >= mean_depth)]
    nan_if_empty(f, x) = isempty(x) ? NaN : f(x)
    return (storage = st, storage_fine = stf, storage_ratio = st / stf,
            wet_area = wa, wet_area_fine = waf, csi = mx > 0 ? mn / mx : NaN,
            level_median = nan_if_empty(median, d), level_mean = nan_if_empty(mean, d),
            level_rms = nan_if_empty(x -> sqrt(mean(abs2, x)), d), n_level = length(d),
            depth_rms = nan_if_empty(x -> sqrt(mean(abs2, x)), dv), n_depth = length(dv))
end

"""
    compare_to_fine(fine::Run, model::Run, t; tol = 300.0, wet = 0.1, cover = 0.9,
                    both_wet = 0.5, owner = nothing, tables = nothing) -> NamedTuple

`compare_aggregated` of `model` against `fine` at time `t` (frames within
`tol` s of it), the fine run aggregated onto the model's cells. Wetted
fractions at local depth `wet` (SGS: `tables`, default `sgs_tables(model)`).
Pass `owner` (`assign_to_cells(fine.mesh, model.mesh)`) and `tables` to reuse
them across calls. The result adds `t_fine`, `t_model` and `owner`.
"""
function compare_to_fine(fine::Run, model::Run, t::Real; tol::Real = 300.0, wet::Real = 0.1,
                         cover::Real = 0.9, both_wet::Real = 0.5, owner = nothing, tables = nothing)
    kf = frame_at(fine, t; tol); km = frame_at(model, t; tol)
    agg = aggregate_onto(fine, model.mesh, kf; wet, cover, owner)
    fm = wetted_fraction(model, km; wet, tables)
    c = compare_aggregated(agg, volume(model, km), wse(model, km), fm, model.mesh.area; both_wet)
    return merge(c, (t_fine = fine.t[kf], t_model = model.t[km], owner = agg.owner))
end
