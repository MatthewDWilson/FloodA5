# mesh/CartesianGrid.jl
# ----------------------
# Pure-Julia Cartesian (4-connected, orthogonal quadrilateral) mesh
# generator — the "cartesian-formulation" branch's second mesh backend,
# alongside A5 (mesh/A5Grid.jl, via the pya5 Python bridge).
#
# See FloodA5_Generalisation_BenchmarkingPlan.md §5.2 and
# FloodA5_NSIDES_Generalisation_Handover.md §7/§8 for the design discussion
# this implements. Summary of the decisions made there:
#
#   - No Python bridge, no external meshing library (Meshes.jl considered
#     and deliberately not used — see the handover's ID-scheme discussion;
#     the same reasoning applies to geometry generation, which is a plain
#     double loop for a regular grid, not a problem a computational-
#     geometry library earns its keep on).
#   - Reuses A5Cell/A5Mesh directly (both already structurally generic —
#     nothing A5-specific in either struct). No new mesh-container type.
#   - Cell IDs are A5Grid._to_hex((UInt64(row)<<32)|UInt64(col)) — 16-char
#     hex, indistinguishable in SHAPE from an A5 ID, trivially invertible,
#     and requires zero changes to any ID-handling code already in the
#     pipeline (the _norm_cell_id/_to_hex_safe hardening from the N_SIDES
#     refactor means a non-hex ID would *also* work now, but there's no
#     reason to choose one — hex-shaped IDs were free and consistent).
#   - Adjacency is populated EXPLICITLY and EXACTLY from (row,col)
#     arithmetic, not left for the shared-vertex geometric fallback to
#     discover — a regular grid's connectivity is known in closed form, so
#     there's no reason to pay for (or risk) floating-point vertex matching.
#   - This generates a FULL RECTANGLE over the AOI's bounding box, not an
#     AOI-polygon-CLIPPED fill the way A5's mesh_for_aoi works. This is
#     deliberate: the primary purpose of this backend is exact structural
#     correspondence with a raster grid (e.g. matching a LISFLOOD-FP ASCII
#     grid cell-for-cell for the Phase 1 benchmark), which requires a
#     complete rectangle — a jagged polygon-clipped boundary would be the
#     wrong shape to compare against a raster. NaN-elevation cells (e.g.
#     outside an irregular real AOI, or outside DEM coverage) are already
#     handled correctly everywhere downstream (Bug 60's isfinite guards,
#     the same mechanism that already handles A5 boundary/NaN cells) — so
#     an irregular AOI can still be represented, just as a rectangle with
#     some cells left NaN, exactly like a real-world ASCII raster's NODATA
#     cells would be.
#   - Geometry is generated directly in decimal degrees (lon/lat), matching
#     every other part of this pipeline. Cell width/height in degrees is
#     derived from a target metre size via a single reference latitude
#     (domain mid-latitude by default) — this makes cells "square in
#     metres" AT that reference latitude only; a domain with significant
#     north-south extent will show the same lon/lat-vs-metres distortion
#     any lon/lat grid has, growing with distance from lat_ref. This is a
#     known, deliberate limitation for this branch — true projected-CRS
#     generation (for exact LISFLOOD-FP raster alignment against a raster
#     that is itself in a projected CRS) is explicitly deferred; see the
#     handover document.

using JSON3

"""
    generate_cartesian_mesh(aoi_path_or_str, dx_m; lat_ref=nothing,
                             output_path=nothing, format=:geoparquet) → A5Mesh

Generate a rectangular, 4-connected Cartesian grid mesh covering the
bounding box of `aoi_path_or_str` (a GeoJSON file path or GeoJSON string —
same accepted forms as `A5Grid.mesh_for_aoi`), with square cells
approximately `dx_m` metres on a side.

**This clips to the AOI's bounding box, not its polygon shape** — see the
module-level comment above for why. If you need the exact bounds of an
existing raster (e.g. to match a LISFLOOD-FP ASCII grid), pass a GeoJSON
whose bounding box equals that raster's extent, or use
[`generate_cartesian_mesh_bounds`](@ref) directly with the raster's own
`(xllcorner, yllcorner, ncols, nrows, cellsize)` header values instead of
going through an AOI file at all — that avoids a bounding-box round-trip
that could reintroduce a fraction-of-a-cell offset.

`lat_ref` (degrees) sets the reference latitude used to convert `dx_m`
into longitude/latitude cell extents. Defaults to the AOI bounding box's
own mid-latitude. See the module-level comment for why cells are only
exactly square in metres at this one latitude.

`output_path`/`format` behave identically to `mesh_for_aoi`: if given, the
generated mesh is also saved to that path (GeoParquet if the path ends in
`.parquet`, else GeoJSON) before being returned.

Returns an `A5Mesh` with `adjacency` already populated (4-connected grid
adjacency, computed directly from row/col arithmetic — not left for the
shared-vertex fallback). `mesh.resolution` is set to `0`
(the A5 resolution-level concept doesn't apply here — see `cell_size_m`
in `static_vars` for the actual cell size, and use `dx_m` at call time for
anything that needs it before the mesh exists).
"""
function generate_cartesian_mesh(aoi_path_or_str::String, dx_m::Float64;
                                  lat_ref::Union{Float64,Nothing}=nothing,
                                  output_path::Union{String,Nothing}=nothing,
                                  format::Symbol=:geoparquet)::A5Mesh
    dx_m > 0.0 || error("generate_cartesian_mesh: dx_m must be positive (got $dx_m)")

    geojson_str = A5Grid.load_aoi(aoi_path_or_str)
    geojson     = JSON3.read(geojson_str, Dict)
    geometry    = A5Grid._aoi_geometry(geojson)
    min_lon, min_lat, max_lon, max_lat = A5Grid._geometry_bbox(geometry)

    mesh = generate_cartesian_mesh_bounds(min_lon, min_lat, max_lon, max_lat, dx_m;
                                           lat_ref=lat_ref)
    mesh.aoi_geojson = geojson_str

    if output_path !== nothing
        A5Grid.save_mesh(mesh, output_path; format=format)   # writes .prj sidecar internally
    end
    return mesh
end

"""
    generate_cartesian_mesh_bounds(min_lon, min_lat, max_lon, max_lat, dx_m;
                                    lat_ref=nothing) → A5Mesh

Lower-level entry point behind [`generate_cartesian_mesh`](@ref): generate
the grid directly from an explicit bounding box, with no AOI-file parsing.
This is the function to call when the bounds must match an existing raster
exactly (e.g. a LISFLOOD-FP ASCII grid's own header) — pass that raster's
own `xllcorner`/`yllcorner`/`ncols`/`nrows`/`cellsize`-derived bounds
directly rather than round-tripping through a GeoJSON AOI file, to avoid
any risk of a fraction-of-a-cell bounding-box mismatch.

Cells are generated for every `row` from `0` to `nrows-1` and every `col`
from `0` to `ncols-1`, `row` counting north from `min_lat`, `col`
counting east from `min_lon` — a complete rectangle, not clipped to any
polygon. `nrows`/`ncols` are chosen so the rectangle covers
`[min_lon,max_lon] × [min_lat,max_lat]` (rounding up — the final row/column
may extend fractionally past `max_lon`/`max_lat` if the span isn't an exact
multiple of the cell size).
"""
function generate_cartesian_mesh_bounds(min_lon::Float64, min_lat::Float64,
                                         max_lon::Float64, max_lat::Float64,
                                         dx_m::Float64;
                                         lat_ref::Union{Float64,Nothing}=nothing)::A5Mesh
    dx_m > 0.0 || error("generate_cartesian_mesh_bounds: dx_m must be positive (got $dx_m)")
    max_lon > min_lon || error("generate_cartesian_mesh_bounds: max_lon must exceed min_lon")
    max_lat > min_lat || error("generate_cartesian_mesh_bounds: max_lat must exceed min_lat")

    lat0 = lat_ref !== nothing ? lat_ref : 0.5 * (min_lat + max_lat)
    R    = A5Grid._EARTH_R
    cos_lat0 = cosd(lat0)
    cos_lat0 > 1e-6 ||
        error("generate_cartesian_mesh_bounds: lat_ref=$lat0 is too close to a pole " *
              "for an equirectangular degree/metre conversion.")

    dlon_deg = rad2deg(dx_m / (R * cos_lat0))
    dlat_deg = rad2deg(dx_m / R)

    ncols = max(1, ceil(Int, (max_lon - min_lon) / dlon_deg))
    nrows = max(1, ceil(Int, (max_lat - min_lat) / dlat_deg))
    n     = nrows * ncols

    @info "Cartesian grid: $nrows rows × $ncols cols = $n cells " *
          "(dx≈$(round(dx_m,digits=1))m at lat_ref=$(round(lat0,digits=4))°, " *
          "dlon=$(round(dlon_deg,digits=6))°, dlat=$(round(dlat_deg,digits=6))°)"

    cells       = Vector{A5Cell}(undef, n)
    cell_size_m = fill(dx_m, n)
    adjacency   = Dict{String,Vector{String}}()
    sizehint!(adjacency, n)

    # id(r,c): row/col packed into a UInt64, then hex-encoded via A5Grid's
    # own helper — see the module-level comment for why this exact scheme.
    _cell_id(r::Int, c::Int) = A5Grid._to_hex((UInt64(r) << 32) | UInt64(c))

    @inbounds for r in 0:(nrows-1), c in 0:(ncols-1)
        idx = r * ncols + c + 1

        lon_w = min_lon + c * dlon_deg
        lon_e = lon_w + dlon_deg
        lat_s = min_lat + r * dlat_deg
        lat_n = lat_s + dlat_deg

        # Closed ring, CCW, first vertex repeated at the end — matches the
        # GeoJSON/A5 boundary convention every downstream function already
        # expects (_shared_edge, _polygon_area_m2, DiamondFlux.jl's
        # centroid computation, etc. all tolerate or expect this form).
        boundary = [
            [lon_w, lat_s],
            [lon_e, lat_s],
            [lon_e, lat_n],
            [lon_w, lat_n],
            [lon_w, lat_s],
        ]
        center_lon = 0.5 * (lon_w + lon_e)
        center_lat = 0.5 * (lat_s + lat_n)

        id = _cell_id(r, c)
        cells[idx] = A5Cell(id, 0, center_lon, center_lat, boundary)

        # Explicit, exact 4-connected adjacency — no shared-vertex geometric
        # matching needed or used (contrast with A5, where adjacency for a
        # loaded-without-precomputed-adjacency mesh falls back to matching
        # shared boundary vertices). A regular grid's connectivity is known
        # in closed form; computing it geometrically here would only add
        # floating-point risk for zero benefit.
        nbrs = String[]
        r > 0        && push!(nbrs, _cell_id(r-1, c))   # south
        r < nrows-1  && push!(nbrs, _cell_id(r+1, c))   # north
        c > 0        && push!(nbrs, _cell_id(r, c-1))   # west
        c < ncols-1  && push!(nbrs, _cell_id(r, c+1))   # east
        adjacency[id] = nbrs
    end

    mesh = A5Mesh(0, "", cells,
                  Dict{String,Vector{Float64}}("cell_size_m" => cell_size_m),
                  Dict{String,Matrix{Float64}}(),
                  adjacency, "cartesian")   # own row/column IDs, not pya5's
    return mesh
end
